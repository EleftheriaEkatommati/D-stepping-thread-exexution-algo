#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <sys/time.h>
/*
Commands to run : gcc -O2 -Wall -pthread delta_stepping_thread.c -o delta
./delta graph_000006_000005_01.00_0.70.txt 0 1 4
*/
/*
 * ============================================================
 *  SERIAL �-STEPPING SINGLE-SOURCE SHORTEST PATHS
 * ============================================================
 *
 * This code is a serial implementation of the
 * �-stepping algorithm (Meyer & Sanders).
 *
 * It is intended to be used as a BASELINE for a POSIX
 * threads parallelization assignment.
 */

#define INITIAL_ADJ_CAPACITY 4
#define INITIAL_BUCKET_CAPACITY 4
struct timeval start, end;
 /* =========================
	Data structures
	========================= */

	/* Directed weighted edge */
typedef struct {
	int to;
	double weight;
} Edge;

/* Simple vertex wrapper (used in buckets) */
typedef struct {
	int id;
} Vertex;

/*
 * Bucket:
 *  - holds vertices whose tentative distances
 *    fall into the same � interval
 */
typedef struct {
	Vertex* v;
	int size;
	int capacity;
} Bucket;

/*
 * Graph representation:
 *  - adjacency lists
 *  - directed edges
 */
typedef struct {
	int n, m;
	Edge** adj;
	int* adj_size; // number of edges in each adjacency list
	int* adj_cap; // capacity of each adjacency list (for dynamic resizing)
} Graph;

typedef struct {
	Graph* g;              // The graph
	double* dist;          // Distance array
	double delta;          // The Δ value
	int num_threads;       // Number of threads

	// Buckets and lists (shared)
	Bucket** B;           // Pointer to main bucket array
	int* Bsize;            // Pointer to size of B
	Bucket* R;             // Current R list
	Bucket* S;             // Current S list

	// Synchronization tools
	pthread_mutex_t lock_dist;   // For the distance array
	pthread_mutex_t lock_buckets; // For bucket_add on R, S and B[i]
	pthread_barrier_t barrier;    // For coordinating phases
    int terminate_inner;
	int current_bucket_idx; // Which bucket we are processing now

	Bucket*** all_local_B;     // [num_threads][bucket_index]
	int* all_local_Bsize;      // size for each thread
}  SharedData ;

typedef struct {
	SharedData* sd;
	int tid;
} ThreadArgs;


//Pthread variables
pthread_mutex_t dist_lock;   




/* =========================
   Graph handling
   ========================= */

   /*
	* Reads a directed graph from a file of the form:
	*
	*   n m
	*   u v w
	*   u v w
	*   ...
	*
	* where:
	* - n : number of vertices in graph
	* - m : number of edges in graph
	* - u v w : An edge from vertice 'u' to vertice 'v' with weight 'w'
	*/
Graph* read_graph(const char* filename) {
	FILE* f = fopen(filename, "r");
	if (!f) { perror("fopen"); exit(EXIT_FAILURE); }

	Graph* g = malloc(sizeof(Graph));
	fscanf(f, "%d %d", &g->n, &g->m);

	g->adj = malloc(g->n * sizeof(Edge*));
	g->adj_size = calloc(g->n, sizeof(int));
	g->adj_cap = malloc(g->n * sizeof(int));

	/* Initialize adjacency lists */
	for (int i = 0; i < g->n; i++) {
		g->adj_cap[i] = INITIAL_ADJ_CAPACITY;
		g->adj[i] = malloc(INITIAL_ADJ_CAPACITY * sizeof(Edge));
	}

	/* Read edges */
	for (int i = 0; i < g->m; i++) {
		int u, v;
		double w;
		fscanf(f, "%d %d %lf", &u, &v, &w);

		/* Grow adjacency list if needed */
		if (g->adj_size[u] == g->adj_cap[u]) {
			g->adj_cap[u] *= 2;
			g->adj[u] = realloc(g->adj[u],
				g->adj_cap[u] * sizeof(Edge));
		}

		g->adj[u][g->adj_size[u]++] = (Edge){ v, w };
	}

	fclose(f);
	return g;
}

/* Frees all graph memory */
void free_graph(Graph* g) {
	for (int i = 0; i < g->n; i++) free(g->adj[i]);
	free(g->adj);
	free(g->adj_size);
	free(g->adj_cap);
	free(g);
}

/* =========================
   Bucket handling
   ========================= */

   /*
	* Creates an empty bucket.
	*/
Bucket* bucket_create(void) {
	Bucket* b = malloc(sizeof(Bucket));
	b->size = 0;
	b->capacity = INITIAL_BUCKET_CAPACITY;
	b->v = malloc(b->capacity * sizeof(Vertex));
	return b;
}

/*
 * Adds a vertex to a bucket.
 *
 * NOTE FOR PARALLEL VERSION:
 *   This function will require synchronization.
 */
void bucket_add(Bucket* b, int v) {
	if (b->size == b->capacity) {
		b->capacity *= 2;
		b->v = realloc(b->v, b->capacity * sizeof(Vertex)); 
	}
	b->v[b->size++].id = v;
}

/*
 * Ensures that the bucket array has at least
 * idx+1 buckets (created lazily).
 */
void ensure_bucket_array(Bucket*** B, int* Bsize, int idx) {
	if (idx >= *Bsize) {
		int new_size = idx + 1;
		*B = realloc(*B, new_size * sizeof(Bucket*));
		for (int i = *Bsize; i < new_size; i++) {
			(*B)[i] = bucket_create();
		}
		*Bsize = new_size;
	}
}

//Each thread will run this method
void * delta_stepping_thread(void* args){
	ThreadArgs* t_args = (ThreadArgs*)args;  // the packet of arguments that we pass to the thread, it contains the shared data and the thread id
	SharedData* sd = t_args->sd;  // the shared data that contains all the common resources that the threads will use
	int tid = t_args->tid; // the thread id, it is used to determine which part of the work each thread will do
	
	int i = 0;

	//Local Buckets
	Bucket** local_B = NULL;
    int local_Bsize = 0;
	sd->all_local_B[tid] = local_B;
    sd->all_local_Bsize[tid] = local_Bsize;

	while (1) {
		pthread_barrier_wait(&sd->barrier);
		// Check for termination: Is the current bucket empty AND no more buckets have work?
		// To keep it simple, we can have Thread 0 check if we've exceeded a reasonable range
		// or if a global "empty" flag is set.
		if (tid == 0) {
			// Find the next non-empty bucket
			while (i < *(sd->Bsize) && sd->B[i]->size == 0) {
				i++;
			}
			sd->current_bucket_idx = i; // Share the found index
		}
		pthread_barrier_wait(&sd->barrier);

		i = sd->current_bucket_idx; // Everyone moves to the next bucket that actually has work

		
		if (i >= *(sd->Bsize)) break; // Termination condition

		// --- PHASE 1: PREPARATION (Leader fills R and S from B[i]) ---
		//R vertices for light relaxation / S vertices for heavy relaxation
		if (tid == 0) {
			for (int j = 0; j < sd->B[i]->size; j++) {
				int v = sd->B[i]->v[j].id;
				bucket_add(sd->R, v);
				bucket_add(sd->S, v);
			}
			sd->B[i]->size = 0;
		}
		pthread_barrier_wait(&sd->barrier);
		 //all the threads wait until the first thread has finished adding.
		 /*
		 * here we will proceed to the light edge relaxation
		 * each thread will take a part of the vertices in R and relax their outgoing light edges
		 * we need to protect the shared variables (dist and buckets) when we relax the edges, since multiple threads can try to update the same vertex or add the same vertex to a bucket at the same time
		 * we can use mutexes to protect these shared variables
		 */
		//R size changes dynamically, We set while(1)->Loop while theres work to do/untill we break
		while (1) {
			int u = -1;
			pthread_mutex_lock(&sd->lock_buckets); // one thread at a time will take a vertex from R to process it
			if (sd->R->size > 0) {
				u = sd->R->v[--sd->R->size].id; 
			}
			pthread_mutex_unlock(&sd->lock_buckets);
			// Relax outgoing light edges of vertex u
            if (u != -1) {
                for (int e = 0; e < sd->g->adj_size[u]; e++) {
                    Edge edge = sd->g->adj[u][e];
                    if (edge.weight <= sd->delta) { 
                        double new_dist = sd->dist[u] + edge.weight;

                        pthread_mutex_lock(&sd->lock_dist); //first mutex for locking the dist array 
                        if (new_dist < sd->dist[edge.to]) { 
                            sd->dist[edge.to] = new_dist; // if the new distance is better, we update it and add the vertex to the appropriate bucket
                            pthread_mutex_unlock(&sd->lock_dist);

                            int bidx = (int)floor(new_dist / sd->delta);
                            
                            if (bidx == i) // if the new bucket index is the same as the curent bucket index 
                            {
                                pthread_mutex_lock(&sd->lock_buckets);
                                bucket_add(sd->R, edge.to); // we add the vertex to R since it will be processed in the current bucket
                                bucket_add(sd->S, edge.to);  // we also add the vertex to S since it will be processed in the current bucket and we need to keep track of it for the heavy edge relaxation
                                pthread_mutex_unlock(&sd->lock_buckets);
                            } 
                            else
                            {	
                                ensure_bucket_array(&local_B, &local_Bsize, bidx);
                                bucket_add(local_B[bidx], edge.to);
                                sd->all_local_B[tid] = local_B;
                                sd->all_local_Bsize[tid] = local_Bsize;
                            }

                        }else{
                            pthread_mutex_unlock(&sd->lock_dist);
                            // if the new distance is not better, we can skip the rest of the loop
                        }
                    }
                }
            }
            pthread_barrier_wait(&sd->barrier);
            if (tid == 0) {
                if (sd->R->size == 0) sd->terminate_inner = 1;
                else sd->terminate_inner = 0;
            }

            pthread_barrier_wait(&sd->barrier);

            if (sd->terminate_inner) break;
		}
		/*
		* After all the threads have finished the light edge relaxtion phase
		* we need to synchronize the threads before we proceed to the heavy edge relxation phase
		* this is because the heavy edge relaxation phase needs to work with the vertices in S,
		* which are shared among all threads and we need to make sure that all threads have finished updating S before we start processing it
		*/
		pthread_barrier_wait(&sd->barrier);

		//Merge the local buckets(Only thread 0)
		if (tid == 0) {

            for (int t = 0; t < sd->num_threads; t++) {

                Bucket** t_local = sd->all_local_B[t];
                int t_size = sd->all_local_Bsize[t];

                for (int j = 0; j < t_size; j++) {

                    if (t_local[j] != NULL && t_local[j]->size > 0) {

                        ensure_bucket_array(&(sd->B), sd->Bsize, j);

                        for (int k = 0; k < t_local[j]->size; k++)
                            bucket_add(sd->B[j], t_local[j]->v[k].id);

                        t_local[j]->size = 0;
                    }
                }
            }
        }
		pthread_barrier_wait(&sd->barrier);
		//S is divided into fixed-size chunks, one per thread 
        int total_s = sd->S->size; // 
        int chunk = (total_s + sd->num_threads - 1) / sd->num_threads; // calculating the number of nodes in s that each thread  will process 
        int start = tid * chunk; // based on the thread id, we calculate the start and the end 
        int end = start + chunk; // for the part of S that each thread will process
        if (end > total_s){  // if the end index is greater than the total size of S, we set it to the total size of S to avoid out of bounds access
			end = total_s; 
		}

        for (int k = start; k < end; k++) {
            int u = sd->S->v[k].id; // we process only the part of S that is assigned to the current thread
            for (int e = 0; e < sd->g->adj_size[u]; e++) { 
                Edge edge = sd->g->adj[u][e]; // we relax only the heavy edges, sd is the shared data among the threads , g is the graph  and adg is the adjacency list of the graph

                if (edge.weight > sd->delta) {
					
                    double new_dist = sd->dist[u] + edge.weight;

                    pthread_mutex_lock(&sd->lock_dist); // lock the dist array 
                    if (new_dist < sd->dist[edge.to]) {
                        sd->dist[edge.to] = new_dist;
                        pthread_mutex_unlock(&sd->lock_dist);

                        int bidx = (int)floor(new_dist / sd->delta);
						ensure_bucket_array(&local_B, &local_Bsize, bidx);
                        bucket_add(local_B[bidx], edge.to);
                        sd->all_local_B[tid] = local_B;
                        sd->all_local_Bsize[tid] = local_Bsize;
                    } else {
                        pthread_mutex_unlock(&sd->lock_dist);
                    }
                }
            }
		}
        sd->all_local_B[tid] = local_B;
        sd->all_local_Bsize[tid] = local_Bsize;

		pthread_barrier_wait(&sd->barrier);
		//Final merge by thread 0
		if (tid == 0) {

            for (int t = 0; t < sd->num_threads; t++) {

                Bucket** t_local = sd->all_local_B[t];
                int t_size = sd->all_local_Bsize[t];

                for (int j = 0; j < t_size; j++) {

                    if (t_local[j] != NULL && t_local[j]->size > 0) {

                        ensure_bucket_array(&(sd->B), sd->Bsize, j);

                        for (int k = 0; k < t_local[j]->size; k++)
                            bucket_add(sd->B[j], t_local[j]->v[k].id);

                        t_local[j]->size = 0;
                    }
                }
            }
        }
		//Clean R and S
		pthread_barrier_wait(&sd->barrier);
        if (tid == 0) {
            sd->R->size = 0;
            sd->S->size = 0;
        }
        pthread_barrier_wait(&sd->barrier);
		//Clear local buckets
		for (int j = 0; j < local_Bsize; j++) {
            if (local_B[j] != NULL)
                local_B[j]->size = 0;
        }
    }

	 for (int j = 0; j < local_Bsize; j++) {
        if (local_B[j] != NULL) {
            free(local_B[j]->v);
            free(local_B[j]);
        }
    }
    free(local_B);

	return NULL;
}

int main(int argc, char** argv) {
	if (argc != 5) {
		fprintf(stderr, "Usage: %s graph_file source DELTA num_threads\n", argv[0]);
		return EXIT_FAILURE;
	}

	Graph* g = read_graph(argv[1]);
	int source = atoi(argv[2]);
	double DELTA = atof(argv[3]);
	int num_threads=atoi(argv[4]);
	
	
	
	//We do all the initializations here instead of the delta_stepping function
	//since each thread needs to have access to the common shared resources
	//here we initialize the dist vector 
	double* dist = malloc(g->n * sizeof(double));
	for (int i = 0; i < g->n; i++)
	{
		dist[i] = DBL_MAX;
	}
	dist[source] = 0.0;

    // each B[I] contains the adress of a bucket
	// ensure_bucket_array(&B, &Bsize, 0); //ensures that we are not going to write to somehwere where there has not been created space
	int Bsize = 100000;
    Bucket** B = malloc(Bsize * sizeof(Bucket*));
    for (int j = 0; j < Bsize; j++) {
        B[j] = bucket_create();
    }
 
	bucket_add(B[0], source);

	// Preparing the shared structure across the threads
	SharedData sd;
	sd.g = g;
	sd.dist = dist;
	sd.delta = DELTA;
	sd.num_threads = num_threads;
	
	//We pass the pointers  in the structures in order in the future to do realloc
	sd.B = B; 
	sd.Bsize = &Bsize;
	sd.R = bucket_create();
	sd.S = bucket_create();
	sd.current_bucket_idx = 0; // counter in order to know in which bucket we are 
    sd.terminate_inner = 0;
	
	//Initialize mutexes/barriers
	pthread_mutex_init(&sd.lock_dist, NULL);
	pthread_mutex_init(&sd.lock_buckets, NULL);
	pthread_barrier_init(&sd.barrier, NULL, num_threads);
	
	sd.all_local_B = malloc(num_threads * sizeof(Bucket**));
	sd.all_local_Bsize = calloc(num_threads, sizeof(int));

	for (int t = 0; t < num_threads; t++) {
    	sd.all_local_B[t] = NULL;
	}	

	//Initialize the threads
	pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
	ThreadArgs* t_args = malloc(num_threads * sizeof(ThreadArgs));
	gettimeofday(&start, NULL);
	for (int i = 0; i < num_threads; i++) {
		t_args[i].sd = &sd;
		t_args[i].tid = i;

		pthread_create(&threads[i], NULL,
					   delta_stepping_thread,
					   &t_args[i]);
	}

	for (int i = 0; i < num_threads; i++) {
			pthread_join(threads[i], NULL);
	}
	gettimeofday(&end, NULL);

	double delta_time = (end.tv_sec - start.tv_sec) + 
						(end.tv_usec - start.tv_usec) / 1000000.0;

	printf("Parallel Delta-Stepping took: %.6f seconds\n", delta_time);
	
	//Print the results
	for (int v = 0; v < g->n; v++) {
		if (dist[v] == DBL_MAX) printf("%d: INF\n", v);
		else printf("%d: %.6f\n", v, dist[v]);
	}

	//Destroy the mutexes/barriers
	pthread_barrier_destroy(&sd.barrier);
	pthread_mutex_destroy(&sd.lock_buckets);
	pthread_mutex_destroy(&sd.lock_dist);

    //Free memory
	free_graph(g);
	free(dist);
	free(threads);
	free(t_args);
	free(sd.R->v);
    free(sd.R);
    free(sd.S->v);
    free(sd.S);
    for (int j = 0; j < Bsize; j++) {
        if (B[j] != NULL) {
            free(B[j]->v); 
            free(B[j]);    
        }
    }
    free(B);
	free(sd.all_local_B);
	free(sd.all_local_Bsize);
	return 0;
}