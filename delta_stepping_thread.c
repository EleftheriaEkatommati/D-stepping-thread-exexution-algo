#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <pthread.h>

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
	Graph* g;              // Ο γράφος
	double* dist;          // Ο πίνακας αποστάσεων
	double delta;          // Η τιμή Δ
	int num_threads;       // Πλήθος threads

	// Buckets και Λίστες (Κοινόχρηστα)
	Bucket* B;           // Δείκτης στον πίνακα κάδων της main
	int* Bsize;            // Δείκτης στο μέγεθος του B
	Bucket* R;             // Η τρέχουσα λίστα R
	Bucket* S;             // Η τρέχουσα λίστα S

	// Εργαλεία Συγχρονισμού
	pthread_mutex_t lock_dist;   // Για τον πίνακα αποστάσεων
	pthread_mutex_t lock_buckets; // Για bucket_add σε R, S και B[i]
	pthread_barrier_t barrier;    // Για τον συντονισμό των φάσεων

	int current_bucket_idx; // Ποιον κάδο επεξεργαζόμαστε τώρα
	int terminate;          // Σημαία για το πότε τελείωσαν όλα
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

/* =========================
   Delta-stepping (correct)
   ========================= */

   /*
	* Computes single-source shortest paths
	* from vertex "source" using �-stepping.
	*/
   //DELTA is the parameter that determines the "light" and "heavy" edges. Edges with weight <= DELTA are considered light, while edges with weight > DELTA are considered heavy.
void *delta_stepping_thread(Graph* g, int source, double DELTA) {
	int n = g->n;

	/*
	 * Process buckets in increasing order.
	 * This loop is inherently sequential.
	 */
	for (int i = 0; i < Bsize; i++) {
		printf("----------------------------------------------------------------------------------\n");
		if (B[i]->size == 0) continue;
		printf("B[%d] size: %d\n", i, B[i]->size); 
		/*
		 * S = vertices settled for bucket i
		 * (snapshot of the bucket contents)
		 */
		Bucket* S = bucket_create(); 
		//here we create a new bucket S to hold the vertices that are currently being processed from bucket B[i]. This allows us to work with a stable set of vertices while we relax edges, even if new vertices are added to B[i] during the process.
		//if we had the source node the
		for (int k = 0; k < B[i]->size; k++) {
			bucket_add(S, B[i]->v[k].id);
		}
		
		/*
		 * R = worklist used for repeated
		 * relaxation of light edges
		 */
		Bucket* R = bucket_create();
		for (int k = 0; k < S->size; k++) {
			bucket_add(R, S->v[k].id);
		}
		//R is holding vertices that are being processed for the current bucket.
		//B is holding the nodes that are in the current bucket 
		//S is holding the nodes that are being processed for the current bucket  and it is filled with the nodes thatt are in the current bucket B, for source node  
		// it will pass the light edge stage and not the heavy edge stage, for the other nodes it will pass both stages,
		// then the bucket S will be empty but on the second stage the bucket S will be filled with the nodes that are in the current bucket B.


		/* Bucket i is now logically empty */
		B[i]->size = 0;

		/* -------- light edges (fixpoint) --------
		 *
		 * MAIN PARALLELIZATION TARGET:
		 *   - Processing vertices from R
		 *   - Relaxing outgoing light edges
		 */


		//Protect shared variables between nodes(Dist[],bucket list,R)  -Ria
		//Each thread will take a part of R 
	
		while (R->size > 0) {
			printf(
				"in the loop ------------------- \n"
			);


			for (int i=0; i<R->size; i++) {
				printf("R->v[%d].id: %d\n", i, R->v[i].id);
			}

			
			printf(" R size: %d\n", R->size);

			// for (int j = 0; j < Bsize; j++) {
			//     printf("B[%d] size: %d\n", j, B[j]->size);
			// }
			
			int u = R->v[--R->size].id;


			 // Debug print
			printf("R->v[0].id: %d\n", R->v[0].id); // Debug print
			printf("Processing vertex %d from bucket %d\n", u, i); // Debug print
		   


			for (int e = 0; e < g->adj_size[u]; e++) {
				Edge edge = g->adj[u][e];
				
				if (edge.weight <= DELTA) {
					

					printf("Relaxing light edge (%d -> %d) with weight %.2f\n",
						u, edge.to, edge.weight); // Debug print


					double nd = dist[u] + edge.weight;  
					if (nd < dist[edge.to]) {
						dist[edge.to] = nd;         //lock -Ria
						int bidx = (int)floor(nd / DELTA);  

						printf("""Relaxing edge (%d -> %d) with weight %.2f: new distance = %.2f, bucket index = %d\n",
							u, edge.to, edge.weight, nd, bidx); // Debug print

							
						//lock  -Ria
						ensure_bucket_array(&B, &Bsize, bidx); //this ensures that the bucket array is large enough to hold the new bucket index

						//lock  -Ria
						if (bidx == i) {
							bucket_add(R, edge.to);    
						} else {
							bucket_add(B[bidx], edge.to); // here the B bucket  is
						}
						
					}
				}
			}
		}


	   
		/* -------- heavy edges (once) --------
		 *
		 * Secondary parallelization opportunity.
		 * Heavy edges only generate future work.
		 */
		for (int k = 0; k < S->size; k++) {
			int u = S->v[k].id;
			for (int e = 0; e < g->adj_size[u]; e++) {
				Edge edge = g->adj[u][e];
				if (edge.weight > DELTA) {
					double nd = dist[u] + edge.weight;
					if (nd < dist[edge.to]) {
						dist[edge.to] = nd;
						int bidx = (int)floor(nd / DELTA);
						ensure_bucket_array(&B, &Bsize, bidx);
						bucket_add(B[bidx], edge.to);
					}
				}
			}
		}

		free(R->v); free(R);
		free(S->v); free(S);
	}

	/* Output final distances */
	for (int v = 0; v < n; v++) {
		if (dist[v] == DBL_MAX) {
			printf("%d: INF\n", v);
		} else {
			printf("%d: %.6f\n", v, dist[v]);
		}
	}


}


int main(int argc, char** argv) {
	if (argc != 4) {
		fprintf(stderr, "Usage: %s graph_file source DELTA\n", argv[0]);
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
	//---------------------------------------------------
	int Bsize = 0; // size of the bucket
	Bucket** B = NULL; // each B[I] contains the adress of a bucket
	ensure_bucket_array(&B, &Bsize, 0); //ensures that we are not going to write to somehwere where there has not been created space
	bucket_add(B[0], source);

	// Preparing the shared structure across the threads
	SharedData sd;
	sd.g = g;
	sd.dist = dist;
	sd.delta = delta;
	sd.num_threads = num_threads;
	
	//We pass the pointers  in the structures in order in the future to do realloc
	sd.B = B; 
	sd.Bsize = &Bsize;
	sd.R = bucket_create();
	sd.S = bucket_create();
	sd.current_bucket_idx = 0; // counter in order to know in which bucket we are 
	
	//Initialize mutexes/barriers
	pthread_mutex_init(&sd.lock_dist, NULL);
	pthread_mutex_init(&sd.lock_buckets, NULL);
	pthread_barrier_init(&sd.barrier, NULL, num_threads);
	

	//Initialize the threads
	pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
	ThreadArgs* t_args = malloc(num_threads * sizeof(ThreadArgs));

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
	return 0;
}
