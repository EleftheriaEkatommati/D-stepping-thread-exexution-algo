import random

# =========================
# USER-ADJUSTABLE PARAMETERS
# =========================

NUM_VERTICES = 1000      # Number of vertices (N)
NUM_EDGES = 5000         # Number of edges (M)

DELTA = 1.0              # Delta parameter for delta-stepping
P_LIGHT = 0.7            # Probability an edge is "light"

LIGHT_WEIGHT_RANGE = (0.1, DELTA)
HEAVY_WEIGHT_RANGE = (DELTA + 0.1, 10.0)

SEED = 42                # Set to None for non-reproducible graphs

# =========================


def make_output_filename():
    """
    Construct the output filename using the agreed naming scheme:

      graph_NUM_VERTICES_NUM_EDGES_DELTA_P_LIGHT.txt

    Formatting rules:
      - NUM_VERTICES, NUM_EDGES: 6 digits, zero-padded
      - DELTA: 2-digit integer part + exactly 2 decimals
      - P_LIGHT: exactly 2 decimals
    """
    delta_str = f"{DELTA:05.2f}"
    p_light_str = f"{P_LIGHT:.2f}"

    return (
        f"graph_"
        f"{NUM_VERTICES:06d}_"
        f"{NUM_EDGES:06d}_"
        f"{delta_str}_"
        f"{p_light_str}.txt"
    )


def generate_graph():
    if SEED is not None:
        random.seed(SEED)

    # ------------------------------------------------------------------
    # We use TWO data structures:
    #
    # 1) edge_set : a set of (u, v) pairs
    #    - Used ONLY to test whether an edge between u and v already exists
    #    - Allows O(1) duplicate detection
    #
    # 2) edges : a list of (u, v, w) triples
    #    - Stores the actual edges including weights
    #    - Preserves all information needed for output
    #
    # This separation avoids O(M) scans when generating large graphs.
    # ------------------------------------------------------------------

    edge_set = set()   # stores (u, v) pairs
    edges = []         # stores (u, v, w) triples

    # --------------------------------------------------------------
    # PHASE 1: Ensure every vertex appears at least once as a source
    # --------------------------------------------------------------
    for u in range(NUM_VERTICES):
        v = random.randint(0, NUM_VERTICES - 1)
        while v == u:  # avoid self-loop
            v = random.randint(0, NUM_VERTICES - 1)

        edge_set.add((u, v))

        # Choose weight type
        if random.random() < P_LIGHT:
            w = random.uniform(*LIGHT_WEIGHT_RANGE)
        else:
            w = random.uniform(*HEAVY_WEIGHT_RANGE)

        edges.append((u, v, w))

    # --------------------------------------------------------------
    # PHASE 2: Add remaining edges randomly
    # --------------------------------------------------------------
    while len(edges) < NUM_EDGES:
        u = random.randint(0, NUM_VERTICES - 1)
        v = random.randint(0, NUM_VERTICES - 1)

        if u == v or (u, v) in edge_set:
            continue  # skip self-loops and duplicates

        edge_set.add((u, v))

        if random.random() < P_LIGHT:
            w = random.uniform(*LIGHT_WEIGHT_RANGE)
        else:
            w = random.uniform(*HEAVY_WEIGHT_RANGE)

        edges.append((u, v, w))

    # --------------------------------------------------------------
    # Optional: Sort edges by (u, v) for readability and determinism
    # --------------------------------------------------------------
    edges.sort(key=lambda e: (e[0], e[1]))

    # --------------------------------------------------------------
    # Write edges to file
    # --------------------------------------------------------------
    output_file = make_output_filename()

    with open(output_file, "w") as f:
        f.write(f"{NUM_VERTICES} {NUM_EDGES}\n")
        for u, v, w in edges:
            f.write(f"{u} {v} {w:.6f}\n")

    # --------------------------------------------------------------
    # Print summary
    # --------------------------------------------------------------
    print(f"Graph written to {output_file}")
    print(f"Vertices: {NUM_VERTICES}")
    print(f"Edges: {NUM_EDGES}")
    print(f"Delta: {DELTA}")
    print(f"Light edge fraction: {P_LIGHT}")


if __name__ == "__main__":
    generate_graph()
