# IncrementalK2HC
This is the implementation of **DYN-KPLL** algorithm presented in "Top-k Distance Queries on Large Time-Evolving Graphs"
submitted at the *21st Symposium of Experimental Algorithms* (SEA2023).
Given a k-2-Hop Cover index \[Akiba et al., AAAI2015\] of a graph, the algorithm is able to update the index after 
incremental graph changes.

## Usage
The executable takes 4 arguments:
1. The input graph, specified by the path to an edge-list format file;
2. K, the parameter of the *k*-Shortest Distances problem;
3. A bool value, 0 if the input graph is undirected, 1 if it is directed (**warning**: the current implementation works only with undirected ones);
4. The number of edges to be randomly inserted.
