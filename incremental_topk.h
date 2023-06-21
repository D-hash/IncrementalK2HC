//
// Created by andrea on 20/12/22.
//#

#ifndef INCREMENTAL_TOPK_H
#define INCREMENTAL_TOPK_H
#include <vector>
#include <tuple>
#include <sys/time.h>
#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include "progressBar.h"
#include "networkit/graph/Graph.hpp"
#include <string>
#include "mytimer.h"

using vertex = int32_t;
using dist = uint32_t;
using edge_id = uint64_t;

class IncrementalTopK{
        struct index_t{
            std::vector<std::pair<vertex,dist>> label_offset;
            std::vector<std::vector<dist>> d_array;
        };

public:
    IncrementalTopK(NetworKit::Graph*, dist, bool, dist);

    ~IncrementalTopK();

    void build();

    // int query(vertex, vertex, dist, std::vector<dist> &);
    void query(vertex, vertex, std::vector<dist> &);
    // int query(vertex, vertex, dist);
    // int query(vertex, vertex);
    

    vertex x;
    vertex y;
    void update_loops();
    void update_lengths();




    double n_reached_nodes();
    // double n_reached_nodes_mbfs();
    void mod_bfs(vertex, vertex, std::vector<dist>&);


    double loops_time;
    double lengths_time;
    vertex loop_entries;
    vertex length_entries;
    vertex aff_hubs;
    vertex aff_cycles;


private:
    void verify_sizes();
    void pruned_bfs (vertex, bool);
    void reset_temp_vars(vertex, const std::vector<vertex> &, bool);
    void set_temp_vars(vertex, bool);
    bool prune(vertex,  dist, bool);
    void compute_loop_entries(vertex);
    void resume_pbfs(vertex, vertex, dist, bool, std::vector<std::tuple<vertex, vertex, dist, dist, bool, vertex>> &);
    void allocate_label(vertex, vertex, dist, dist, bool);
    void extend_label(vertex, vertex, dist, dist, bool, size_t);
    void extend_label_repair(vertex, vertex, dist, dist, bool);

    static const dist null_distance;
    static const vertex null_vertex;

    std::vector<vertex> ordering;
    std::vector<vertex> reverse_ordering;
    std::vector<std::pair<double,vertex> > ordering_rank;

    NetworKit::Graph * graph;
    dist K;
    bool directed;
    dist ordering_type;



    std::vector<vertex> reached_nodes;
    std::vector<vertex> reached_mbfs;


    std::vector<std::vector<dist> > loop_labels;
    std::vector<index_t> length_labels[2];


    std::vector<bool> tmp_pruned;
    std::vector<vertex> tmp_offset;
    std::vector<vertex> tmp_count;
    std::vector<dist>  tmp_dist_count[2];
    std::vector<dist>  tmp_s_offset;
    std::vector<std::vector<dist> > tmp_s_count;
    std::vector<dist> visited_in_update_loops;



    void deallocate_aux();


    

};
#endif //INCREMENTAL_TOPK_H
