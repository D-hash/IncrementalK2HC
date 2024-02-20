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
#include <queue>
#include <set>
#include <limits>
#include "progressBar.h"
#include "networkit/graph/Graph.hpp"
#include <string>
#include "mytimer.h"

using vertex = uint32_t;
using dist = uint16_t;
using edge_id = uint32_t;

class IncrementalTopK{
        struct index_t{
            std::vector<std::pair<vertex,dist>> label_offset;
            std::vector<std::vector<dist>> d_array;
            index_t(){
                label_offset.clear();
                d_array.clear();
            }
            ~index_t(){
                label_offset.clear();
                d_array.clear();
            }
        };

public:
    IncrementalTopK(NetworKit::Graph*, dist, bool, dist, bool);

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


    uint64_t compute_index_size();

    double n_reached_nodes();
    double n_reached_nodes_mbfs();
    // void mod_bfs(vertex, vertex, std::vector<dist>&);


    double loops_time;
    double lengths_time;
    vertex loop_entries;
    vertex length_entries;
    vertex aff_hubs;
    vertex aff_cycles;

    void deallocate_aux();

private:
    std::vector<dist> dists;
    std::vector<vertex> updated;
    std::queue<std::tuple<vertex, vertex, dist, dist, bool, vertex>> new_labels;
    std::set<vertex> vertices_to_update;
    bool is_from_scratch_only;
    std::vector<dist> tmp_v;
    void verify_sizes();
    void pruned_bfs (vertex, bool);
    void reset_temp_vars(vertex, bool);
    void set_temp_vars(vertex, bool);
    bool prune(vertex,  dist, bool);
    void compute_loop_entries(vertex);
    void resume_pbfs(vertex, vertex, dist, bool);
    void allocate_label(vertex, vertex, dist, dist, bool);
    void extend_label(vertex, vertex, dist, dist, bool, size_t);
    void extend_label_repair(vertex, vertex, dist, dist, bool);

    static const dist null_distance;
    static const vertex null_vertex;
    std::queue<vertex> * node_que;

    vertex* ordering;
    vertex* reverse_ordering; 
    std::pair<double,vertex>* ordering_rank;
    vertex total;
    NetworKit::Graph * graph;
    dist K;
    bool directed;
    dist ordering_type;

    std::vector<std::pair<vertex,dist>> old_label_a;
    std::vector<std::pair<vertex,dist>> old_label_b;
    std::vector<std::vector<dist>> old_distances_a;
    std::vector<std::vector<dist>> old_distances_b;


    std::vector<vertex> reached_nodes;
    std::vector<vertex> reached_mbfs;


    std::vector<dist>* loop_labels;
    index_t ** length_labels;


    bool* tmp_pruned;
    dist* tmp_offset;
    vertex* tmp_count;
    std::vector<dist>*  tmp_dist_count;
    dist*  tmp_s_offset;
    std::vector<dist>* tmp_s_count;
    dist* visited_in_update_loops;

};
#endif //INCREMENTAL_TOPK_H
