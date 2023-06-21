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

class IncrementalTopK{

        struct index_t{
            std::vector<std::pair<uint32_t, uint8_t>>   label_offset;
            std::vector<std::vector<uint8_t>>           d_array;
        };

public:
    IncrementalTopK() :
        V(0), K(0), directed(0), loop_count_time(0), indexing_time(0){
            index[0].clear();
            index[1].clear();
        }
    ~IncrementalTopK();
    int KDistanceQuery(int s, int t, uint8_t k, std::vector<int> &ret);
    int KDistanceQuery(int s, int t, std::vector<int> &ret){ return KDistanceQuery(s, t, K, ret); }
    int KDistanceQuery(int s, int t, uint8_t k);
    int KDistanceQuery(int s, int t){ return KDistanceQuery(s, t, K); }
    bool   ConstructIndex(NetworKit::Graph* graph, size_t K, bool directed);
    void   AddEdge(uint32_t a, uint32_t b);
    void   RemoveEdge(uint32_t a, uint32_t b);
    void   UpdateLoops(std::pair<int, int>);
    void   UpdateIndex(std::pair<int, int>);

    double IndexingTime()  const { return indexing_time; }
    double LoopCountTime() const { return loop_count_time; }
    uint32_t AffectedHubs() const { return aff_hubs; }
    double ReachedNodes() const { uint64_t sum= 0; for(auto n: reached_nodes) sum+=n; return sum / reached_nodes.size();}
    uint32_t AffectedCycles() const { return aff_cycles; }
    double ReachedMBFS() const { uint64_t sum= 0; for(auto n: reached_mbfs) sum+=n; return reached_mbfs.size()>0 ? sum / reached_mbfs.size() : 0;}
    void modBFS(uint32_t s, uint32_t t, std::vector<int> &ret);
    size_t NumOfVertex ();
    size_t IndexSize();
    double AvgIndexSize();

    std::vector<uint32_t> ordering;
    std::vector<uint32_t> reverse_ordering;
    NetworKit::Graph* graph;

    uint32_t loopcounter;
    uint32_t labelscounter;
    std::vector<std::vector<uint8_t> > loop_count;
    std::vector<index_t> index[2];

private:
    size_t V;
    uint8_t K;
    bool directed;
    // We assume that the diameter of a given network is less than 128.
    static const uint8_t INF8;

    double loop_count_time;
    double indexing_time;
    uint32_t aff_hubs;
    std::vector<uint32_t> reached_nodes;
    uint32_t aff_cycles;
    std::vector<uint32_t> reached_mbfs;


    // index[0] corresponds to L_in. index[1] corresponds to L_out

    std::vector<bool>     tmp_pruned;
    std::vector<uint32_t> tmp_offset;
    std::vector<uint32_t> tmp_count;
    std::vector<uint8_t>  tmp_dist_count[2];
    std::vector<uint8_t>  tmp_s_offset;
    std::vector<std::vector<uint8_t> > tmp_s_count;
    std::vector<uint8_t> visited_in_update_loops;

    void Init();
    void Free();
    void FreeAuxiliary();
    bool Labeling();
    void CountLoops(uint32_t s, bool &status);
    void PrunedBfs (uint32_t s, bool dir, bool &status);
    void ResumePBfs (uint32_t s, uint32_t t, uint8_t  d, bool dir, bool &status,
                     std::vector<std::tuple<u_int32_t, u_int32_t, uint8_t, u_int8_t, bool, u_int32_t>> &new_labels);
    inline void SetStartTempVars(uint32_t s, bool rev);
    inline void ResetTempVars(uint32_t s, const std::vector<uint32_t> &updated, bool rev);
    inline void AllocLabel(uint32_t v, uint32_t s, uint8_t d, uint8_t dc, bool rev);
    inline void ExtendLabel(uint32_t v, uint32_t s, uint8_t d, uint8_t dc, bool rev, size_t pos);
    inline void  ExtendLabelRepair(uint32_t v, uint32_t start, uint8_t dist, uint8_t count, bool dir);
    bool Pruning(uint32_t v,  uint8_t d, bool rev);

    void KHubsDistanceQuery(uint32_t s, uint32_t t, uint8_t k, std::vector<std::pair<u_int32_t,int>> &ret);
};
#endif //INCREMENTAL_TOPK_H
