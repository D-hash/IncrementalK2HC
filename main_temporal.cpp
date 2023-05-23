// Created by andrea
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cassert>
#include <map>
#include <algorithm>
#include "incremental_topk.h"
#include <string>
#include <bits/stdc++.h>
#include "networkit/io/EdgeListReader.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/centrality/Closeness.hpp"
#include "networkit/graph/GraphTools.hpp"

using namespace std;

double median(vector<double> &a) {
    int n = a.size();
    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(),
                    a.begin() + n / 2,
                    a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(),
                    a.begin() + (n - 1) / 2,
                    a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double) (a[(n - 1) / 2]
                         + a[n / 2])
               / 2.0;
    }

        // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(),
                    a.begin() + n / 2,
                    a.end());

        // Value at index (N/2)th
        // is the median
        return (double) a[n / 2];
    }
}


double GetCurrentTimeInSec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

double average(std::vector<double> const& v){
    if(v.empty()){
        return 0;
    }

    auto const count = static_cast<double>(v.size());
    double sum = 0;
    for(double time: v) sum += time;
    return sum / count;
}

uint32_t read_graph(const string graph_file, vector<pair<uint64_t, pair<uint32_t , uint32_t>> > &es){
    uint32_t max_id = 0;
    ifstream ifs(graph_file);
    if (!ifs.good()){
        cerr << "Error: Cannot open " << graph_file << "." << endl;
        exit(EXIT_FAILURE);
    }
    uint32_t vertices, edges, weighted, directed;
    ifs >> vertices >> edges >> weighted >> directed;
    es.clear();
    uint32_t  u, v, d;
    for (uint64_t t; ifs >> t >> u >> v >> d;){
        es.emplace_back(t, make_pair(u, v));
        max_id = max(max_id, u);
        max_id = max(max_id, v);
    }
    //assert(es.size() == edges);
    //assert(vertices == max_id + 1);

    ifs.close();
    return max_id;
}

int main(int argc, char **argv) {
    srand (time(NULL));
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " (graph_file) (K) (D) (I)" << endl;
        exit(EXIT_FAILURE);
    }

    string  graph_file = argv[1];
    const size_t  K          = atoi(argv[2]);
    const bool    directed   = atoi(argv[3]) != 0;
    long long int    input_ins   = atoi(argv[4]);
    vector<pair<uint64_t, pair<uint32_t , uint32_t>> > edges;
    uint32_t V = read_graph(graph_file, edges);
    std:cout << "Starting on " << graph_file << " with K= " << K << "\n";
    std::cout << "Sorting temporal edges...\n";
    sort(edges.begin(), edges.end());
    std::cout << "Sorting done!\n";
    // NetworKit::EdgeListReader* elr = new NetworKit::EdgeListReader(' ', 0);
    //NetworKit::Graph raw_g = elr->read(graph_file);
    NetworKit::Graph raw_g;
    raw_g.addNodes(V+1);
    for(auto e: edges)
        raw_g.addEdge(e.second.first, e.second.second);
    raw_g.removeSelfLoops();
    raw_g.removeMultiEdges();
    uint32_t  initial_number_of_nodes = raw_g.numberOfNodes();
    std::cout << "Graph with " << raw_g.numberOfNodes() << " vertices and " << raw_g.numberOfEdges() << " edges.\n";
    long long int num_insertions = std::min((long long int)(input_ins), (long long int)(raw_g.numberOfNodes()*(raw_g.numberOfNodes()-1)/2 - raw_g.numberOfEdges()));
    std::cout << "Number of insertions " << num_insertions << "\n";
    std::cout << "Removing " << num_insertions << " edges\n";
    vector<pair<uint32_t, uint32_t> > removed_edges;
    long long int ni = 0;
    for(size_t i = edges.size()-1; ni < num_insertions; i--){
        if(find(removed_edges.begin(), removed_edges.end(), make_pair(edges[i].second.first, edges[i].second.second)) != removed_edges.end() ||
                find(removed_edges.begin(), removed_edges.end(), make_pair(edges[i].second.second, edges[i].second.first)) != removed_edges.end())
            continue;
        if(edges[i].second.first == edges[i].second.second) continue;
        removed_edges.emplace_back(edges[i].second.first, edges[i].second.second);
        raw_g.removeEdge(edges[i].second.first, edges[i].second.second);
        ni++;
    }
    edges.clear();
    std::cout << "Edges after removal " << raw_g.numberOfEdges() << '\n';
    IncrementalTopK* kpll = new IncrementalTopK();
    kpll->ConstructIndex(raw_g, K, directed);
    raw_g.~Graph();
    std::cout << "First Labeling Loop time: " << kpll->LoopCountTime() << "s | First Labeling Indexing time:" << kpll->IndexingTime()
              << "\n";
    std::cout << "Number Vertices: " << kpll->NumOfVertex() << "\n";

    std::ofstream ofs;
    ofs.open(graph_file+"_"+std::to_string(K)+"_"+std::to_string(num_insertions)+"_temporal.csv");
    ofs << "Graph,Vertices,Edges,K,Insertions,"
           "ULLoopTime,ULLabelingTime,ULSize,ULMeanQueryTime,"
           "ULMedianQueryTime,AffectedHubs,ReachedNodes\n";
    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << 0 << ","
        << kpll->LoopCountTime() << "," << kpll->IndexingTime() << "," << kpll->IndexSize() << ","
        << 0 << ","
        << 0 << "," << 0 << "," << 0 <<"\n";
    int num_queries = 100000;

    std::vector<double> update_loops;
    std::vector<double> update_lengths;
    std::vector<double> avg_index_size;
    std::vector<size_t> index_size;
    std::vector<uint32_t> affected_hubs;
    std::vector<double> reached_nodes;
    int i = 0;
    for(auto edge: removed_edges){
        uint32_t a = edge.first;
        uint32_t b = edge.second;
        std::cout << "New edge " << a << " " << b << "\n";
        std::cout << "Updating first labeling.. \n";
        double ul_loops = -GetCurrentTimeInSec();
        kpll->UpdateLoops(make_pair(a,b));
        ul_loops += GetCurrentTimeInSec();
        double ul_labeling = -GetCurrentTimeInSec();
        kpll->UpdateIndex(make_pair(a,b));
        ul_labeling += GetCurrentTimeInSec();
        std::cout << "Update index time: " << ul_labeling << " | Update loops time: " << ul_loops << "\n";
        std::cout << i+1 << "-th insertion correct!" << "\n";
        update_loops.push_back(ul_loops);
        update_lengths.push_back(ul_labeling);
        // avg_index_size.push_back(kpll->AvgIndexSize());
        index_size.push_back(0);
        affected_hubs.push_back(kpll->AffectedHubs());
        reached_nodes.push_back(kpll->ReachedNodes());
        i++;
    }

    IncrementalTopK scratch_kpll;
    scratch_kpll.ConstructIndex(kpll->graph, K, directed);
    std::cout << "Scratch LB Loop time: " << scratch_kpll.LoopCountTime() << "s | Scratch LB Indexing time:"
              << scratch_kpll.IndexingTime()
              << "\n";
    vector<double> sl_time;
    vector<double> khl_time;
    ProgressStream query_bar(num_queries);
    query_bar.label() << "Queries";
    for(int j=0; j<num_queries; j++){
        int32_t u = NetworKit::GraphTools::randomNode(kpll->graph);
        int32_t v = NetworKit::GraphTools::randomNode(kpll->graph);
        vector<int> up_dist;
        vector<int> sc_dist;
        double khl_query_time = -GetCurrentTimeInSec();
        kpll->KDistanceQuery(u, v, up_dist);
        khl_query_time += GetCurrentTimeInSec();
        khl_time.push_back(khl_query_time);
        double scratch = -GetCurrentTimeInSec();
        scratch_kpll.KDistanceQuery(u, v, sc_dist);
        scratch += GetCurrentTimeInSec();
        sl_time.push_back(scratch);
        assert(up_dist.size() == sc_dist.size());
        for(size_t l=0; l < up_dist.size(); l++){
            if(up_dist[l] != sc_dist[l]){
                std::cout << "Error bw " << u << "-" << v << "\n";
                std::cout << "Updated labeling distance: " << up_dist[l] << "\n";
                std::cout << "Scratch labeling distance: " << sc_dist[l] << "\n";
                for(size_t id=0; id < up_dist.size(); id++){
                    std:: cout << "Up " << up_dist[id] << " | Scratch " << sc_dist[id] << "\n";
                }
                assert(false);
            }
        }
        ++query_bar;
    }
    std::cout << "Writing on csv file\n";
    for(size_t j = 0; j < num_insertions; j++) {
        ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << ","
            << j + 1 << ","
            << update_loops[j] << "," << update_lengths[j] << "," << index_size[j] << ","
            << 0 << ","
            << 0 << "," << affected_hubs[j] << "," << reached_nodes[j] << "\n";
    }
    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << i+1 << ","
        << "final" << "," << "final" << "," << kpll->IndexSize() << ","
        << average(khl_time) << ","
        << median(khl_time) << "," << kpll->AffectedHubs() << "," << kpll->ReachedNodes() <<"\n";
    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << num_insertions << ","
        << scratch_kpll.LoopCountTime() << "," << scratch_kpll.IndexingTime() << "," << scratch_kpll.IndexSize() << ","
        << average(sl_time) << ","
        << median(sl_time) << ",scratch,scratch\n";
    std::cout << "Writing done!\n";
    ofs.close();
    delete kpll;

    exit(EXIT_SUCCESS);
}
