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

void read_graph(const string graph_file, vector<pair<uint32_t , uint32_t> > &es){
    ifstream ifs(graph_file);
    if (!ifs.good()){
        cerr << "Error: Cannot open " << graph_file << "." << endl;
        exit(EXIT_FAILURE);
    }

    es.clear();
    for (uint32_t u, v; ifs >> u >> v;){
        es.push_back(make_pair(u, v));
    }

    ifs.close();
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
    std:cout << "Starting on " << graph_file << " with K= " << K << "\n";
    NetworKit::EdgeListReader* elr = new NetworKit::EdgeListReader(' ', 0);
    NetworKit::Graph raw_g = elr->read(graph_file);
    raw_g.removeSelfLoops();
    raw_g.removeMultiEdges();
    NetworKit::ConnectedComponents *bic = new NetworKit::ConnectedComponents(raw_g);
    bic->run();
    auto g = bic->extractLargestConnectedComponent(raw_g,true);
    raw_g.~Graph();
    std::cout << "Graph with " << g.numberOfNodes() << " vertices and " << g.numberOfEdges() << " edges.\n";
    long long int num_insertions = std::min((long long int)(input_ins), (long long int)(g.numberOfNodes()*(g.numberOfNodes()-1)/2 - g.numberOfEdges()));
    std::cout << "Number of insertions " << num_insertions << "\n";
    //std::cout << "Removing " << num_insertions << " edges\n";
//    vector<pair<uint32_t, uint32_t> > removed_edges;
//    for(size_t i = 0; i < num_insertions; i++){
//        uint32_t a = NetworKit::GraphTools::randomNode(g);
//        uint32_t b = NetworKit::GraphTools::randomNeighbor(g,a);
//
//        g.removeEdge(a,b);
//        bic = new NetworKit::ConnectedComponents(g);
//        bic->run();
//        while(bic->numberOfComponents()>1){
//            g.addEdge(a,b);
//            a = NetworKit::GraphTools::randomNode(g);
//            b = NetworKit::GraphTools::randomNeighbor(g,a);
//            g.removeEdge(a,b);
//            bic->run();
//        }
//        removed_edges.emplace_back((uint32_t) a, (uint32_t) b);
//    }
//    std::cout << "Edges after removal " << g.numberOfEdges() << '\n';

    IncrementalTopK* kpll = new IncrementalTopK();
    kpll->ConstructIndex(g, K, directed);
    std::cout << "First Labeling Loop time: " << kpll->LoopCountTime() << "s | First Labeling Indexing time:" << kpll->IndexingTime()
              << "\n";
    std::cout << "Number Vertices: " << kpll->NumOfVertex() << "\n";

    g.~Graph();

    std::ofstream ofs;
    ofs.open(graph_file+"_"+std::to_string(K)+"_"+std::to_string(num_insertions)+"_random_only_dyn.csv");
    ofs << "Graph,Vertices,Edges,K,Insertions,NewEdgeX,NewEdgeY,"
           "ULLoopTime,ULLabelingTime,ULSize,ULMeanQueryTime,"
           "ULMedianQueryTime,AffectedHubs,ReachedNodes\n";
    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << 0 << ","
        << 0 << "," << 0 << ","  << kpll->LoopCountTime() << "," << kpll->IndexingTime() << "," << kpll->IndexSize() << ","
        << 0 << ","
        << 0 << "," << 0 << "," << 0 <<"\n";
    int num_queries = 100000;
    std::vector<double> update_loops;
    std::vector<double> update_lengths;
    std::vector<double> avg_index_size;
    std::vector<size_t> index_size;
    std::vector<uint32_t> affected_hubs;
    std::vector<double> reached_nodes;
    std::vector<pair<uint32_t, uint32_t>> added_edges;
    for(int i=0; i < num_insertions; i++){
        uint32_t a = NetworKit::GraphTools::randomNode(g);
        uint32_t b = NetworKit::GraphTools::randomNode(g);

        while(kpll->graph.hasEdge(a,b) || a == b){
            a = NetworKit::GraphTools::randomNode(g);
            b = NetworKit::GraphTools::randomNode(g);
        }
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
        index_size.push_back(0);
        affected_hubs.push_back(kpll->AffectedHubs());
        reached_nodes.push_back(kpll->ReachedNodes());
        added_edges.emplace_back(a,b);
        i++;
    }

    vector<double> khl_time;
    ProgressStream query_bar(num_queries);
    query_bar.label() << "Queries";
    for(int j=0; j<num_queries; j++){
        int32_t u = NetworKit::GraphTools::randomNode(kpll->graph);
        int32_t v = NetworKit::GraphTools::randomNode(kpll->graph);
        vector<int> up_dist;
        double khl_query_time = -GetCurrentTimeInSec();
        kpll->KDistanceQuery(u, v, up_dist);
        khl_query_time += GetCurrentTimeInSec();
        khl_time.push_back(khl_query_time);
        ++query_bar;
    }
    std::cout << "Writing on csv file\n";
    for(size_t j = 0; j < num_insertions; j++) {
        ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << j << ","
            << added_edges[j].first << "," << added_edges[j].second << ","  << update_loops[j] << "," << update_lengths[j] << "," << index_size[j] << ","
            << 0 << ","
            << 0 << "," << affected_hubs[j] << "," << reached_nodes[j] <<"\n";
    }

    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << 20000 << ","
        << "none" << "," << "none" << ","  << "final" << "," << "final" << "," << kpll->IndexSize() << ","
        << average(khl_time) << ","
        << median(khl_time) << "," << kpll->AffectedHubs() << "," << kpll->ReachedNodes() <<"\n";
    std::cout << "Writing done!\n";
    ofs.close();


    exit(EXIT_SUCCESS);
}
