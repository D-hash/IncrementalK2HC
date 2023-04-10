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

void bridgeUtil(NetworKit::Graph& gr, uint32_t u, bool visited[], uint32_t disc[], uint32_t low[], int32_t parent[],
                vector<pair<uint32_t, uint32_t>>& bridge_edges){
    static int time = 0;
    visited[u] = true;
    disc[u] = low[u] = ++time;
    assert(gr.degree(u) > 0);
    assert(gr.hasNode(u));
    gr.forNeighborsOf(u, [&](NetworKit::node v) {
        if(!visited[v]){
            parent[v] = u;
            bridgeUtil(gr, v, visited, disc, low, parent, bridge_edges);
            low[u]  = min(low[u], low[v]);
            if (low[v] > disc[u])
                bridge_edges.emplace_back(u,v);
        }
        else {
            if(v != parent[u])
                low[u]  = min(low[u], disc[v]);
        }
    });
}

void bridges(NetworKit::Graph& gr, vector<pair<uint32_t, uint32_t>>& bridge_edges){
    bool* visited = new bool[gr.numberOfNodes()];
    uint32_t* disc = new uint32_t[gr.numberOfNodes()];
    uint32_t* low = new uint32_t[gr.numberOfNodes()];
    int32_t * parent = new int32_t[gr.numberOfNodes()];
    for (uint32_t i = 0; i < gr.numberOfNodes(); i++)
    {
        parent[i] = -1;
        visited[i] = false;
    }
    for(uint32_t i = 0; i < gr.numberOfNodes(); i++){
        if(!visited[i]){
            bridgeUtil(gr, i, visited, disc, low, parent, bridge_edges);
        }
    }
}

void read_graph(const string graph_file, vector<pair<int, int> > &es){
    ifstream ifs(graph_file);
    if (!ifs.good()){
        cerr << "Error: Cannot open " << graph_file << "." << endl;
        exit(EXIT_FAILURE);
    }

    es.clear();
    for (int u, v; ifs >> u >> v;){
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
    std::cout << "Graph with " << g.numberOfNodes() << " vertices and " << g.numberOfEdges() << " edges.\n";
    vector<pair<uint32_t, uint32_t> > removed_edges;
    //vector<pair<uint32_t, uint32_t> > be;
    //bridges(g, be);
    long long int num_insertions = std::min((long long int)input_ins,
                                            (long long int)
                                                    (g.numberOfEdges()));
    std::cout << "Removing " << num_insertions << " edges\n";
//    removed_edges.emplace_back(2,5);
//    removed_edges.emplace_back(6,9);
//    g.removeEdge(2,5);
//    g.removeEdge(6,9);
    for(size_t i = 0; i < num_insertions; i++){
        uint32_t a = NetworKit::GraphTools::randomNode(g);
        uint32_t b = NetworKit::GraphTools::randomNeighbor(g,a);

        g.removeEdge(a,b);
//        while(find(be.begin(), be.end(), make_pair(a,b)) != be.end() || find(be.begin(), be.end(), make_pair(b,a)) != be.end()){
//            g.addEdge(a,b);
//
//            a = NetworKit::GraphTools::randomNode(g);
//            b = NetworKit::GraphTools::randomNeighbor(g,a);
//            g.removeEdge(a,b);
//        }
//        bic = new NetworKit::ConnectedComponents(g);
//        bic->run();
//        assert(bic->numberOfComponents() == 1);
        //be.clear();
        //bridges(g, be);
        bic = new NetworKit::ConnectedComponents(g);
        bic->run();
        while(bic->numberOfComponents()>1){
            g.addEdge(a,b);
            a = NetworKit::GraphTools::randomNode(g);
            b = NetworKit::GraphTools::randomNeighbor(g,a);
            g.removeEdge(a,b);
            bic = new NetworKit::ConnectedComponents(g);
            bic->run();
        }
        removed_edges.emplace_back((uint32_t) a, (uint32_t) b);
    }
    std::cout << "Edges after removal " << g.numberOfEdges() << '\n';
//    bic = new NetworKit::ConnectedComponents(g);
//    bic->run();
//    assert(bic->numberOfComponents() == 1);
    IncrementalTopK* kpll = new IncrementalTopK();
    kpll->ConstructIndex(g, K, directed);
    std::cout << "First Labeling Loop time: " << kpll->LoopCountTime() << "s | First Labeling Indexing time:" << kpll->IndexingTime()
              << "\n";
    std::cout << "Number Vertices: " << kpll->NumOfVertex() << "\n";

    std::cout << "Number of insertions " << num_insertions << "\n";

    std::ofstream ofs;
    ofs.open(graph_file+"_"+std::to_string(K)+"_"+std::to_string(num_insertions)+"_03_prefetch.csv");
//    ofs << "Graph,Vertices,Edges,K,Insertions,NewEdgeX,NewEdgeY,SLLoopTime,"
//           "SLLabelingTime,SLSize,ULLoopTime,ULLabelingTime,ULSize,DiffAvgIndexSize,ULMeanQueryTime,SLMeanQueryTime,"
//           "ULMedianQueryTime,SLMedianQueryTime,AffectedHubs,ReachedNodes,AffectedCycles,ReachedMBFS\n";
        ofs << "Graph,Vertices,Edges,K,Insertions,NewEdgeX,NewEdgeY,"
           "ULLoopTime,ULLabelingTime,ULSize,AvgIndexSize,ULMeanQueryTime,"
           "ULMedianQueryTime,AffectedHubs,ReachedNodes\n";
    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << 0 << ","
        << 0 << "," << 0 << ","  << kpll->LoopCountTime() << "," << kpll->IndexingTime() << "," << kpll->IndexSize() << ","
        << kpll->AvgIndexSize() << "," << 0 << "," << 0 << "," << 0 << "," << 0 <<"\n";
    uint16_t i = 0;
    int num_queries = 100000;
    while(!removed_edges.empty()) {
        i++;
        pair<uint32_t, uint32_t> new_edge = removed_edges.back();
        removed_edges.pop_back();
        uint32_t a = new_edge.first;
        uint32_t b = new_edge.second;
        uint32_t t = min(a, b);
        b = max(a, b);
        a = t;
        std::cout << "New edge " << a << " " << b << "\n";
        g.addEdge(a,b);
//        IncrementalTopK scratch_kpll;
//        scratch_kpll.ConstructIndex(g, K, directed);
//        std::cout << "Scracth Loop time: " << scratch_kpll.LoopCountTime() << "s | Scratch Indexing time:"
//                  << scratch_kpll.IndexingTime()
//                  << "\n";
        std::cout << "Updating first labeling.. \n";
        double ul_loops = -GetCurrentTimeInSec();
        kpll->UpdateLoops(make_pair(a,b));
        ul_loops += GetCurrentTimeInSec();
        double ul_labeling = -GetCurrentTimeInSec();
        kpll->UpdateIndex(make_pair(a,b));
        ul_labeling += GetCurrentTimeInSec();
        std::cout << "Update index time: " << ul_labeling << " | Update loops time: " << ul_loops << "\n";

        vector<double> sl_time, ul_time, bfs_time;
        ProgressStream query_bar(num_queries);
        query_bar.label() << "Queries";
        for(int j=0; j<num_queries; j++){
//            int32_t u = 3;
//            int32_t v = 9;
            int32_t u = rand() % kpll->NumOfVertex();
            int32_t v = rand() % kpll->NumOfVertex();
            vector<int> up_dist;
            vector<int> sc_dist;
            double up = -GetCurrentTimeInSec();
            kpll->KDistanceQuery(u, v, up_dist);
            up += GetCurrentTimeInSec();
            ul_time.push_back(up);
            double scratch = -GetCurrentTimeInSec();
            //scratch_kpll.KDistanceQuery(u, v, sc_dist);
            scratch += GetCurrentTimeInSec();
            sl_time.push_back(scratch);
//            assert(up_dist.size() == sc_dist.size());
//            for(size_t l=0; l < up_dist.size(); l++){
//                if(up_dist[l] != sc_dist[l]){
//                    std::cout << "Error bw " << u << "-" << v << "\n";
//                    std::cout << "Updated labeling distance: " << up_dist[l] << "\n";
//                    std::cout << "Scratch labeling distance: " << sc_dist[l] << "\n";
//                    for(size_t id=0; id < up_dist.size(); id++){
//                        std:: cout << "Up " << up_dist[id] << " | Scratch " << sc_dist[id] << "\n";
//                    }
//                    assert(false);
//                }
//            }
            ++query_bar;
        }
        std::cout << i << "-th insertion correct!" << "\n";
        std::cout << "UL avg query time " << average(ul_time) << "\n";
//        std::cout << "SL avg query time " << average(sl_time) << "\n";
//        ofs << graph_file << "," << kpll->NumOfVertex() << "," << g.numberOfEdges() << "," << K << "," << i+1 << ","
//            << a << "," << b << "," << scratch_kpll.LoopCountTime() << "," << scratch_kpll.IndexingTime() << ","
//            << scratch_kpll.IndexSize() << "," << ul_loops << "," << ul_labeling << "," << kpll->IndexSize() << ","
//            << kpll->AvgIndexSize() - scratch_kpll.AvgIndexSize() << "," << average(ul_time) << ","
//            << average(sl_time) << ","  << median(ul_time) << "," << median(sl_time) << "," << kpll->AffectedHubs()
//            << "," << kpll->ReachedNodes() << "," << kpll->AffectedCycles() << "," << kpll->ReachedMBFS() << "\n";
        ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << i << ","
            << a << "," << b << ","  << ul_loops << "," << ul_labeling << "," << kpll->IndexSize() << ","
            << kpll->AvgIndexSize() << "," << average(ul_time) << ","
            << median(ul_time) << "," << kpll->AffectedHubs() << "," << kpll->ReachedNodes() <<"\n";
    }

    IncrementalTopK scratch_kpll;
    scratch_kpll.ConstructIndex(g, K, directed);
    std::cout << "Scratch LB Loop time: " << scratch_kpll.LoopCountTime() << "s | Scratch LB Indexing time:"
              << scratch_kpll.IndexingTime()
              << "\n";
    vector<double> sl_time;
    ProgressStream query_bar(num_queries);
    query_bar.label() << "Queries";
    for(int j=0; j<num_queries; j++){
        int32_t u = rand() % kpll->NumOfVertex();
        int32_t v = rand() % kpll->NumOfVertex();
        vector<int> up_dist;
        vector<int> sc_dist;
        kpll->KDistanceQuery(u, v, up_dist);
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
    ofs << graph_file << "," << kpll->NumOfVertex() << "," << kpll->graph.numberOfEdges() << "," << K << "," << num_insertions << ","
        << "scratch" << "," << "scratch" << ","  << scratch_kpll.LoopCountTime() << "," << scratch_kpll.IndexingTime() << "," << scratch_kpll.IndexSize() << ","
        << scratch_kpll.AvgIndexSize() << "," << average(sl_time) << ","
        << median(sl_time) << ",scratch,scratch\n";

    ofs.close();

    exit(EXIT_SUCCESS);
}
