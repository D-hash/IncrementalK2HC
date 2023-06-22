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
#include <boost/program_options.hpp>
#include "networkit/components/ConnectedComponents.hpp"
#include <networkit/graph/GraphTools.hpp>
#include "mytimer.h"
#include "GraphManager.hpp"
#include <networkit/distance/Diameter.hpp>

using namespace std;
double median(std::vector<double>& arr) { //SORTS
    size_t n = arr.size() / 2;
    if (n % 2 == 0) {
        std::nth_element(arr.begin(),arr.begin() + n/2,arr.end());
        std::nth_element(arr.begin(),arr.begin() + (n - 1) / 2,arr.end());
        return (double) (arr[(n-1)/2]+ arr[n/2])/2.0;
    }

    else{
        std::nth_element(arr.begin(),arr.begin() + n / 2,arr.end());
        return (double) arr[n/2];
    }
    assert(false);
}



double average(std::vector<double> & arr) {

    auto const count = static_cast<double>(arr.size());
    double sum = 0;
    for(double value: arr) sum += value;
    return sum / count;
}


int main(int argc, char **argv) {
    srand (time(NULL));
    
    //declare supported options
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	
	desc.add_options()
	("graph_location,g", po::value<std::string>(), "Input Graph File Location")
	("k_paths,k", po::value<int>(), "Number of Top Paths to Compute")
	("num_insertions,n", po::value<int>(), "Number of Insertions to Be Performed")
    ("num_queries,q", po::value<int>(), "Number of Queries to Be Performed")
	("directed,d",po::value<int>(), "[FALSE(0) TRUE(1)]")
	("ordering,o",po::value<int>(), "Type of Node Ordering [DEGREE(0) APPROX-BETWEENESS(1) k-PATH(2)]")
	;

    po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
    if(vm.empty()){
		std::cout << desc << "\n";
		throw std::runtime_error("empty argument array");
	}
	std::string graph_location;
	
    if(vm.count("graph_location")){
		graph_location = vm["graph_location"].as<std::string>();
    }

    
    int K = -1;

	if (vm.count("k_paths")){
		K = vm["k_paths"].as<int>();
    }

	if (K < 2){
		std::cout << desc << "\n";
		throw std::runtime_error("K must be at least 2");
	}

    int ordering = -1; 

	if (vm.count("ordering")){
		ordering = vm["ordering"].as<int>();
	}
	if(ordering != 0 && ordering != 1 && ordering != 2){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong ordering selection (0 or 1 or 2)");
	}

    int num_insertions = -1;
	
    if (vm.count("num_insertions")){
		num_insertions = vm["num_insertions"].as<int>();
	}
	if(num_insertions < 2){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong num_insertions");
	}

    int num_queries = -1;
	if (vm.count("num_queries")){
		num_queries = vm["num_queries"].as<int>();
	}
	if(num_queries < 2){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong num queries");
	}

	int directed = -1;
	if (vm.count("directed")){
		directed = vm["directed"].as<int>();
	}
	if(directed != 0 && directed !=1){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong directed selection");
	}
    
    if(directed==0){
        std::cout << "Graph is undirected\n";
    }
    else{ 
        throw std::runtime_error("not yet implemented");

    }
    
    std::cout << "Reading " << graph_location << " building with K = " << K << " num_queries = "<<num_queries<<" ordering "<<ordering<<"\n";


    // LETTURA
    NetworKit::Graph *graph = new NetworKit::Graph();


	if(graph_location.find(".hist") != std::string::npos){
		GraphManager::read_hist(graph_location,&graph);
	}
	else if(graph_location.find(".nde") != std::string::npos || graph_location.find(".el") != std::string::npos){
		GraphManager::read_nde(graph_location,&graph);

	}
	else{
		throw std::runtime_error("Unsupported graph file format");
	}

	*graph = NetworKit::GraphTools::toUnweighted(*graph);
	*graph = NetworKit::GraphTools::toUndirected(*graph);

	const NetworKit::Graph& graph_handle = *graph;
	NetworKit::ConnectedComponents *cc = new NetworKit::ConnectedComponents(graph_handle);
	*graph = cc->extractLargestConnectedComponent(graph_handle, true);
    graph->shrinkToFit();
    graph->indexEdges();
	std::cout<<"Graph after CC has "<<graph->numberOfNodes()<<" vertices and "<<graph->numberOfEdges()<<" edges\n";
    double density = NetworKit::GraphTools::density(*graph);
    NetworKit::Diameter *dm = new NetworKit::Diameter(graph_handle,NetworKit::DiameterAlgo::EXACT,0.0);
    dm->run();

    
    double diameter = dm->getDiameter().first;
    delete dm;
    
	std::cout<<"Density: "<<density<<"\n";
	std::cout<<"Diameter: "<<diameter<<"\n";
    
    if(num_insertions>=graph->numberOfEdges()){
        num_insertions=graph->numberOfEdges();
    }
    std::cout << "Number of insertions " << num_insertions << "\n";
    
    std::cout << "Removal of " << num_insertions << " edges\n";

    edge_id original_n_edges = graph->numberOfEdges();
    std::vector<std::pair<vertex,vertex>> removed_edges;
    
    uint16_t attempts = 0;
    for(size_t i = 0; i < num_insertions; i++){

        vertex a = NetworKit::GraphTools::randomNode(*graph);
        vertex b = NetworKit::GraphTools::randomNeighbor(*graph,a);

        graph->removeEdge(a,b);
        attempts = 0;
        cc->run();
        while(cc->numberOfComponents()>1){
            attempts+=1;
            graph->addEdge(a,b);
            a = NetworKit::GraphTools::randomNode(*graph);
            b = NetworKit::GraphTools::randomNeighbor(*graph,a);
            graph->removeEdge(a,b);
            cc->run();
            if(attempts>graph->numberOfEdges())
                throw new std::runtime_error("experiment fails, too many removals, cc cannot be preserved");
        }
        removed_edges.push_back(std::make_pair(a,b));
    }
    assert(graph->numberOfEdges()==original_n_edges-removed_edges.size());
    std::cout << "Edges after removal " << graph->numberOfEdges() << '\n';

    IncrementalTopK* kpll = new IncrementalTopK(graph, K, directed,ordering, false);

    kpll->build();

    std::cout << "First Labeling Loop time: " << kpll->loops_time << " s | First Labeling Indexing time: " << kpll->lengths_time<< " s\n";
    std::cout << "First Labeling Loop entries: " << kpll->loop_entries<<" First Labeling Length Entries: "<<kpll->length_entries<< "\n";
    std::cout << "Number of Vertices: " << graph->numberOfNodes() << "\n";

    std::string order_string;

    switch(ordering){
        case (0):
            order_string = "DEG";
            break;
        case (1):
            order_string = "BET";
            break;
        case (2):
            order_string = "KPT";
            break;
        default:
            throw new std::runtime_error("problem");

    
    }
    //timestring
    std::time_t rawtime;
    std::tm* timeinfo;
    char buffer [80];

    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    std::strftime(buffer,80,"%d-%m-%Y-%H-%M-%S",timeinfo);
    std::string tmp_time = buffer;   
	std::string shortenedName = graph_location.substr(0, 16);
    stringstream string_object_name;
    string_object_name<<K;
    std::string k_string,n_insert_string;
    string_object_name>>k_string;
    string_object_name<<num_insertions;
    string_object_name>>n_insert_string;
	std::string timestampstring = shortenedName+"_"+k_string+"_"+n_insert_string+"_"+order_string+"_"+tmp_time;

	std::string logFile = timestampstring +"_remove_add.csv";
    

    std::ofstream ofs(logFile);


    ofs << "G,V,E,K,Insertions,x,y,UTLoops,UTLengths,ULoopsSize,ULengthSize,UAvgQT,UMedQT,aff,reached\n";
    ofs << graph_location << "," << graph->numberOfNodes() << "," << graph->numberOfEdges() << "," << K << "," << 0 << ","
        << 0 << "," << 0 << ","  << kpll->loops_time << "," << kpll->lengths_time << "," << kpll->loop_entries << "," << kpll->length_entries << ","
        << 0 << "," << 0 << "," << 0 << "," << 0 <<"\n";
    

    std::vector<double> update_loops_times;
    std::vector<double> update_lengths_times;
    

    std::vector<size_t> index_loops_size;
    std::vector<size_t> index_lengths_size;

    
    std::vector<vertex> affected_hubs;
    std::vector<double> reached_nodes;
    std::vector<pair<vertex, vertex>> added_edges;


    mytimer time_counter;

    for(size_t t=0;t<removed_edges.size();t++){

        kpll->x = removed_edges[t].first;
        kpll->y = removed_edges[t].second;
        std::cout << "New edge " << kpll->x << " " << kpll->y << "\n";

        std::cout << "Updating loops...";
        time_counter.restart();
        kpll->update_loops();
        update_loops_times.push_back(time_counter.elapsed());
        std::cout<<"done! \nUpdate loops time: " << time_counter.elapsed() << "\n"<<std::flush;
        std::cout << "Updating lengths...";
        time_counter.restart();
        kpll->update_lengths();
        update_lengths_times.push_back(time_counter.elapsed());
        std::cout << "done! \nUpdate lengths time: " << time_counter.elapsed()<<"\n"<<std::flush;
        std::cout << t+1 << "-th insertion correct!" << "\n";

        index_loops_size.push_back(kpll->loop_entries);
        index_lengths_size.push_back(kpll->length_entries);

        affected_hubs.push_back(kpll->aff_hubs);
        reached_nodes.push_back(kpll->n_reached_nodes());
        added_edges.push_back(removed_edges[t]);
        
    }
    kpll->deallocate_aux();

    assert(added_edges.size()==removed_edges.size());

    IncrementalTopK* scratch_kpll = new IncrementalTopK(graph, K, directed, ordering, true);

    scratch_kpll->build();
    //OK, no updates
    scratch_kpll->deallocate_aux();

    std::cout << "From Scratch Loop time: " << scratch_kpll->loops_time << " s | From Scratch Indexing time: "<< scratch_kpll->lengths_time<< " s\n";
    std::cout << "From Scratch Loop entries: " << scratch_kpll->loop_entries<<" From Scratch Length Entries: "<<scratch_kpll->length_entries<< "\n";
    std::vector<double> sl_time;
    std::vector<double> khl_time;
    ProgressStream query_bar(num_queries);
    
    query_bar.label() << "Queries";
    size_t l;
    std::vector<dist> update_distances;
    std::vector<dist> from_scratch_distances;

    for(uint64_t j=0; j<num_queries; j++){
        vertex u = NetworKit::GraphTools::randomNode(*graph);
        vertex v = NetworKit::GraphTools::randomNode(*graph);

        
        time_counter.restart();
        kpll->query(u, v, update_distances);
        khl_time.push_back(time_counter.elapsed());
        time_counter.restart();

        scratch_kpll->query(u, v, from_scratch_distances);
        sl_time.push_back(time_counter.elapsed());

        if(update_distances.size() != from_scratch_distances.size()){
            throw new std::runtime_error("cardinality problem");
        }
        for(l=0; l < update_distances.size(); l++){
            if(update_distances[l] != from_scratch_distances[l]){
                std::cout << "Error bw " << u << "-" << v << "\n";
                std::cout << "Updated labeling distance: " << update_distances[l] << "\n";
                std::cout << "Scratch labeling distance: " << from_scratch_distances[l] << "\n";
                for(size_t id=0; id < update_distances.size(); id++){
                    std:: cout << "Up " << update_distances[id] << " | Scratch " << from_scratch_distances[id] << "\n";
                }
                throw new std::runtime_error("correctness problem");
            }
        }
        ++query_bar;
    }
    std::cout << "Writing CSV file...";
    for(size_t j = 0; j < num_insertions; j++) {
        ofs << graph_location << "," << graph->numberOfNodes() << "," << graph->numberOfEdges() << "," << K << "," << j << ","
            << added_edges[j].first << "," << added_edges[j].second << ","  << update_loops_times[j] << "," << update_lengths_times[j] << "," << index_loops_size[j] << "," << index_lengths_size[j] << ","
            << 0 << "," << 0 << "," << affected_hubs[j] << "," << reached_nodes[j] <<"\n";
    }
    ofs << graph_location << "," 
        << graph->numberOfNodes() << "," 
        << graph->numberOfEdges() << "," 
        << K << "," 
        << removed_edges.size()+1 << ","
        << "none" << "," << "none" << ","  << "final" << "," << "final" << "," << kpll->loop_entries << "," << kpll->length_entries << ","
        << average(khl_time) << ","
        << median(khl_time) << "," 
        << kpll->aff_hubs << "," 
        << kpll->n_reached_nodes() <<"\n";

    ofs << graph_location << "," 
        << graph->numberOfNodes() << "," 
        << graph->numberOfEdges() << "," 
        << K << "," 
        << num_insertions << ","
        << "scratch" << "," << "scratch" << ","  << scratch_kpll->loops_time << "," << scratch_kpll->lengths_time << "," << scratch_kpll->loop_entries << "," << scratch_kpll->length_entries << ","
        << average(sl_time) << ","
        << median(sl_time) << ",scratch,scratch\n";
    std::cout << "done!\n";
    ofs.close();    
    delete graph;
    delete kpll;
    delete scratch_kpll;


    exit(EXIT_SUCCESS);
}
