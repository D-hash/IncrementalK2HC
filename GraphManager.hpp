#include <fstream>
#include "progressBar.h"
#include <networkit/graph/GraphTools.hpp>           
#include <networkit/graph/Graph.hpp>
#include <networkit/io/EdgeListReader.hpp>
#ifndef GRAPHMANAGER_HPP_
#define GRAPHMANAGER_HPP_

class GraphManager
{
public:
    static void read_hist(std::string source, NetworKit::Graph **g, std::vector<std::pair<uint64_t, std::pair<uint32_t , uint32_t>> > *edge_structure)
    {
        std::ifstream ifs(source);
        if (!ifs)
            throw std::runtime_error("Error opening File ");

        int vertices = -1, edges = -1, weighted = -1, directed = -1;

        ifs >> vertices >> edges >> weighted >> directed;

        assert((weighted == 0 || weighted == 1) && (directed == 0 || directed == 1) && (vertices >= 0 && edges >= 0));

        ProgressStream reader(edges);
        std::string t1 = weighted == 0 ? " unweighted" : " weighted";
        std::string t2 = directed == 0 ? " undirected" : " directed";

        reader.label() << "Reading" << t1 << t2 << " graph in " << source << " (HIST FORMAT) containing " << vertices << " vertices and " << edges << " edges ";
        NetworKit::Graph *graph = new NetworKit::Graph(vertices, weighted, directed);
        int time, v1, v2, weight;

        while (true)
        {

            ifs >> time >> v1 >> v2 >> weight;
            if (ifs.eof())
                break;

            ++reader;

            //assert(weighted == 1 || weight == 1 || weight == -1);

            if (v1 == v2)
                continue;
            assert(graph->hasNode(v1) && graph->hasNode(v2));
            if (graph->hasEdge(v1, v2))
                // std::cout<<"SKIPPING ROW"<<std::endl;
                ;
            else
            {
                graph->addEdge(v1, v2, weight);
                edge_structure->emplace_back(time, std::make_pair(v1,v2));
#ifndef NDEBUG
                if (!directed)
                {
                    if (!graph->hasEdge(v1, v2) && !graph->hasEdge(v2, v1))
                        throw std::runtime_error("wrong edge insertion during construction");
                }
                else
                {
                    if (!graph->hasEdge(v1, v2))
                        throw std::runtime_error("wrong edge insertion during construction");
                }
#endif
            }
        }

        ifs.close();
        graph->indexEdges();
        *g = graph;
    }

    static void read_nde(std::string source, NetworKit::Graph **g)
    {
        std::ifstream ifs(source);

        if (!ifs)
            throw std::runtime_error("Error opening File ");
        
        
        NetworKit::Graph *graph = new NetworKit::Graph();
        const char x = ' ';
        NetworKit::EdgeListReader reader(x,0,"#",true,false);
        *graph = reader.read(source);
        *g = graph;
        graph->indexEdges();
        std::cout << "Read graph in " << source << " (EDGELIST FORMAT) containing " << graph->numberOfNodes() << " vertices and " << graph->numberOfEdges() << " edges\n";

    }
};

#endif /* GRAPHMANAGER_HPP_ */
