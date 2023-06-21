//
// created by anonym 20/06/23
//#
#include "incremental_topk.h"
#include <queue>
#include <set>
#include <algorithm>
#include <cassert>
#include <climits>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include "networkit/centrality/DegreeCentrality.hpp"

using namespace std;

const vertex IncrementalTopK::null_vertex = round(std::numeric_limits<vertex>::max()/2);
const dist IncrementalTopK::null_distance = round(std::numeric_limits<dist>::max()/2);

IncrementalTopK::IncrementalTopK(NetworKit::Graph* g, dist k_value, bool dir, dist ord){

    this->graph = g;
    this->K = k_value;
    this->directed = dir;
    this->ordering_type = ord;
    this->length_labels[0].clear();
    this->length_labels[1].clear();
    this->loops_time = 0.0;
    this->lengths_time = 0.0;
    this->ordering_rank.resize(graph->numberOfNodes(),{null_vertex,null_vertex});
    this->reverse_ordering.resize(graph->numberOfNodes(),null_vertex);
    this->ordering.resize(graph->numberOfNodes(),null_vertex);
    double centr_time = 0.0;
    this->loop_entries = 0;
    this->length_entries = 0;

    if(this->ordering_type==0){

        INFO("BY DEGREE");       
        mytimer local_constr_timer;
        local_constr_timer.restart();
        const NetworKit::Graph& hand = *graph;
        NetworKit::DegreeCentrality* rank = new NetworKit::DegreeCentrality(hand);
        rank->run();
        this->graph->forNodes([&] (vertex i){
            assert(graph->hasNode(i));
            this->ordering_rank[i]=std::make_pair(rank->score(i),i);

        });
        delete rank;
        centr_time = local_constr_timer.elapsed();


    }
    else if(this->ordering_type==1){
        INFO("BY APX BETW");
        mytimer local_constr_timer;
        double max_time = 30.0;
        double cumulative_time = 0.0;
        double fract = 0.66;
        double n_samples =  round(std::pow((double)graph->numberOfNodes(),fract));

        const NetworKit::Graph& hand = *graph;
        while(cumulative_time<max_time && n_samples<(double)graph->numberOfNodes()){
            local_constr_timer.restart();

            std::cout<<"fract: "<<fract<<" "<<n_samples<<" SAMPLES\n";
            NetworKit::EstimateBetweenness* rank = new NetworKit::EstimateBetweenness(hand,n_samples,false,true);

            
            rank->run();
            
            this->graph->forNodes([&] (vertex i){
                assert(graph->hasNode(i));
                assert(i<graph->numberOfNodes());
                this->ordering_rank[i]=std::make_pair(rank->score(i),i);

            });
            delete rank;
            cumulative_time+=local_constr_timer.elapsed();
            n_samples*=2;
        }
        centr_time = cumulative_time;


    }
    else{
        assert(this->ordering_type==2);
        INFO("BY kPATH");
        mytimer local_constr_timer;
        local_constr_timer.restart();
        const NetworKit::Graph& hand = *graph;

        NetworKit::KPathCentrality* rank = new NetworKit::KPathCentrality(hand,0.0,round(std::pow((double)graph->numberOfNodes(),0.3)));

            
        rank->run();
            
        this->graph->forNodes([&] (vertex i){
            assert(graph->hasNode(i));
            assert(i<graph->numberOfNodes());
            this->ordering_rank[i]=std::make_pair(rank->score(i),i);

        });
        delete rank;

        centr_time = local_constr_timer.elapsed();
        //TODO : USE IT


    }
        

    std::sort(this->ordering_rank.begin(), this->ordering_rank.end(), [](const std::pair<double,vertex>  &a, const std::pair<double,vertex>  &b) {
        if(a.first == b.first)
            return a.second > b.second;
        else{
            return a.first > b.first;
        }
    });
    

    for(size_t count = 0; count < graph->numberOfNodes();count++){
        this->reverse_ordering[count]=this->ordering_rank[count].second;
        this->ordering[this->ordering_rank[count].second]=count;
        

    }

    for(size_t count = 0; count < 10 ;count++)
        std::cout<<"In position "<<count<<" we have vertex "<<this->reverse_ordering[count]<<" rank "<<this->ordering_rank[count].first<<std::endl;

    this->ordering_rank.clear();

};
IncrementalTopK::~IncrementalTopK(){};

void IncrementalTopK::build(){

    
    tmp_pruned.resize(graph->numberOfNodes(), false);
    tmp_offset.resize(graph->numberOfNodes(), null_distance);
    tmp_count .resize(graph->numberOfNodes(), 0);

    tmp_s_offset.resize(graph->numberOfNodes(), null_distance); 
    tmp_s_offset.push_back(0);
    tmp_s_count .resize(graph->numberOfNodes());
    
    visited_in_update_loops.resize(graph->numberOfNodes(), null_distance);

    for (size_t j = 0; j < 2; j++){
        tmp_dist_count[j].resize(graph->numberOfNodes(), 0);
    }

    loop_labels.resize(graph->numberOfNodes());
    #ifndef NDEBUG
        for (size_t v = 0; v < graph->numberOfNodes(); v++)
        assert(loop_labels[v].size()==0);
    #endif

    for (size_t dir = 0; dir < 1 + directed; dir++){
        length_labels[dir].resize(graph->numberOfNodes());

        for (size_t v = 0; v < graph->numberOfNodes(); v++){
            length_labels[dir][v].label_offset.clear();
            length_labels[dir][v].d_array.clear();
            length_labels[dir][v].d_array.clear();
        }
    }

    reached_mbfs.clear();

    
    
    mytimer build_timer;
    build_timer.restart();
    ProgressStream loop_bar(graph->numberOfNodes());

    loop_bar.label() << "Loops construction";
    
    for(size_t v = 0; v < graph->numberOfNodes(); v++){
        this->compute_loop_entries(v);
        ++loop_bar;
    }


    loops_time = build_timer.elapsed();

    build_timer.restart();
    ProgressStream index_bar(graph->numberOfNodes());
    index_bar.label() << "Index construction";
    
    for(size_t v = 0; v < graph->numberOfNodes(); v++){

        this->pruned_bfs(v, false);
        ++index_bar;
        if (directed){
            this->pruned_bfs(v, true);
        }
    }
    lengths_time = build_timer.elapsed();

    #ifndef NDEBUG
        verify_sizes();
    #endif

}

void IncrementalTopK::query(vertex s, vertex t, std::vector<dist> & container){
    container.clear();

    s = ordering[s];
    t = ordering[t];
    size_t pos1 = 0;
    size_t pos2 = 0;

    std::vector<dist> count(30, 0);
    const index_t &ids = length_labels[directed][s];
    const index_t &idt = length_labels[0][t];
    vertex W;
    dist d_tmp, c_tmp;
    for (;;){
        if (pos1 >= ids.label_offset.size()){
            break;
        }
        if (pos2 >= idt.label_offset.size()){
            break;
        }
        if (ids.label_offset[pos1].first == idt.label_offset[pos2].first){
            W = ids.label_offset[pos1].first;
            
            for (size_t i = 0; i < ids.d_array[pos1].size(); i++){
                for (size_t j = 0; j < idt.d_array[pos2].size(); j++){
                    for (size_t m = 0; m < loop_labels[W].size(); m++){
                        d_tmp = ids.label_offset[pos1].second + idt.label_offset[pos2].second + i + j + m;
                        c_tmp = loop_labels[W][m] - (m ? loop_labels[W][m-1] : 0);
                        if (count.size() <= d_tmp){
                            count.resize(d_tmp + 1, 0);
                        }
                        count[d_tmp] += (vertex)ids.d_array[pos1][i] * idt.d_array[pos2][j] * c_tmp;
                    }
                }
            }
            pos1++;
            pos2++;
        } 
        else {
            if (ids.label_offset[pos1].first < idt.label_offset[pos2].first){
                pos1++;
            } 
            else {
                pos2++;
            }
        }
    }

    for (size_t i = 0; i < count.size(); i++){
        while (container.size() < this->K && count[i]-- > 0){
            container.push_back(i);
        }
    }

    //return container.size() < this->K ? INT_MAX : 0;
}

void IncrementalTopK::verify_sizes(){
    
    vertex sz = 0;
    for(vertex i=0;i<graph->numberOfNodes();i++){
        assert(graph->hasNode(i));
        sz += loop_labels[i].size(); // loopcount

    }
    assert(sz==loop_entries);
    sz = 0;
    for (size_t dir = 0; dir < 1 + directed; dir++){
        for(vertex i=0;i<graph->numberOfNodes();i++){
            assert(graph->hasNode(i));
            sz += length_labels[dir][i].label_offset.size();
            for (size_t j = 0; j < length_labels[dir][i].d_array.size(); j++){
                sz += length_labels[dir][i].d_array[j].size();
            }
        }

    }
    // std::cout<<sz<<" "<<length_entries<<"\n";
    assert(sz==length_entries);
    
}


void IncrementalTopK::deallocate_aux() {
    visited_in_update_loops.clear();
    tmp_pruned.clear();
    tmp_offset.clear();
    tmp_count .clear();

    for (int i = 0; i < 2; i++){
        tmp_dist_count[i].clear();
    }
    tmp_s_offset.clear();
    tmp_s_count .clear();
}
void IncrementalTopK::compute_loop_entries(vertex s){

    size_t  count = 0;
    vertex     curr  = 0;
    vertex     next  = 1;
    dist distance  = 0;

    std::queue<vertex> node_que[2];
    
    this->updated.clear();

    node_que[curr].push(s);
    this->updated.push_back(s);
    tmp_dist_count[curr][s] = 1;
    
    vertex currently_reached_nodes = 0;
    vertex to_v;
    vertex v;
    dist c;

    for (;;){


        while (!node_que[curr].empty() && count < this->K){
            v = node_que[curr].front(); node_que[curr].pop();
            c = tmp_dist_count[curr][v]; // the number of path from s to v with dist hops.
            tmp_dist_count[curr][v] = 0;
            if (c == 0) {
                continue;
            }

            if (v == s){
                int old_size = loop_labels[s].size();
                
                loop_labels[s].resize(distance + 1, 0);
                loop_labels[s][distance] += c;


                loop_entries+=(loop_labels[s].size()-old_size);
                count += c;
            }
            currently_reached_nodes++;

            for(vertex u : graph->neighborRange(reverse_ordering[v])){

                to_v = ordering[u];
                if (tmp_count[to_v] == 0){
                    this->updated.push_back(to_v);
                }

                if (to_v >= s && tmp_count[to_v] < this->K){
                    tmp_count[to_v] += c;
                    node_que[next].push(to_v);
                    tmp_dist_count[next][to_v] += c;
                }
            }
        }
        if(node_que[next].empty() || count >= K){
            break;
        }
        std::swap(curr, next);
        distance++;
    }

    for(size_t i = 1; i < loop_labels[s].size(); i++){
        loop_labels[s][i] += loop_labels[s][i-1];
    }
    assert(loop_labels[s][0] == 1);
    
    reached_mbfs.push_back(currently_reached_nodes);
    
    reset_temp_vars(s, false);
}

void IncrementalTopK::pruned_bfs(vertex s, bool reversed){
    
    set_temp_vars(s, reversed);

    vertex curr = 0;
    vertex     next = 1;
    dist distance = 0;

    std::queue<vertex> node_que[2];
    this->updated.clear();

    node_que[curr].push(s);
    tmp_dist_count[curr][s] = 1;
    this->updated.push_back(s);

    for (;;){

        while (!node_que[curr].empty()){

            vertex v = node_que[curr].front(); node_que[curr].pop();
            dist c = tmp_dist_count[curr][v];
            tmp_dist_count[curr][v] = 0;

            if(c == 0 || tmp_pruned[v]){ 
                continue;
            }
            tmp_pruned[v] = prune(v, distance, reversed);

            if(tmp_pruned[v]){
                continue;
            }

            if(tmp_offset[v] == null_distance){
                tmp_offset[v] = distance;
                allocate_label(v, s, distance, c, reversed);
            }
            else{
                extend_label(v, s, distance, c, reversed, 0);
            }

            graph->forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
                vertex to = ordering[u];
                if(tmp_count[to] == 0){
                    this->updated.push_back(to);
                }

                if(to > s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            });
        }

        if (node_que[next].empty()){
            break;
            }
        swap(curr, next);
        distance++;
    }
    reset_temp_vars(s, reversed);
};

void IncrementalTopK::update_loops() {
    
    std::set<vertex> to_update;
    std::set<vertex> reset_visited;
    std::queue<vertex> q;

    vertex dequeued_v;
    q.push(ordering[this->x]);
    visited_in_update_loops[ordering[this->x]] = 0;
    while(!q.empty()){
        dequeued_v = q.front();
        q.pop();
        if(visited_in_update_loops[dequeued_v] > K){
            continue;
        }
        aff_cycles++;
        to_update.insert(dequeued_v);
        graph->forNeighborsOf(reverse_ordering[dequeued_v], [&](NetworKit::node u) {
            vertex to = ordering[u];
            if(visited_in_update_loops[to] == null_distance && to != ordering[this->y]){
                q.push(to);
                visited_in_update_loops[to] = visited_in_update_loops[dequeued_v] + 1;
                if(reset_visited.size() != graph->numberOfNodes()){
                    reset_visited.insert(to);
                }
            }
        });
    }
    std::queue<vertex> empty;
    std::swap( q, empty );
    
    vertex to_v;
    
    q.push(ordering[this->y]);

    visited_in_update_loops[ordering[this->y]] = 0;

    while(!q.empty()){
        dequeued_v = q.front(); 
        q.pop();
        if(visited_in_update_loops[dequeued_v] > K){
            continue;
        }
        aff_cycles++;
        to_update.insert(dequeued_v);
        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){

            to_v = ordering[u];
            if(visited_in_update_loops[to_v] == null_distance){
                q.push(to_v);
                visited_in_update_loops[to_v] = visited_in_update_loops[dequeued_v] + 1;
                if(reset_visited.size() != graph->numberOfNodes()){
                    reset_visited.insert(to_v);
                }
            }
        }
    }
    assert(!graph->hasEdge(this->x,this->y));

    graph->addEdge(this->x,this->y);

    vertex min_order = min(ordering[this->x], ordering[this->y]);
    
    aff_cycles = to_update.size();
    
    reached_mbfs.clear();

    vertex ordered_degree = 0;
    for(vertex u: to_update){
        if(u > min_order){
            continue;
        }

        ordered_degree = 0;
        for(vertex neighbor : graph->neighborRange(reverse_ordering[u])){

            if(u < ordering[neighbor]){
                ordered_degree++;
            }
            if(ordered_degree >= K){
                continue;
            }
        }

        if(ordered_degree >= K){
            continue;
        }
        loop_labels[u].clear();
        
        this->compute_loop_entries(u);
    }

    for(auto v: reset_visited) {
        visited_in_update_loops[v] = null_distance;
    }
    visited_in_update_loops[ordering[this->x]] = null_distance;
    visited_in_update_loops[ordering[this->y]] = null_distance;
    #ifndef NDEBUG
        for(auto v: visited_in_update_loops) 
            assert(v == null_distance);
    #endif
}

void IncrementalTopK::update_lengths() {




    aff_hubs = 0;
    reached_nodes.clear();

    const index_t &idva = length_labels[0][ordering[this->x]];
    const index_t &idvb = length_labels[0][ordering[this->y]];
    
    

    
    // std::cout<< "Updating affected hubs ...";
    
    this->old_label_a.clear();
    this->old_label_b.clear();
    this->old_distances_a.clear();
    this->old_distances_b.clear();
    vertex w_a, w_b;

    for(auto& elem: idva.label_offset){
        this->old_label_a.push_back(elem);
    }
    
    for(auto& elem: idvb.label_offset){
        this->old_label_b.push_back(elem);
    }
    
    this->old_distances_a.resize(idva.d_array.size());
    this->old_distances_b.resize(idvb.d_array.size());
    
    for(size_t i = 0; i < idva.d_array.size(); i++){
        for(size_t j = 0; j < idva.d_array[i].size(); j++){
            this->old_distances_a[i].push_back(idva.d_array[i][j]);
        }
    }
    for(size_t i = 0; i < idvb.d_array.size(); i++){
        for(size_t j = 0; j < idvb.d_array[i].size(); j++){
            this->old_distances_b[i].push_back(idvb.d_array[i][j]);
        }
    }
    size_t pos_a = 0;
    size_t pos_b = 0;
    while (pos_a != this->old_label_a.size() && pos_b != this->old_label_b.size()){
        
        new_labels.clear();
        
        w_a = pos_a < this->old_label_a.size() ? this->old_label_a[pos_a].first : graph->numberOfNodes();
        w_b = pos_b < this->old_label_b.size() ? this->old_label_b[pos_b].first : graph->numberOfNodes();
        
        if(w_a < w_b){
            if(w_a < ordering[this->y]){
                aff_hubs++;
                for(size_t i = 0; i < this->old_distances_a[pos_a].size(); i++){
                    for(size_t j = 0; j < this->old_distances_a[pos_a][i]; j++){
                        resume_pbfs(w_a,ordering[this->y], i+this->old_label_a[pos_a].second+1, false);
                    }
                }
                reset_temp_vars(w_a, directed);
            }
            pos_a++;

        }
        else if(w_b < w_a){
            if(w_b < ordering[this->x]){
                aff_hubs++;
                for(size_t i = 0; i < this->old_distances_b[pos_b].size(); i++){
                    for(size_t j = 0; j < this->old_distances_b[pos_b][i]; j++){
                        resume_pbfs(w_b,ordering[this->x], i+this->old_label_b[pos_b].second+1, false);
                    }
                }
                reset_temp_vars(w_b, directed);
            }
            pos_b++;

        }
        else {
            aff_hubs++;

            if(w_a < ordering[this->y]){
                for(size_t i = 0; i < this->old_distances_a[pos_a].size(); i++){
                    for(size_t j = 0; j < this->old_distances_a[pos_a][i]; j++){
                        resume_pbfs(w_a,ordering[this->y], i+this->old_label_a[pos_a].second+1, false);
                    }
                }
                reset_temp_vars(w_a, directed);
            }
            if(w_b < ordering[this->x]){
                for(size_t i = 0; i < this->old_distances_b[pos_b].size(); i++){
                    for(size_t j = 0; j < this->old_distances_b[pos_b][i]; j++){
                        resume_pbfs(w_b,ordering[this->x], i+this->old_label_b[pos_b].second+1, false);
                    }
                }
                reset_temp_vars(w_b, directed);
            }
            pos_a++;
            pos_b++;
        }
        
        for(size_t t=0;t<new_labels.size();t++){
            extend_label_repair(std::get<0>(new_labels[t]), std::get<1>(new_labels[t]), std::get<2>(new_labels[t]), std::get<3>(new_labels[t]), std::get<4>(new_labels[t]));
        }
    }
    // std::cout<<"done!\n";
}





inline bool IncrementalTopK::prune(vertex v,  dist d, bool rev){
    const index_t &idv = length_labels[rev][v];

    size_t pcount = 0;

    // cerr << "prune start" << endl;
    for (size_t pos = 0; pos < idv.label_offset.size(); pos++){
        vertex w = idv.label_offset[pos].first;

        if (tmp_s_offset[w] == null_distance) continue;

        const vector<dist> &dcs = tmp_s_count[w];

        int l = dcs.size() - 1;
        int c = d - tmp_s_offset[w] - idv.label_offset[pos].second;

        // By using precomputed table tmp_s_count, compute the number of path with a single loop.
        for (int i = 0; i <= c && i < idv.d_array[pos].size(); i++){
            pcount += (int)dcs[std::min(c - i, l)] * idv.d_array[pos][i];
        }

        if (pcount >= K) return true;
    }
    return false;
}
void IncrementalTopK::resume_pbfs(vertex s, vertex t, dist d, bool rev) {
    set_temp_vars(s, rev);

    vertex     curr = 0;
    vertex     next = 1;
    dist distance = d;

    std::queue<vertex> node_que[2];

    this->updated.clear();

    node_que[curr].push(t);
    tmp_dist_count[curr][t] = 1;
    
    this->updated.push_back(t);
    
    vertex currently_reached_nodes = 0;
    std::vector<dist> dists;

    for (;;){

        while (!node_que[curr].empty()){

            vertex v = node_que[curr].front(); node_que[curr].pop();
            dist  c = tmp_dist_count[curr][v];
            tmp_dist_count[curr][v] = 0;

            if(c == 0 || tmp_pruned[v]){
                continue;
            }
            this->query(reverse_ordering[s], reverse_ordering[v], dists);
            currently_reached_nodes += 1;
            tmp_pruned[v] = dists.size() == K && *dists.rbegin() <= distance;

            if(tmp_pruned[v]){continue;}

            new_labels.emplace_back(v, s, distance, c, rev, tmp_offset[v]);
            graph->forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
                vertex to = ordering[u];
                if(tmp_count[to] == 0){
                    this->updated.push_back(to);
                }

                if(to > s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            });
        }

        if (node_que[next].empty()){
            break;
        }
        swap(curr, next);
        distance++;
    }
    reached_nodes.push_back(currently_reached_nodes);
    reset_temp_vars(t, rev);
}

inline void IncrementalTopK::set_temp_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];
    vertex w;
    vector<dist> tmp_v;

    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        tmp_s_offset[w] = ids.label_offset[pos].second;

        tmp_v.clear();    
        for(size_t i = 0; i < ids.d_array[pos].size(); i++){
            tmp_v.push_back(ids.d_array[pos][i]);
        }
        tmp_s_count[w].resize(tmp_v.size() + loop_labels[w].size() - 1, 0);

        for(size_t i = 0; i < tmp_v.size(); i++){
            for(size_t j = 0; j < loop_labels[w].size(); j++){
                tmp_s_count[w][i+j] += tmp_v[i] * loop_labels[w][j];
            }
        }
    }
}

inline void IncrementalTopK::reset_temp_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];
    vertex w;
    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        tmp_s_offset[w] = null_distance;
        tmp_s_count[w].clear();
    }

    for(size_t i = 0; i < this->updated.size(); i++){
        tmp_count [this->updated[i]] = 0;
        tmp_offset[this->updated[i]] = null_distance;
        tmp_pruned[this->updated[i]] = false;
        for(size_t j = 0; j < 2; j++){
            tmp_dist_count[j][this->updated[i]] = 0;
        }
    }
}

inline void IncrementalTopK::allocate_label(vertex v, vertex start, dist distance, dist count, bool dir){
    index_t &idv = length_labels[dir][v];
    
    idv.label_offset.push_back(make_pair(start,distance));
    length_entries++;
    int old_size = idv.d_array.size();
    idv.d_array.resize(idv.d_array.size()+1);
    idv.d_array[idv.label_offset.size()-1].resize(1,count);
    
    length_entries+=(idv.d_array.size()-old_size);

}
double IncrementalTopK::n_reached_nodes(){
    double sum = 0.0; 
    for(auto& element:reached_nodes){
        sum+=(double)element;
    }
    return reached_nodes.size()>0 ? sum / (double)reached_nodes.size() : 0.0;
}
// double IncrementalTopK::n_reached_nodes_mbfs(){ 
//     double sum = 0.0; 
//     for(auto& element:reached_mbfs){
//         sum+=(double)element;
//     }
//     return reached_mbfs.size()>0 ? sum / (double) reached_mbfs.size() : 0.0;
// }



inline void IncrementalTopK::extend_label(vertex v, vertex start, dist distance, dist count, bool dir, size_t pos){
    index_t &idv = length_labels[dir][v];
    if(pos == 0){
        for(; pos < idv.label_offset.size(); pos++){
            if (idv.label_offset[pos].first == start){
                break;
            }
        }
    }
    dist offset = distance - idv.label_offset[pos].second;
    if (idv.d_array[pos].size() > offset){
        idv.d_array[pos][offset] += count;
    }
    else{
        while(idv.d_array[pos].size() != offset){
            idv.d_array[pos].push_back(0);
            length_entries++;
        }
        idv.d_array[pos].push_back(count);
        length_entries++;
    }
    int ol_size = 0;
    for(size_t p = 0; p < idv.label_offset.size(); p++){
        vertex tot_count = 0;
        for(size_t i = 0; i < idv.d_array[p].size(); i++){
            tot_count += idv.d_array[p][i];
            if(tot_count >= K){
                ol_size = idv.d_array[p].size();
                idv.d_array[p].resize(i+1);
                length_entries+=(idv.d_array[p].size()-ol_size);
                break;
            }
        }
    }
}

inline void IncrementalTopK::extend_label_repair(vertex v, vertex start, dist distance, dist count, bool dir){
    index_t &idv = length_labels[dir][v];

    size_t last = 0;
    for(; last < idv.label_offset.size(); last++) {
        if (idv.label_offset[last].first == start) {
            break;
        }
        if (idv.label_offset[last].first > start){
            break;
        }
    }
    if (last == idv.label_offset.size()){
        allocate_label(v, start, distance, count, dir);
        return;
    }
    else if(idv.label_offset[last].first == start){
        if(idv.label_offset[last].second > distance){
            idv.d_array[last].resize(idv.d_array[last].size() + idv.label_offset[last].second-distance, 0); // todo consider dequeue
            size_t rev = 0;
            for(rev = idv.d_array[last].size()-1; rev > 0; rev--){
                idv.d_array[last][rev] = idv.d_array[last][rev-(idv.label_offset[last].second-distance)];
            }
            for(rev = 1; rev != idv.label_offset[last].second-distance; rev++){
                idv.d_array[last][rev] = 0;
            }
            idv.d_array[last][0] = count;
            idv.label_offset[last].second = distance;
        } 
        else{
            extend_label(v, start, distance, count, dir, last);
        }
    } else {
        idv.label_offset.insert(idv.label_offset.begin()+last, std::make_pair(start,distance));
        idv.d_array.insert(idv.d_array.begin() + last, {count});
        for(size_t p = 0; p < idv.label_offset.size(); p++){
            vertex tot_count = 0;
            for(size_t i = 0; i < idv.d_array[p].size(); i++){
                tot_count += idv.d_array[p][i];
                if(tot_count >= K){
                    idv.d_array[p].resize(i+1);
                    break;
                }
            }
        }
    };

}

void IncrementalTopK::mod_bfs(vertex s, vertex t, std::vector<dist> &ret) {
    s = ordering[s];
    t = ordering[t];
    std:vector<std::vector<dist> > distanze(graph->numberOfNodes());
    std::priority_queue<pair<dist, vertex > > que;

    que.push(make_pair(0, s));
    vertex v;
    vertex dist_negative;
    while (!que.empty()) {
        v = que.top().second;
        dist_negative = -que.top().first;
        que.pop();

        if(distanze[t].size() >= this->K){
            break;
        }
        if (distanze[v].size() >= this->K){
            continue;
        }

        distanze[v].push_back(dist_negative);
        graph->forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
            
            que.push(make_pair(-(1 + dist_negative), ordering[u]));
        });
    }
    ret = distanze[t];
}

