//
// created by anonym 20/06/23
//#
#include "incremental_topk.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include "networkit/centrality/DegreeCentrality.hpp"

using namespace std;

const vertex IncrementalTopK::null_vertex = round(std::numeric_limits<vertex>::max()/2);
const dist IncrementalTopK::null_distance = round(std::numeric_limits<dist>::max()/2);

IncrementalTopK::IncrementalTopK(NetworKit::Graph* g, dist k_value, bool dir, dist ord, bool is_from_scratch){

    this->graph = g;
    this->K = k_value;
    this->directed = dir;
    this->ordering_type = ord;
    this->is_from_scratch_only = is_from_scratch;


    
    this->loops_time = 0.0;
    this->lengths_time = 0.0;
    
    this->ordering_rank = new std::pair<double,vertex>[graph->numberOfNodes()];
    this->ordering = new vertex[graph->numberOfNodes()];
    this->reverse_ordering = new vertex[graph->numberOfNodes()];

    this->graph->parallelForNodes([&] (vertex i){
        assert(graph->hasNode(i));
        this->ordering[i]=null_vertex;
        this->reverse_ordering[i]=null_vertex;
        this->ordering_rank[i] = {null_vertex,null_vertex};


    });
  
    
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
        

    std::sort(this->ordering_rank, this->ordering_rank+graph->numberOfNodes(), [](const std::pair<double,vertex>  &a, const std::pair<double,vertex>  &b) {
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

    
    delete[] this->ordering_rank;

};
IncrementalTopK::~IncrementalTopK(){
    delete[] this->ordering;
    delete[] this->reverse_ordering;
    delete[] this->loop_labels;
    


    for (size_t dir = 0; dir < 1 + directed; dir++){
        delete[] this->length_labels[dir];

    }
    delete[] this->length_labels;


};

void IncrementalTopK::deallocate_aux() {
    
    delete[] this->tmp_pruned;
    delete[] this->tmp_offset;
    delete[] this->tmp_count;
    delete[] this->tmp_s_offset;
    delete[] this->tmp_s_count;
    if(!this->is_from_scratch_only){
        delete[] this->visited_in_update_loops;
    }
    delete[] this->tmp_dist_count;



}

uint64_t IncrementalTopK::compute_index_size() {
    uint64_t sz = 0;
    for (size_t dir = 0; dir < 1 + directed; dir++){
        for(vertex i=0;i<graph->numberOfNodes();i++){
            sz += this->loop_labels[i].size() * sizeof(dist); // loopcount
            sz += this->length_labels[dir][i].label_offset.size() * (sizeof(vertex)+sizeof(dist)); // label_offset count
            sz += this->length_labels[dir][i].d_array.size() * sizeof(dist); // array of distances count
        }
    }
    return sz;
}

void IncrementalTopK::build(){

    this->tmp_pruned = new bool[graph->numberOfNodes()];
    this->tmp_offset = new dist[graph->numberOfNodes()];
    this->tmp_count = new vertex[graph->numberOfNodes()];
    this->tmp_dist_count = new std::vector<dist>[2];
    for (size_t j = 0; j < 2; j++){
        this->tmp_dist_count[j].resize(graph->numberOfNodes(), 0);
    }

    this->tmp_s_offset = new dist[graph->numberOfNodes()+1];
    this->tmp_s_count = new std::vector<dist>[graph->numberOfNodes()];

    this->graph->parallelForNodes([&] (vertex i){
        assert(graph->hasNode(i));
        this->tmp_pruned[i]=false;        
        this->tmp_offset[i]=null_distance;
        this->tmp_count[i]=0;

        this->tmp_s_offset[i] = null_distance;
        this->tmp_s_count[i].clear();

    });
    this->tmp_s_offset[graph->numberOfNodes()] = 0;
    
    if(!this->is_from_scratch_only){
        this->visited_in_update_loops = new dist[graph->numberOfNodes()];
        this->graph->parallelForNodes([&] (vertex i){
            assert(graph->hasNode(i));
            this->visited_in_update_loops[i]=null_distance;
        });
    }

    this->reached_mbfs.clear();

    this->updated.clear();

    this->length_labels = new index_t*[2];

    
    this->loop_labels = new std::vector<dist>[graph->numberOfNodes()];
    


    for (size_t dir = 0; dir < 1 + directed; dir++){
        this->length_labels[dir] = new index_t[graph->numberOfNodes()];

        for (size_t v = 0; v < graph->numberOfNodes(); v++){
            this->loop_labels[v].clear();
            this->length_labels[dir][v].label_offset.clear();
            this->length_labels[dir][v].d_array.clear();
            this->length_labels[dir][v].d_array.clear();
        }
    }

    #ifndef NDEBUG
        for (size_t v = 0; v < graph->numberOfNodes(); v++)
        assert(loop_labels[v].size()==0);
    #endif
    
    
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

inline void IncrementalTopK::compute_loop_entries(vertex s){

    size_t  count = 0;
    vertex     curr  = 0;
    vertex     next  = 1;
    dist distance  = 0;

    this->node_que = new std::queue<vertex>[2];
    
    assert(this->updated.empty());

    this->node_que[curr].push(s);
   
    this->updated.push_back(s);
    this->tmp_dist_count[curr][s] = 1;
    
    vertex num_reached = 0;
    vertex to_v;
    vertex v;
    dist c;

    for (;;){


        while (!this->node_que[curr].empty() && count < this->K){
            v = this->node_que[curr].front(); 
            this->node_que[curr].pop();
            c = this->tmp_dist_count[curr][v]; // the number of path from s to v with dist hops.
            this->tmp_dist_count[curr][v] = 0;

            if (c == 0) {
                continue;
            }

            if (v == s){
                int old_size = loop_labels[s].size();
                
                this->loop_labels[s].resize(distance + 1, 0);
                this->loop_labels[s][distance] += c;


                this->loop_entries+=(loop_labels[s].size()-old_size);
                count += c;
            }
            num_reached++;

            for(vertex u : graph->neighborRange(reverse_ordering[v])){

                to_v = ordering[u];
                if (this->tmp_count[to_v] == 0){
                    // assert(std::find(this->updated.begin(),this->updated.end(),to_v)==this->updated.end());
                    this->updated.push_back(to_v);
                }

                if (to_v >= s && this->tmp_count[to_v] < this->K){
                    this->tmp_count[to_v] += c;
                    this->node_que[next].push(to_v);
                    this->tmp_dist_count[next][to_v] += c;
                }
            }
        }
        if(this->node_que[next].empty() || count >= K){
            break;
        }
        std::swap(curr, next);
        distance++;
    }
    delete[] this->node_que;
    
    for(size_t i = 1; i < loop_labels[s].size(); i++){
        loop_labels[s][i] += loop_labels[s][i-1];
    }
    assert(loop_labels[s][0] == 1);
    
    this->reached_mbfs.push_back(num_reached);
    
    reset_temp_vars(s, false);
}


inline void IncrementalTopK::pruned_bfs(vertex s, bool reversed){

    allocate_label(s, s, 0, 1, reversed);
    this->tmp_offset[s] = 0;
    set_temp_vars(s, reversed);

    vertex curr = 0;
    vertex     next = 1;
    dist distance = 0;

    this->node_que = new std::queue<vertex>[2];
    this->updated.clear();

    this->node_que[curr].push(s);
    this->tmp_dist_count[curr][s] = 1;
    this->updated.push_back(s);
    vertex v,to_vert;
    dist c;

    for (;;){

        while (!this->node_que[curr].empty()){

            v = this->node_que[curr].front();
            this->node_que[curr].pop();
            c = this->tmp_dist_count[curr][v];
            this->tmp_dist_count[curr][v] = 0;

            if(c == 0 || tmp_pruned[v]){
                continue;
            }
            this->tmp_pruned[v] = prune(v, distance, reversed);

            if(this->tmp_pruned[v]){
                continue;
            }
            if(s!=v) {
                if (this->tmp_offset[v] == null_distance) {
                    this->tmp_offset[v] = distance;
                    allocate_label(v, s, distance, c, reversed);
                } else {
                    extend_label(v, s, distance, c, reversed, 0);
                }
            }
            for(vertex u : graph->neighborRange(reverse_ordering[v])){
                to_vert = ordering[u];
                if(this->tmp_count[to_vert] == 0){
                    this->updated.push_back(to_vert);
                }

                if(to_vert > s && this->tmp_count[to_vert] < K){
                    this->tmp_count[to_vert] += c;
                    this->node_que[next].push(to_vert);
                    this->tmp_dist_count[next][to_vert] += c;
                }
            }
        }

        if (this->node_que[next].empty()){
            break;
        }
        std::swap(curr, next);
        distance++;
    }
    delete[] node_que;

    reset_temp_vars(s, reversed);
};


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
                            // count.shrink_to_fit();
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

inline void IncrementalTopK::verify_sizes(){
    
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


void IncrementalTopK::update_loops() {
    
    this->vertices_to_update.clear();
    std::set<vertex> reset_visited;
    std::queue<vertex> *q = new std::queue<vertex>();

    vertex dequeued_v;
    vertex to_v;
    q->push(ordering[this->x]);
    this->visited_in_update_loops[ordering[this->x]] = 0;
    
    while(!q->empty()){
        dequeued_v = q->front();
        q->pop();
        if(this->visited_in_update_loops[dequeued_v] > K){
            continue;
        }
        this->vertices_to_update.insert(dequeued_v);
        
        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){
            to_v = ordering[u];
            if(this->visited_in_update_loops[to_v] == null_distance && to_v != ordering[this->y]){
                q->push(to_v);
                this->visited_in_update_loops[to_v] = this->visited_in_update_loops[dequeued_v] + 1;
                if(reset_visited.size() != graph->numberOfNodes()){
                    reset_visited.insert(to_v);
                }
            }
        }
    }
    assert(q->empty());
    
    
    
    q->push(ordering[this->y]);

    this->visited_in_update_loops[ordering[this->y]] = 0;

    while(!q->empty()){
        dequeued_v = q->front(); 
        q->pop();
        if(this->visited_in_update_loops[dequeued_v] > K){
            continue;
        }
        this->vertices_to_update.insert(dequeued_v);
        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){

            to_v = ordering[u];
            if(this->visited_in_update_loops[to_v] == null_distance){
                q->push(to_v);
                this->visited_in_update_loops[to_v] = this->visited_in_update_loops[dequeued_v] + 1;
                if(reset_visited.size() != graph->numberOfNodes()){
                    reset_visited.insert(to_v);
                }
            }
        }
    }
    assert(!graph->hasEdge(this->x,this->y));
    delete q;

    graph->addEdge(this->x,this->y);

    
    aff_cycles = this->vertices_to_update.size();
    
    this->reached_mbfs.clear(); // tracing nodes visited while updating loops

    
    vertex ordered_degree = 0;
    vertex u;
    std::set<vertex>::iterator it;

    for(it=this->vertices_to_update.begin();it!=this->vertices_to_update.end();it++){
    // for(vertex u: to_update){
        u = *it;
        if(u > std::min(ordering[this->x], ordering[this->y])){
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
        loop_entries -= this->loop_labels[u].size();
        this->loop_labels[u].clear();
        // this->loop_labels[u].shrink_to_fit();
        this->compute_loop_entries(u);
    }
    for(it=reset_visited.begin();it!=reset_visited.end();it++){
        this->visited_in_update_loops[*it] = null_distance;
    }
    reset_visited.clear();
    
    this->visited_in_update_loops[ordering[this->x]] = null_distance;
    this->visited_in_update_loops[ordering[this->y]] = null_distance;
    #ifndef NDEBUG
    this->graph->parallelForNodes([&] (vertex i){
            assert(graph->hasNode(i));
            assert(this->visited_in_update_loops[i]==null_distance);
        });
    #endif
}

void IncrementalTopK::update_lengths() {
    this->aff_hubs = 0;
    this->reached_nodes.clear();

    const index_t &idva = length_labels[0][ordering[this->x]];
    const index_t &idvb = length_labels[0][ordering[this->y]];

    
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
    // this->old_distances_a.shrink_to_fit();
    // this->old_distances_b.shrink_to_fit();
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
        
        
        assert(this->new_labels.empty());
        
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
                assert(this->updated.empty());
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
                assert(this->updated.empty());
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
                assert(this->updated.empty());
                reset_temp_vars(w_a, directed);
            }
            if(w_b < ordering[this->x]){
                for(size_t i = 0; i < this->old_distances_b[pos_b].size(); i++){
                    for(size_t j = 0; j < this->old_distances_b[pos_b][i]; j++){
                        resume_pbfs(w_b,ordering[this->x], i+this->old_label_b[pos_b].second+1, false);
                    }
                }
                assert(this->updated.empty());
                reset_temp_vars(w_b, directed);
            }
            pos_a++;
            pos_b++;
        }
        
        while(!this->new_labels.empty()){
            extend_label_repair(std::get<0>(this->new_labels.front()), std::get<1>(this->new_labels.front()), std::get<2>(this->new_labels.front()), std::get<3>(this->new_labels.front()), std::get<4>(this->new_labels.front()));
            this->new_labels.pop();
        }
    }
    // std::cout<<"done!\n";
}

inline bool IncrementalTopK::prune(vertex v,  dist d, bool rev){
    const index_t &idv = length_labels[rev][v];

    size_t pcount = 0;
    vertex w = 0;
    // cerr << "prune start" << endl;
    for (size_t pos = 0; pos < idv.label_offset.size(); pos++){
        w = idv.label_offset[pos].first;

        if (tmp_s_offset[w] == null_distance) continue;

        const vector<dist> &dcs = tmp_s_count[w];

        int l = dcs.size() - 1;
        int c = d - tmp_s_offset[w] - idv.label_offset[pos].second;

        for (int i = 0; i <= c && i < idv.d_array[pos].size(); i++){
            pcount += (int)dcs[std::min(c - i, l)] * idv.d_array[pos][i];
        }

        if (pcount >= K) return true;
    }
    return false;
}

inline void IncrementalTopK::resume_pbfs(vertex s, vertex t, dist d, bool rev) {
    
    set_temp_vars(s, rev);

    vertex     curr = 0;
    vertex     next = 1;
    dist distance = d;

    this->node_que = new std::queue<vertex>[2];



    assert(this->updated.empty());
    this->node_que[curr].push(t);
    this->tmp_dist_count[curr][t] = 1;
    
    this->updated.push_back(t);
    
    vertex currently_reached_nodes = 0;
    this->dists.clear();
    // this->dists.shrink_to_fit();
    vertex to_v;
    vertex v;
    dist c;
    for (;;){

        while (!this->node_que[curr].empty()){

            v = this->node_que[curr].front(); 
            this->node_que[curr].pop();
            c = this->tmp_dist_count[curr][v];
            this->tmp_dist_count[curr][v] = 0;

            if(c == 0 || this->tmp_pruned[v]){
                continue;
            }
            this->query(reverse_ordering[s], reverse_ordering[v], dists);
            currently_reached_nodes += 1;
            this->tmp_pruned[v] = this->dists.size() == this->K && *dists.rbegin() <= distance;

            if(this->tmp_pruned[v]){
                continue;
            }

            this->new_labels.push(std::tuple(v, s, distance, c, rev, tmp_offset[v]));

            for(vertex u : graph->neighborRange(reverse_ordering[v])){

                to_v = ordering[u];
                if(this->tmp_count[to_v] == 0){
                    this->updated.push_back(to_v);
                }

                if(to_v > s && tmp_count[to_v] < K){
                    this->tmp_count[to_v] += c;
                    this->node_que[next].push(to_v);
                    this->tmp_dist_count[next][to_v] += c;
                }
            }
        }

        if (this->node_que[next].empty()){
            break;
        }
        swap(curr, next);
        distance++;
    }
    delete[] this->node_que;
    this->reached_nodes.push_back(currently_reached_nodes);
    
    reset_temp_vars(s, rev);
}

inline void IncrementalTopK::set_temp_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];

    vertex w;

    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        
        this->tmp_s_offset[w] = ids.label_offset[pos].second;

        this->tmp_v.clear();

        for(size_t i = 0; i < ids.d_array[pos].size(); i++){
            this->tmp_v.push_back(ids.d_array[pos][i]);
        }
        
        this->tmp_s_count[w].resize(tmp_v.size() + loop_labels[w].size() - 1, 0);
        // this->tmp_s_count[w].shrink_to_fit();
        for(size_t i = 0; i < tmp_v.size(); i++){
            for(size_t j = 0; j < loop_labels[w].size(); j++){
                this->tmp_s_count[w][i+j] += this->tmp_v[i] * this->loop_labels[w][j];
            }
        }
    }
}

inline void IncrementalTopK::reset_temp_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];
    vertex w;
    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        this->tmp_s_offset[w] = null_distance;
        this->tmp_s_count[w].clear();
        // this->tmp_s_count[w].shrink_to_fit();
    }

    for(std::vector<vertex>::iterator i=this->updated.begin(); i != this->updated.end(); i++){
        w = *i;
        this->tmp_count[w] = 0;
        this->tmp_offset[w] = null_distance;
        this->tmp_pruned[w] = false;
        for(size_t j = 0; j < 2; j++){
            this->tmp_dist_count[j][w] = 0;
        }
    }

    this->updated.clear();

    // #ifndef NDEBUG
    //     this->graph->parallelForNodes([&] (vertex i){
    //         assert(graph->hasNode(i));
    //         assert(this->tmp_pruned[i]==false);
    //         assert(this->tmp_offset[i]==null_distance);
    //         assert(this->tmp_count[i]==0);
    //         assert(this->tmp_s_offset[i]==null_distance);

    //         assert(this->tmp_s_count[i].empty());

    //         assert(this->visited_in_update_loops[i]==null_distance);
    //         for (size_t j = 0; j < 2; j++)
    //             assert(this->tmp_dist_count[j][i]==0);
    
    //     });
    //     assert(this->tmp_s_offset[graph->numberOfNodes()]==0);
    // #endif
}

inline void IncrementalTopK::allocate_label(vertex v, vertex start, dist distance, dist count, bool dir){
    index_t &idv = length_labels[dir][v];
    
    idv.label_offset.emplace_back(start,distance);
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

 double IncrementalTopK::n_reached_nodes_mbfs(){
     double sum = 0.0;
     for(auto& element:reached_mbfs){
         sum+=(double)element;
     }
     return reached_mbfs.size()>0 ? sum / (double) reached_mbfs.size() : 0.0;
 }

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
        total = 0;
        for(size_t i = 0; i < idv.d_array[p].size(); i++){
            total += idv.d_array[p][i];
            if(total >= K){
                ol_size = idv.d_array[p].size();
                idv.d_array[p].resize(i+1);
                // idv.d_array[p].shrink_to_fit();

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
            // idv.d_array[last].shrink_to_fit();

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
            total = 0;
            for(size_t i = 0; i < idv.d_array[p].size(); i++){
                total += idv.d_array[p][i];
                if(total >= K){
                    idv.d_array[p].resize(i+1);
                    // idv.d_array[p].shrink_to_fit();
                    break;
                }
            }
        }
    };

}

// void IncrementalTopK::mod_bfs(vertex s, vertex t, std::vector<dist> &ret) {
//     s = ordering[s];
//     t = ordering[t];
//     std:vector<std::vector<dist> > distanze(graph->numberOfNodes());
//     std::priority_queue<pair<dist, vertex > > que;

//     que.push(make_pair(0, s));
//     vertex v;
//     vertex dist_negative;
//     while (!que.empty()) {
//         v = que.top().second;
//         dist_negative = -que.top().first;
//         que.pop();

//         if(distanze[t].size() >= this->K){
//             break;
//         }
//         if (distanze[v].size() >= this->K){
//             continue;
//         }

//         distanze[v].push_back(dist_negative);
//         for(vertex u : graph->neighborRange(reverse_ordering[v])){

//         // graph->forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
            
//             que.push(make_pair(-(1 + dist_negative), ordering[u]));
//         }
//     }
//     ret = distanze[t];
// }

