//
// Created by andrea on 20/12/22.
//
#include "incremental_topk.h"
#include <queue>
#include <set>
#include <algorithm>
#include <cassert>
#include <climits>

#include "networkit/centrality/DegreeCentrality.hpp"

using namespace std;

const uint8_t IncrementalTopK::INF8 = std::numeric_limits<uint8_t>::max() / 2;


double GetCurrentTimeSec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

size_t IncrementalTopK::
NumOfVertex()        { return V; }

bool IncrementalTopK::
ConstructIndex(NetworKit::Graph graph, size_t K, bool directed){
    Free();

    this->graph = graph;
    this->V = graph.numberOfNodes();
    auto deg = new NetworKit::DegreeCentrality(graph);
    deg->run();
    ordering.resize(V);
    reverse_ordering.resize(V);
    auto rank = deg->ranking();
    for(uint32_t s = 0; s < V; s++){
        ordering[rank[s].first] = s;
        reverse_ordering[s] = rank[s].first;
        assert(rank[s].second == graph.degree(rank[s].first));
    }

    this->K = K;
    this->directed = directed;

    Init();

    return Labeling();
}

IncrementalTopK::
~IncrementalTopK(){
    Free();
}

int IncrementalTopK::
KDistanceQuery(int s, int t, uint8_t k){
    vector<int> dists;
    return KDistanceQuery(s, t, k, dists);
}

int IncrementalTopK::
KDistanceQuery(int s, int t, uint8_t k, vector<int> &ret){
    ret.clear();
    s = ordering[s];
    t = ordering[t];
    size_t pos1 = 0;
    size_t pos2 = 0;

    vector<int> count(30, 0);
    const index_t &ids = index[directed][s];
    const index_t &idt = index[0][t];

    for (;;){
        if (pos1 >= ids.label_offset.size()) break;
        if (pos2 >= idt.label_offset.size()) break;
        if (ids.label_offset[pos1].first == idt.label_offset[pos2].first){
            uint32_t W = ids.label_offset[pos1].first;

            for (size_t i = 0; i < ids.d_array[pos1].size(); i++){
                for (size_t j = 0; j < idt.d_array[pos2].size(); j++){
                    for (size_t m = 0; m < loop_count[W].size(); m++){
                        uint8_t d_tmp = ids.label_offset[pos1].second + idt.label_offset[pos2].second + i + j + m;
                        uint8_t c_tmp = loop_count[W][m] - (m ? loop_count[W][m-1] : 0);
                        if (count.size() <= d_tmp) count.resize(d_tmp + 1, 0);
                        count[d_tmp] += (int)ids.d_array[pos1][i] * idt.d_array[pos2][j] * c_tmp;
                    }
                }
            }
            pos1++, pos2++;
        } else {
            if (ids.label_offset[pos1].first < idt.label_offset[pos2].first){
                pos1++;
            } else {
                pos2++;
            }
        }
    }

    for (size_t i = 0; i < count.size(); i++){
        while (ret.size() < k && count[i]-- > 0){
            ret.push_back(i);
        }
    }

    return ret.size() < k ? INT_MAX : 0;
}

size_t IncrementalTopK::
IndexSize(){
    size_t sz = 0;

    for(size_t v = 0; v < V; v++){
        sz += sizeof(uint8_t ) * loop_count[v].size(); // loopcount
    }

    // index's size
    for (int dir = 0; dir < 1 + directed; dir++){
        for(size_t v = 0; v < V; v++){
            sz += index[dir][v].label_offset.size();
            for (size_t a = 0; a < index[dir][v].d_array.size(); a++)
                sz += index[dir][v].d_array[a].size();
        }
    }
    return sz;
}

double IncrementalTopK::
AvgIndexSize() {
    size_t sz = 0;

    for (int dir = 0; dir < 1 + directed; dir++){
        for(size_t v = 0; v < V; v++){
            sz += index[dir][v].label_offset.size();
            for (size_t a = 0; a < index[dir][v].d_array.size(); a++)
                sz += index[dir][v].d_array[a].size();
        }
    }
    return (double)((double)sz/(double)V);
}

void IncrementalTopK::
Init(){
    tmp_pruned.resize(V, false);
    tmp_offset.resize(V, INF8);
    tmp_count .resize(V, 0);
    tmp_s_offset.resize(V, INF8); tmp_s_offset.push_back(0);
    tmp_s_count .resize(V);

    for (int j = 0; j < 2; j++) tmp_dist_count[j].resize(V, 0);

    loop_count.resize(V);

    for (int dir = 0; dir < 1 + directed; dir++){
        index[dir].resize(V);

        for (size_t v = 0; v < V; v++){
            index[dir][v].label_offset.clear();
            index[dir][v].d_array.clear();
            index[dir][v].d_array.clear();
        }
    }
    reached_mbfs.clear();
}

void IncrementalTopK::
Free(){
    V = K = loop_count_time = indexing_time = 0;
    directed = false;

    ordering.clear();
    reverse_ordering.clear();
    loop_count.clear();

    tmp_pruned.clear();
    tmp_offset.clear();
    tmp_count .clear();

    for (int i = 0; i < 2; i++){
        tmp_dist_count[i].clear();
    }


    for (int dir = 0; dir < 1 + directed; dir++){
        for (size_t v = 0; v < V; v++){
            index_t &idv = index[dir][v];

            idv.label_offset.clear();

            if (!idv.d_array.empty()){
                for (size_t i = 0; i < idv.label_offset.size(); i ++){
                    idv.d_array.clear();
                }
            }
        }
        index[dir].clear();
    }
}

bool IncrementalTopK::
Labeling(){

    bool status = true;

    loop_count_time = -GetCurrentTimeSec();
    ProgressStream loop_bar(V);
    loop_bar.label() << "Loops construction";
    for(size_t v = 0; v < V; v++){
        CountLoops(v, status);
        ++loop_bar;
    }
    loop_count_time += GetCurrentTimeSec();

    indexing_time = -GetCurrentTimeSec();
    ProgressStream index_bar(V);
    index_bar.label() << "Index construction";
    for(size_t v = 0; v < V; v++){

        // compute L_in
        PrunedBfs(v, false, status);
        ++index_bar;
        if (directed){
            //  compute L_out
            PrunedBfs(v, true, status);
        }
    }
    indexing_time += GetCurrentTimeSec();

    return status;
}

void IncrementalTopK::
CountLoops(uint32_t s, bool &status){
    size_t  count = 0;
    int     curr  = 0;
    int     next  = 1;
    uint8_t dist  = 0;

    std::queue<uint32_t> node_que[2];
    vector<uint32_t>     updated;

    node_que[curr].push(s);
    updated.push_back(s);
    tmp_dist_count[curr][s] = 1;
    uint32_t currently_reached_nodes = 0;
    for (;;){
        if (dist == INF8 && status){
            cerr << "Warning: Self loops become too long." << endl;
            status = false;
        }

        while (!node_que[curr].empty() && count < K){
            uint32_t v = node_que[curr].front(); node_que[curr].pop();
            uint8_t  c = tmp_dist_count[curr][v]; // the number of path from s to v with dist hops.
            tmp_dist_count[curr][v] = 0;
            if (c == 0) continue;

            if (v == s){
                loop_count[s].resize(dist + 1, 0);
                loop_count[s][dist] += c;
                count += c;
            }
            currently_reached_nodes++;
            graph.forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
                uint32_t to = ordering[u];
                if (tmp_count[to] == 0){
                    updated.push_back(to);
                }

                if (to >= s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            });
        }
        if(node_que[next].empty() || count >= K) break;
        swap(curr, next);
        dist++;
    }

    for(size_t i = 1; i < loop_count[s].size(); i++){
        loop_count[s][i] += loop_count[s][i-1];
    }
    assert(loop_count[s][0] == 1);
    reached_mbfs.push_back(currently_reached_nodes);
    ResetTempVars(s, updated, false);
}

void IncrementalTopK::
PrunedBfs(uint32_t s, bool rev, bool &status){
    SetStartTempVars(s, rev);

    int     curr = 0;
    int     next = 1;
    uint8_t dist = 0;

    std::queue<uint32_t> node_que[2];
    vector<uint32_t>     updated;

    node_que[curr].push(s);
    tmp_dist_count[curr][s] = 1;
    updated.push_back(s);

    for (;;){
        if (dist == INF8 && status){
            cerr << "Warning: Distance from a source node becomes too long." << endl;
            status = false;
        }

        while (!node_que[curr].empty()){

            uint32_t v = node_que[curr].front(); node_que[curr].pop();
            uint8_t  c = tmp_dist_count[curr][v];
            tmp_dist_count[curr][v] = 0;

            if(c == 0 || tmp_pruned[v]) continue;
            vector<int> ret;
            int check = KDistanceQuery(reverse_ordering[s], reverse_ordering[v], K, ret);
            tmp_pruned[v] = (check == 0 && ret.back() <= dist);

            if(tmp_pruned[v]) continue;

            if(tmp_offset[v] == INF8){
                // Make new label for a node v
                tmp_offset[v] = dist;
                AllocLabel(v, s, dist, c, rev);
            }else{
                ExtendLabel(v, s, dist, c, rev, 0);
            }
            graph.forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
                uint32_t to = ordering[u];
                if(tmp_count[to] == 0){
                    updated.push_back(to);
                }

                if(to > s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            });
        }

        if (node_que[next].empty()) break;
        swap(curr, next);
        dist++;
    }
    ResetTempVars(s, updated, rev);
};

void IncrementalTopK::
UpdateLoops(std::pair<int, int> new_edge) {
    std::vector<uint8_t> visited;
    visited.resize(V, INF8);
    std::set<uint32_t> to_update;
    std::queue<uint32_t> q;
    uint32_t a = new_edge.first;
    uint32_t b = new_edge.second;
    // RemoveEdge(a,b);
    q.push(ordering[a]);
    visited[ordering[a]] = 0;
    while(!q.empty()){
        uint32_t vertex = q.front();
        q.pop();
        if(visited[vertex] > K) continue;
        aff_cycles++;
        to_update.insert(vertex);
        graph.forNeighborsOf(reverse_ordering[vertex], [&](NetworKit::node u) {
            uint32_t to = ordering[u];
            if(visited[to] == INF8 && to != ordering[b]){
                q.push(to);
                visited[to] = visited[vertex] + 1;
            }
        });
    }
    std::queue<uint32_t> empty;
    std::swap( q, empty );
    q.push(ordering[b]);
    visited[ordering[b]] = 0;
    while(!q.empty()){
        uint32_t vertex = q.front(); q.pop();
        if(visited[vertex] > K) continue;
        aff_cycles++;
        to_update.insert(vertex);
        graph.forNeighborsOf(reverse_ordering[vertex], [&](NetworKit::node u) {
            uint32_t to = ordering[u];
            if(visited[to] == INF8){
                q.push(to);
                visited[to] = visited[vertex] + 1;
            }
        });
    }
    AddEdge(a,b);
    uint32_t min_order = min(ordering[a], ordering[b]);
    aff_cycles = to_update.size();
    reached_mbfs.clear();
    for(uint32_t u: to_update){
        if(u > min_order) continue;
        loop_count[u].clear();
        bool status = true;
        CountLoops(u, status);
    }
}

void IncrementalTopK::
UpdateIndex(std::pair<int, int> new_edge) {
    uint32_t a = new_edge.first;
    uint32_t b = new_edge.second;

    // AddEdge(a,b);

    aff_hubs = 0;
    reached_nodes.clear();
    bool status = true;
    const index_t &idva = index[0][ordering[a]];
    const index_t &idvb = index[0][ordering[b]];
    size_t max_length = max(idva.label_offset.size(),idvb.label_offset.size());
    size_t pos_a = 0;
    size_t pos_b = 0;
    ProgressStream up_bar(max_length);
    up_bar.label() << "Updating affected hubs..";
    vector<pair<uint32_t,uint8_t>> old_label_a;
    vector<pair<uint32_t,uint8_t>> old_label_b;
    vector<std::vector<uint8_t>> old_distances_a;
    vector<std::vector<uint8_t>> old_distances_b;
    for(auto elem: idva.label_offset) old_label_a.push_back(elem);
    for(auto elem: idvb.label_offset) old_label_b.push_back(elem);
    old_distances_a.resize(idva.d_array.size());
    old_distances_b.resize(idvb.d_array.size());
    for(size_t i = 0; i < idva.d_array.size(); i++){
        for(size_t j = 0; j < idva.d_array[i].size(); j++)
            old_distances_a[i].push_back(idva.d_array[i][j]);
    }
    for(size_t i = 0; i < idvb.d_array.size(); i++){
        for(size_t j = 0; j < idvb.d_array[i].size(); j++)
            old_distances_b[i].push_back(idvb.d_array[i][j]);
    }
    while (pos_a != old_label_a.size() && pos_b != old_label_b.size()){
        vector<tuple<u_int32_t, u_int32_t, uint8_t, u_int8_t, bool, u_int32_t>> new_labels;
        new_labels.clear();
        uint32_t w_a = pos_a < old_label_a.size() ? old_label_a[pos_a].first : V;
        uint32_t w_b = pos_b < old_label_b.size() ? old_label_b[pos_b].first : V;
        if(w_a < w_b){
            if(w_a < ordering[b]){
                aff_hubs++;
                for(size_t i = 0; i < old_distances_a[pos_a].size(); i++)
                    for(size_t j = 0; j < old_distances_a[pos_a][i]; j++)
                        ResumePBfs(w_a,ordering[b], i+old_label_a[pos_a].second+1, false, status, new_labels);
            }
            pos_a++;
        }
        else if(w_b < w_a){
            if(w_b < ordering[a]){
                aff_hubs++;
                for(size_t i = 0; i < old_distances_b[pos_b].size(); i++)
                    for(size_t j = 0; j < old_distances_b[pos_b][i]; j++)
                        ResumePBfs(w_b,ordering[a], i+old_label_b[pos_b].second+1, false, status, new_labels);
            }
            pos_b++;
        }
        else {
            aff_hubs++;
            if(w_a < ordering[b])
                for(size_t i = 0; i < old_distances_a[pos_a].size(); i++)
                    for(size_t j = 0; j < old_distances_a[pos_a][i]; j++)
                        ResumePBfs(w_a,ordering[b], i+old_label_a[pos_a].second+1, false, status, new_labels);
            if(w_b < ordering[a])
                for(size_t i = 0; i < old_distances_b[pos_b].size(); i++)
                    for(size_t j = 0; j < old_distances_b[pos_b][i]; j++)
                        ResumePBfs(w_b,ordering[a], i+old_label_b[pos_b].second+1, false, status, new_labels);
            pos_a++; pos_b++;
        }
        for(auto nl: new_labels)
            ExtendLabelRepair(std::get<0>(nl), std::get<1>(nl), std::get<2>(nl), std::get<3>(nl), std::get<4>(nl));
        ++up_bar;
    }
}

void IncrementalTopK::
AddEdge(uint32_t a, uint32_t b){
    assert(!graph.hasEdge(a,b));
    graph.addEdge(a,b);
}

void IncrementalTopK::
RemoveEdge(uint32_t a, uint32_t b){
    graph.removeEdge(a,b);
}

void IncrementalTopK::
ResumePBfs(uint32_t s, uint32_t t, uint8_t d, bool dir, bool &status,
           std::vector<std::tuple<u_int32_t, u_int32_t, uint8_t, u_int8_t, bool, u_int32_t>> &new_labels) {
    SetStartTempVars(s, dir);

    int     curr = 0;
    int     next = 1;
    uint8_t dist = d;

    std::queue<uint32_t> node_que[2];
    vector<uint32_t>     updated;

    node_que[curr].push(t);
    tmp_dist_count[curr][t] = 1;
    updated.push_back(t);
    uint32_t currently_reached_nodes = 0;
    for (;;){
        if (dist == INF8 && status){
            cerr << "Warning: Distance from a source node becomes too long." << endl;
            status = false;
        }


        while (!node_que[curr].empty()){

            uint32_t v = node_que[curr].front(); node_que[curr].pop();
            uint8_t  c = tmp_dist_count[curr][v];
            tmp_dist_count[curr][v] = 0;

            if(c == 0 || tmp_pruned[v]) continue;
            vector<int> dists;
            KDistanceQuery(reverse_ordering[s], reverse_ordering[v], dists);
            currently_reached_nodes += 1;
            tmp_pruned[v] = dists.size() == K && *dists.rbegin() <= dist;

            if(tmp_pruned[v]) continue;

            new_labels.emplace_back(v, s, dist, c, dir, tmp_offset[v]);
            graph.forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
                uint32_t to = ordering[u];
                if(tmp_count[to] == 0){
                    updated.push_back(to);
                }

                if(to > s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            });
        }

        if (node_que[next].empty()) break;
        swap(curr, next);
        dist++;
    }
    reached_nodes.push_back(currently_reached_nodes);
    ResetTempVars(t, updated, dir);
}

inline void IncrementalTopK::
SetStartTempVars(uint32_t s, bool rev){
    const index_t &ids = index[directed && !rev][s];

    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        int w = ids.label_offset[pos].first;
        tmp_s_offset[w] = ids.label_offset[pos].second;

        vector<uint8_t> tmp_v;
        for(size_t i = 0; i < ids.d_array[pos].size(); i++){
            tmp_v.push_back(ids.d_array[pos][i]);
        }
        tmp_s_count[w].resize(tmp_v.size() + loop_count[w].size() - 1, 0);

        for(size_t i = 0; i < tmp_v.size(); i++){
            for(size_t j = 0; j < loop_count[w].size(); j++){
                tmp_s_count[w][i+j] += tmp_v[i] * loop_count[w][j];
            }
        }
    }
}

inline void IncrementalTopK::
ResetTempVars(uint32_t s, const vector<uint32_t> &updated, bool rev){
    const index_t &ids = index[directed && !rev][s];

    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        int w = ids.label_offset[pos].first;
        tmp_s_offset[w] = INF8;
        tmp_s_count[w].clear();
    }

    for(size_t i = 0; i < updated.size(); i++){
        tmp_count [updated[i]] = 0;
        tmp_offset[updated[i]] = INF8;
        tmp_pruned[updated[i]] = false;
        for(int j = 0; j < 2; j++) tmp_dist_count[j][updated[i]] = 0;
    }
}

inline void IncrementalTopK::
AllocLabel(uint32_t v, uint32_t start, uint8_t dist, uint8_t count, bool dir){
    index_t &idv = index[dir][v];

    idv.label_offset.push_back(make_pair(start,dist));
    idv.d_array.resize(idv.d_array.size()+1);
    idv.d_array[idv.label_offset.size()-1].resize(1,count);

}

inline void IncrementalTopK::
ExtendLabel(uint32_t v, uint32_t start, uint8_t dist, uint8_t count, bool dir, size_t pos){
    index_t &idv = index[dir][v];
    if(pos == 0)
        for(; pos < idv.label_offset.size(); pos++)
            if (idv.label_offset[pos].first == start) break;
    uint8_t offset = dist - idv.label_offset[pos].second;
    if (idv.d_array[pos].size() > offset)
        idv.d_array[pos][offset] += count;
    else{
        while(idv.d_array[pos].size() != offset) idv.d_array[pos].push_back(0);
        idv.d_array[pos].push_back(count);
    }
    for(size_t p = 0; p < idv.label_offset.size(); p++){
        int tot_count = 0;
        for(size_t i = 0; i < idv.d_array[p].size(); i++){
            tot_count += idv.d_array[p][i];
            if(tot_count >= K){
                idv.d_array[p].resize(i+1);
                break;
            }
        }
    }
}

inline void IncrementalTopK::
ExtendLabelRepair(uint32_t v, uint32_t start, uint8_t dist, uint8_t count, bool dir){
    index_t &idv = index[dir][v];

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
        AllocLabel(v, start, dist, count, dir);
        return;
    }
    else if(idv.label_offset[last].first == start){
        if(idv.label_offset[last].second > dist){
            idv.d_array[last].resize(idv.d_array[last].size() + idv.label_offset[last].second-dist, 0); // todo consider dequeue
            for(size_t rev = idv.d_array[last].size()-1; rev > 0; rev--)
                idv.d_array[last][rev] = idv.d_array[last][rev-(idv.label_offset[last].second-dist)];
            for(size_t rev = 1; rev != idv.label_offset[last].second-dist; rev++)
                idv.d_array[last][rev] = 0;
            idv.d_array[last][0] = count;
            idv.label_offset[last].second = dist;
        } else{
            ExtendLabel(v, start, dist, count, dir, last);
        }
    } else {
        idv.label_offset.insert(idv.label_offset.begin()+last, std::make_pair(start,dist));
        idv.d_array.insert(idv.d_array.begin() + last, {count});
        for(size_t p = 0; p < idv.label_offset.size(); p++){
            int tot_count = 0;
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

void IncrementalTopK::modBFS(uint32_t s, uint32_t t, std::vector<int> &ret) {
    s = ordering[s];
    t = ordering[t];
    vector<vector<int> > dist(V);
    priority_queue<pair<int, int > > que;
    que.push(make_pair(0, s));

    while (!que.empty()) {
        int v = que.top().second;
        int c = -que.top().first;
        que.pop();

        if(dist[t].size() >= K) break;
        if (dist[v].size() >= K)  continue;

        dist[v].push_back(c);
        graph.forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
            uint32_t to = ordering[u];
            que.push(make_pair(-(1 + c), to));
        });
    }
    ret = dist[t];
}

