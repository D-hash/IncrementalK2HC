//
// Created by andrea on 20/12/22.
//
#include "incremental_topk.h"
#include <queue>
#include <set>
#include <algorithm>
#include <cassert>
#include <climits>
using namespace std;

const uint8_t IncrementalTopK::INF8 = std::numeric_limits<uint8_t>::max() / 2;
template <typename T> inline bool ReAlloc(T*& ptr, size_t nmemb){
    ptr = (T*)realloc(ptr, nmemb * sizeof(T));
    return ptr != NULL;
}

double GetCurrentTimeSec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

size_t IncrementalTopK::
NumOfVertex()        { return V; }

bool IncrementalTopK::
ConstructIndex(const vector<pair<int, int> > &es, size_t K, bool directed){
    Free();

    this->V = 0;
    this->K = K;
    this->directed = directed;

    for (size_t i = 0; i < es.size(); i++){
        V = std::max(V, (size_t)std::max(es[i].first, es[i].second) + 1);
    }

    for (int dir = 0; dir < 1 + directed; dir++){
        graph[dir].resize(V);
    }

    // renaming
    {
        vector<std::pair<int, int> > deg(V, std::make_pair(0, 0));

        for (size_t i = 0; i < V; i++) deg[i].second = i;

        for (size_t i = 0; i < es.size(); i++){
            deg[es[i].first ].first++;
            deg[es[i].second].first++;
        }

        sort(deg.begin(), deg.end(), greater<pair<int, int> >());
        ordering.resize(V);
        reverse_ordering.resize(V);
        for (size_t i = 0; i < V; i++) {
            ordering[deg[i].second] = i;
            reverse_ordering[i] = deg[i].second;
        }

        for (size_t i = 0; i < es.size(); i++){
            graph[0][ordering[es[i].first]].push_back(ordering[es[i].second]);

            if (directed){
                graph[1][ordering[es[i].second]].push_back(ordering[es[i].first]);
            } else {
                graph[0][ordering[es[i].second]].push_back(ordering[es[i].first]);
            }
        }
    }

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
    // cerr << directed << " " << s << " " << t << endl;
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

void IncrementalTopK::
KHubsDistanceQuery(uint32_t s, uint32_t t, uint8_t k, vector<pair<u_int32_t,int>> &ret){
    ret.clear();
    size_t pos1 = 0;
    size_t pos2 = 0;

    // cerr << directed << " " << s << " " << t << endl;
    const index_t &ids = index[directed][s];
    const index_t &idt = index[0][t];
    vector<vector<u_int32_t>> count(30);
    for (;;){
        if (pos1 >= ids.label_offset.size()) break;
        if (pos2 >= idt.label_offset.size()) break;

        if (ids.label_offset[pos1].first == idt.label_offset[pos2].first){
            uint32_t W = ids.label_offset[pos1].first;

            for (size_t i = 0; i < ids.d_array[pos1].size(); i++){
                for (size_t j = 0; j < idt.d_array[pos2].size(); j++){
                    uint8_t d_tmp = ids.label_offset[pos1].second + idt.label_offset[pos2].second + i + j;
                    if (count.size() <= d_tmp) count.resize(d_tmp + 1);
                    for(int idx = 0; idx < (int)ids.d_array[pos1][i] *  idt.d_array[pos2][j]; idx++){
                        count[d_tmp].push_back(W);
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

    for (size_t i = 0; i < count.size(); i++) {
        for (size_t j = 0; j < count[i].size(); j++){
            ret.emplace_back(count[i][j], i);
            if (ret.size() == k) return;
        }
    }
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

    for (int i = 0; i < 2; i++){
        graph[i].clear();
    }

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
    const vector<vector<uint32_t> > &fgraph = graph[0];

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
            for (size_t i = 0; i < fgraph[v].size(); i++){
                uint32_t to = fgraph[v][i];

                if (tmp_count[to] == 0){
                    updated.push_back(to);
                }

                if (to >= s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            }
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
    const vector<vector<uint32_t> > &graph_ = graph[rev];

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
            tmp_pruned[v] = Pruning(v, dist, rev);
            // cerr << "Pruning done" << endl;

            if(tmp_pruned[v]) continue;

            if(tmp_offset[v] == INF8){
                // Make new label for a node v
                tmp_offset[v] = dist;
                AllocLabel(v, s, dist, c, rev);
            }else{
                // assert(s != 3);
                ExtendLabel(v, s, dist, c, rev, 0);
            }

            for(size_t i = 0; i < graph_[v].size(); i++){
                uint32_t to  = graph_[v][i];
                if(tmp_count[to] == 0){
                    updated.push_back(to);
                }

                if(to > s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            }
        }

        if (node_que[next].empty()) break;
        swap(curr, next);
        dist++;
    }
    // cerr <<"#visited nodes: " << num_of_labeled_vertices[s] << endl;
    ResetTempVars(s, updated, rev);
};

void IncrementalTopK::
UpdateLoops(std::pair<int, int> new_edge) {
    std::vector<bool> visited;
    visited.resize(V, false);
    std::set<uint32_t> to_update;
    std::queue<uint32_t> q;
    uint32_t a = new_edge.first;
    uint32_t b = new_edge.second;
    RemoveEdge(a,b);
    const vector<vector<uint32_t> > &graph_ = graph[0];
    q.push(ordering[a]);
    int distance = 0;
    while(!q.empty() and distance <= K){
        uint32_t vertex = q.front();
        q.pop();
        aff_cycles++;
        visited[vertex] = true;
        to_update.insert(vertex);
        for(size_t i = 0; i < graph_[vertex].size(); i++){
            uint32_t to  = graph_[vertex][i];
            if(visited[to]) continue;
            q.push(to);
        }
        distance ++;
    }
    std::queue<uint32_t> empty;
    std::swap( q, empty );
    q.push(ordering[b]);
    distance = 0;
    while(!q.empty() and distance <= K){
        uint32_t vertex = q.front(); q.pop();
        visited[vertex] = true;
        to_update.insert(vertex);
        for(size_t i = 0; i < graph_[vertex].size(); i++){
            uint32_t to  = graph_[vertex][i];
            if(visited[to]) continue;
            q.push(to);
        }
        distance ++;
    }
    AddEdge(a,b);
    uint32_t min_order = min(ordering[a], ordering[b]);
    aff_cycles = to_update.size();
    reached_mbfs.clear();
    for(uint32_t u: to_update){
        if(ordering[u] > min_order) continue;
        loop_count[u].clear();
        bool status = true;
        CountLoops(u, status);
    }
}

void IncrementalTopK::
UpdateIndex(std::pair<int, int> new_edge) {
    uint32_t a = new_edge.first;
    uint32_t b = new_edge.second;

    AddEdge(a,b);

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
    while (pos_a != idva.label_offset.size() && pos_b != idvb.label_offset.size()){
        vector<tuple<u_int32_t, u_int32_t, uint8_t, u_int8_t, bool, u_int32_t>> new_labels;
        uint32_t w_a = pos_a < idva.label_offset.size() ? idva.label_offset[pos_a].first : V;
        uint32_t w_b = pos_b < idvb.label_offset.size() ? idvb.label_offset[pos_b].first : V;
        if(w_a < w_b){
            if(w_a < ordering[b]){
                aff_hubs++;
                for(size_t i = 0; i < idva.d_array[pos_a].size(); i++)
                    for(size_t j = 0; j < idva.d_array[pos_a][i]; j++)
                        ResumePBfs(w_a,ordering[b], i+idva.label_offset[pos_a].second+1, false, status, new_labels);
            }
            pos_a++;
        }
        else if(w_b < w_a){
            if(w_b < ordering[a]){
                aff_hubs++;
                for(size_t i = 0; i < idvb.d_array[pos_b].size(); i++)
                    for(size_t j = 0; j < idvb.d_array[pos_b][i]; j++)
                        ResumePBfs(w_b,ordering[a], i+idvb.label_offset[pos_b].second+1, false, status, new_labels);
            }
            pos_b++;
        }
        else {
            aff_hubs++;
            if(w_a < ordering[b])
                for(size_t i = 0; i < idva.d_array[pos_a].size(); i++)
                    for(size_t j = 0; j < idva.d_array[pos_a][i]; j++)
                        ResumePBfs(w_a,ordering[b], i+idva.label_offset[pos_a].second+1, false, status, new_labels);
            if(w_b < ordering[a])
                for(size_t i = 0; i < idvb.d_array[pos_b].size(); i++)
                    for(size_t j = 0; j < idvb.d_array[pos_b][i]; j++)
                        ResumePBfs(w_b,ordering[a], i+idvb.label_offset[pos_b].second+1, false, status, new_labels);
            pos_a++; pos_b++;
        }
        for(auto nl: new_labels)
            ExtendLabelRepair(std::get<0>(nl), std::get<1>(nl), std::get<2>(nl), std::get<3>(nl), std::get<4>(nl));
        ++up_bar;
    }
}

void IncrementalTopK::
AddEdge(uint32_t a, uint32_t b){
    graph[0][ordering[a]].push_back(ordering[b]);
    graph[0][ordering[b]].push_back(ordering[a]);
}

void IncrementalTopK::
RemoveEdge(uint32_t a, uint32_t b){
    graph[0][ordering[a]].pop_back();
    graph[0][ordering[b]].pop_back();
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
    const vector<vector<uint32_t> > &graph_ = graph[dir];
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
            // cerr << "Pruning done" << endl;

            if(tmp_pruned[v]) continue;

            new_labels.emplace_back(v, s, dist, c, dir, tmp_offset[v]);

            for(size_t i = 0; i < graph_[v].size(); i++){
                uint32_t to  = graph_[v][i];
                if(tmp_count[to] == 0){
                    updated.push_back(to);
                }

                if(to > s && tmp_count[to] < K){
                    tmp_count[to] += c;
                    node_que[next].push(to);
                    tmp_dist_count[next][to] += c;
                }
            }
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
    // cerr << rev << " " << s << " " << V << endl;
    const index_t &ids = index[directed && !rev][s];

    // cerr << ids.length << " " << ids.label[0] << endl;
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

inline bool IncrementalTopK::
Pruning(uint32_t v,  uint8_t d, bool rev){
    const index_t &idv = index[rev][v];

    size_t pcount = 0;

    // cerr << "Pruning start" << endl;
    for (size_t pos = 0; pos < idv.label_offset.size(); pos++){
        uint32_t w = idv.label_offset[pos].first;

        if (tmp_s_offset[w] == INF8) continue;

        const vector<uint8_t> &dcs = tmp_s_count[w];

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

