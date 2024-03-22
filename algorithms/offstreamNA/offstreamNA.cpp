#include "offstreamNA.hpp"

#include <utility>

using namespace std;

//固定随机数
// 构造函数
OffstreamNAPartitioner::OffstreamNAPartitioner(BaseGraph& baseGraph, string input, const string &algorithm,
                             size_t num_partitions)
        : EdgePartitioner(baseGraph, algorithm, num_partitions), input(std::move(input)), gen(985) {
    config_output_files();
    current_partition = 0;
    assigned_edges = 0;
    stream_part = baseGraph.stream_part;
    off_part = baseGraph.off_part;
    partial_degree = baseGraph.partial_degree;
    // capacity = num_edges * BALANCE_RATIO / num_partitions + 1;
    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundaries.assign(num_partitions, dense_bitset(num_vertices));
    dis.param(
            uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    window.reserve(WINDOW_SIZE);

    vertex_partitions.assign(num_vertices, vector<size_t>(num_partitions, 0));
    vp_set.assign(num_vertices, dense_bitset(num_partitions));

    capacity =  BALANCE_RATIO * (off_part.size() / (double)num_partitions) + 1;
    max_partition_load = (uint64_t) BALANCE_RATIO * num_edges / num_partitions;


}

//最后一个子图就是剩下边组合而成
void OffstreamNAPartitioner::assign_remaining() {
    auto &is_boundary = is_boundaries[num_partitions - 1], &is_core = is_cores[num_partitions - 1];
    repv(u, num_vertices) for (auto &i: adj_out[u])
            if (edges[i.v].valid()) {
                vid_t v = edges[i.v].second;
                assign_edge(num_partitions - 1, u, v);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(v);
            }
    repv(i, num_vertices) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            //在其他子图中是核心集的不予理睬，不能设置为本子图的核心集
            rep(j, num_partitions - 1) if (is_cores[j].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
        }
    }
}

size_t OffstreamNAPartitioner::count_mirrors() {
    size_t result = 0;
    rep(i, num_partitions) result += is_boundaries[i].popcount();
    return result;
}

void OffstreamNAPartitioner::split() {
    total_time.start();
    // 初始化最小堆，用于存储S\C的顶点信息
    min_heap.reserve(num_vertices);

    LOG(INFO) << "Start offstreamNA partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    size_t off_size = off_part.size();
    size_t on_size = stream_part.size();
    ss << "offstreamNA" << endl
       << "BALANCE RATIO: " << BALANCE_RATIO
       << " | Edges Ratio: " << EDGE_RATIO
       << " | Buffer Windows Size: " << WINDOW_SIZE
       << " | Offline edge size: " << off_size
       << " | Stream edge size: " << on_size
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());
    // TODO 先内存划分，再流式划分
    // 前p-1个分区
    Timer offline_time;
    offline_time.start();
    for (current_partition = 0; current_partition < num_partitions - 1; current_partition++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid;
            if (!min_heap.get_min(degree, vid)) { // 当S\C为空时，从V\C中随机选择顶点
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
//                     LOG(INFO) << "No vertex in V, current_partition: " << current_partition;
//                        repv(v, num_vertices) {
//                            if ((adj_out[v].size() + adj_in[v].size()) > 0 && !is_cores[current_partition].get(v)) {
//                                LOG(INFO) << v;
//                                LOG(INFO) << adj_out[v].size() + adj_in[v].size();
//                            }
//                        }
                    break;
                }
                degree = adj_out[vid].size() + adj_in[vid].size();
            } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
                min_heap.remove(vid);
            }
            // 将顶点加入到核心集C
            occupy_vertex(vid, degree);
        }
        min_heap.clear();

        rep(direction, 2) repv(vid, num_vertices) {
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        i++;
                    } else {
                        swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
    }

//    rep(direction, 2) repv(vid, num_vertices) {
//        if (adj_out[vid].size() + adj_in[vid].size()) {
//            LOG(INFO) << vid;
//            LOG(INFO) << adj_out[vid].size() + adj_in[vid].size();
//        }
//    }
    current_partition = num_partitions - 1;
    // 把剩余的边放入最后一个分区
    assign_remaining();
    offline_time.stop();
    repv(j, num_partitions) {
        LOG(INFO) << "Partition " << j << " Edge Count: " << occupied[j];
    }

    LOG(INFO) << "Offline partition cost time: " << offline_time.get_time();
    LOG(INFO) << "Allocate edges count: " << assigned_edges;
    LOG(INFO) << "Capacity: " << capacity;
    LOG(INFO) << "Capacity * partition: " << capacity * num_partitions;
    LOG(INFO) << "Window size: " << window.size();
    CHECK_EQ(assigned_edges + window.size(), off_size);

    // 使用贪心来划分

    min_load = *min_element(occupied.begin(), occupied.end());
    max_load = *max_element(occupied.begin(), occupied.end());

    LOG(INFO) << "Start adwise partitioning" << endl;
    LOG(INFO) << "max_partition_load: " << max_partition_load; // 50809
    // 下面进入adwise阶段

    Timer stream_time;
    stream_time.start();
    size_t i = 0;
    size_t remove_count = 0;
    // TODO 只移除了一半的边
    while(i < on_size) {
        auto& edge = stream_part[i++];
        if (window.size() == WINDOW_SIZE) {
            remove_edge_from_window();
        }
        add_edge_to_window(edge);
    }
    LOG(INFO) << "Window assign edges when stream end: " << remove_count;
    LOG(INFO) << "Remaining edges in window: " << window.size();
    while(!window.empty()) { // 每次移除一条边就重复计算剩余所有边
        remove_edge_from_window();
    }

//    for (auto &edge: stream_part) {
//        int partition = find_max_score_partition(edge);
//        is_mirrors[edge.first].set_bit_unsync(partition);
//        is_mirrors[edge.second].set_bit_unsync(partition);
//        occupied[partition]++;
//        assigned_edges++;
//        update_min_max_load(partition);
//    }
    stream_time.stop();



    LOG(INFO) << "Stream partition cost time: " << stream_time.get_time();

    CHECK_EQ(assigned_edges, num_edges);
    total_time.stop();
}


bool OffstreamNAPartitioner::get_free_vertex(vid_t &vid) {
    //随机选择一个节点
    size_t index = dis(gen);
    vid  = index;
    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices &&
           (adj_out[vid].size() + adj_in[vid].size() == 0
           || adj_out[vid].size() + adj_in[vid].size() > 2 * avg_degree
            || is_cores[current_partition].get(vid))) { // TODO 为什么要加这个判断条件
        vid = (index + ++count) % num_vertices;
    }
    if (count == num_vertices) {
        return false;
    }
//    bool re_get = false;
//    if (partial_degree[vid] == 1) {
//        // LOG(INFO) << "Vid size is 1: " << vid;
//        if (adj_out[vid].size() == 1) {
//            // 获取对端顶点
//            adjlist_t &neighbors = adj_out[vid];
//            vid_t &u = edges[neighbors[0].v].second;
//            if (partial_degree[u] == 1) {
//                LOG(INFO) << "U size is 1: " << u;
//                // TODO 把顶点剩余度数清空，把边不可用
//                if (window.size() < WINDOW_SIZE) {
//                    // TODO 拷贝
//                    add_edge_to_window(edges[neighbors[0].v]);
//                    LOG(INFO) << window.size();
//                    swap(neighbors[0], neighbors.back());
//                    neighbors.pop_back();
//                    adjlist_t &u_neighbors = adj_in[u];
//                    swap(u_neighbors[0], u_neighbors.back());
//                    u_neighbors.pop_back();
//                    edges[neighbors[0].v].remove();
//                    LOG(INFO) << "Push edge into window.";
//                    re_get = true;
//                }else {
//                    LOG(INFO) << "Window is full";
//                }
//            }
//            else vid = u;
//        } else if (adj_in[vid].size() == 1){
//            // 获取对端顶点
//            adjlist_t &neighbors = adj_in[vid];
//            vid_t &u = edges[neighbors[0].v].first;
//            if (partial_degree[u] == 1) {
//                LOG(INFO) << "U size is 1: " << u;
//                if (window.size() < WINDOW_SIZE) {
//                    add_edge_to_window(edges[neighbors[0].v]);
//                    // LOG(INFO) << window.size();
//                    swap(neighbors[0], neighbors.back());
//                    neighbors.pop_back();
//                    adjlist_t &u_neighbors = adj_out[u];
//                    swap(u_neighbors[0], u_neighbors.back());
//                    u_neighbors.pop_back();
//                    edges[neighbors[0].v].remove(); // TODO remove会改变边，不能直接remove
//                    LOG(INFO) << "Push edge into window, window size: " << window.size();
//                    re_get = true;
//                } else {
//                    LOG(INFO) << "Window is full";
//                }
//            }
//            else vid = u;
//        } else {
//            LOG(INFO) << "Error, vid: " << vid << ", partial degree: " << adj_in[vid].size() + adj_out[vid].size();
//        }
//
//    }
//    if (re_get) return get_free_vertex(vid);
    return true;
}

double OffstreamNAPartitioner::calculate_lb_score(size_t partition_id) {

    double  lb_score = (double)max_load - occupied[partition_id];
    if (min_load != UINT64_MAX) {
        lb_score /= (epsilon + (double)max_load - (double)min_load);
    }

    return lb_score;
}

double OffstreamNAPartitioner::calculate_rf_score(vid_t u, vid_t v, size_t partition_id) {
    double gu = 0, gv = 0;
    size_t u_degree = partial_degree[u];
    size_t v_degree = partial_degree[v];
    size_t sum = v_degree + u_degree;
    // 归一化
    if (is_mirrors[u].get(partition_id)) {
        gu = (double)u_degree;
        gu /= (double)sum;
        gu = 1 + (1 - gu);
    }
    if (is_mirrors[v].get(partition_id)) {
        gv = (double)v_degree;
        gv /= (double)sum;
        gv = 1 + (1 - gv);
    }
    double rf_score = gu + gv;
    return rf_score;
}

double OffstreamNAPartitioner::calculate_cs_score(vid_t u, vid_t v, size_t partition_id) {
    // 拿到 v -> p
    vector<size_t>& v_neighbor_list = vertex_partitions[v];
    vector<size_t>& u_neighbor_list = vertex_partitions[u];
    size_t v_neighbors = 0;
    size_t u_neighbors = 0;
    for(size_t i = 0; i < num_partitions; i++) {
        v_neighbors += v_neighbor_list[i];
        u_neighbors += u_neighbor_list[i];
    }
    // TODO 不去重
    size_t total_neighbors = v_neighbors + u_neighbors;
    // TODO 为什么会大于0
    double cs_score = (double)(v_neighbor_list[partition_id] + u_neighbor_list[partition_id]) / (double)total_neighbors;
    return cs_score;
}

void OffstreamNAPartitioner::occupy_vertex(vid_t vid, vid_t degree) {
    CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
    is_cores[current_partition].set_bit_unsync(vid);
    if (degree == 0) return;
    add_boundary(vid);

    for (auto &i: adj_out[vid]) {
        if (edges[i.v].valid())
            add_boundary(edges[i.v].second);
    }
    adj_out[vid].clear();

    for (auto &i: adj_in[vid]) {
        if (edges[i.v].valid())
            add_boundary(edges[i.v].first);
    }
    adj_in[vid].clear();
}

void OffstreamNAPartitioner::assign_edge(size_t partition_id, vid_t from, vid_t to) {
    // assign_edge更新分区相关的信息
    is_mirrors[from].set_bit_unsync(partition_id);
    is_mirrors[to].set_bit_unsync(partition_id);

    vp_set[from].set_bit_unsync(partition_id);
    vp_set[to].set_bit_unsync(partition_id);
    vertex_partitions[from][partition_id] += 1;
    vertex_partitions[to][partition_id] += 1;

    assigned_edges++;
    occupied[partition_id]++;
}

void OffstreamNAPartitioner::remove_edge_from_window() {
    Timer buffer_calculate_cost;
    buffer_calculate_cost.start();

    double global_max_score = -0.2;
    size_t global_max_partition = 0;
    edge_t global_max_edge;
    size_t global_max_index = 0;
    for (int idx = 0; idx < window.size(); idx++) {
        auto& edge = window[idx];

        vid_t u = edge.first;
        vid_t v = edge.second;

        double local_max_score = -0.1;
        size_t local_max_partition = 0;
        size_t local_max_index = 0;
        for(int p = 0; p < num_partitions; p++) {
            if (occupied[p] >= max_partition_load) {
                continue;
            }
            if (local_max_partition == num_partitions) {
                local_max_partition = p;
            }

            // 计算三部分的分数
            double lb_score = calculate_lb_score(p);
            double rf_score = calculate_rf_score(u,v,p);
            double cs_score = calculate_cs_score(u, v, p);
            // LOG(INFO) << "lb_score: " << lb_score << " rf_score: " << rf_score << " | cs_score: " << cs_score;
            double local_score = lambda * lb_score + rf_score + cs_score;
            if (local_score > local_max_score) {
                local_max_partition = p;
                local_max_score = local_score;
                local_max_index = idx;
            }
        }
        if (local_max_score > global_max_score) {
            global_max_score = local_max_score;
            global_max_partition = local_max_partition;
            global_max_edge = edge;
            global_max_index = local_max_index;
        }
    }
    // 将边分配到分区
    assign_edge(global_max_partition, global_max_edge.first, global_max_edge.second);
    // 从缓存窗口移除边
    update_min_max_load(global_max_partition);
    window[global_max_index] = window.back();
    window.pop_back();

    buffer_calculate_cost.stop();
    // LOG(INFO) << "Buffer calculate cost time: " << buffer_calculate_cost.get_time();
}

void OffstreamNAPartitioner::add_edge_to_window(edge_t& edge) {
    window.push_back(edge);
    // 先更新度数
    vid_t u = edge.first;
    vid_t v = edge.second;
    partial_degree[u]++;
    partial_degree[v]++;
}

void OffstreamNAPartitioner::add_boundary(vid_t vid) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[current_partition];
    auto &is_boundary = is_boundaries[current_partition];

    if (is_boundary.get(vid)) {
        return;
    }
    is_boundary.set_bit_unsync(vid);

    if (!is_core.get(vid)) {
        min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
    }

    rep (direction, 2) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        // 遍历顶点的邻边
        for (size_t i = 0; i < neighbors.size();) {
            // 判断邻居edges[neighbors[i].v]是否已经分配
            if (edges[neighbors[i].v].valid()) {
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                if (is_core.get(u)) { // 如果顶点在核心集中，不考虑负载
                    assign_edge(current_partition, direction ? vid : u,
                                direction ? u : vid);

                    min_heap.decrease_key(vid); // 默认移除一条边
                    edges[neighbors[i].v].remove();
                    swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else if (is_boundary.get(u) &&
                           occupied[current_partition] < capacity) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                    // 将边加入边集
                    assign_edge(current_partition, direction ? vid : u, direction ? u : vid);
                    min_heap.decrease_key(vid);
                    min_heap.decrease_key(u);
                    edges[neighbors[i].v].remove();
                    swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else
                    i++;
            } else {
                std::swap(neighbors[i], neighbors.back());
                neighbors.pop_back();
            }
        }
    }
}

int OffstreamNAPartitioner::find_max_score_partition(edge_t &e) {
    auto degree_u = partial_degree[e.first];
    auto degree_v = partial_degree[e.second];

    uint32_t sum;
    double max_score = 0;
    // TODO 这里默认分区为0
    uint32_t max_p = 0;
    double bal, gv, gu;

    for (int j = 0; j < num_partitions; j++) {
        // 如果分区已经超过了最大分区负载，直接跳过
        if (occupied[j] >= max_partition_load) {
            continue;
        }
        if (max_p == num_partitions) {
            max_p = j;
        }
        // 以下对应着hdrf算法的核心实现
        gu = 0, gv = 0;
        sum = degree_u + degree_v;
        // 归一化
        if (is_mirrors[e.first].get(j)) {
            gu = degree_u;
            gu /= sum;
            gu = 1 + (1 - gu);
        }
        if (is_mirrors[e.second].get(j)) {
            gv = degree_v;
            gv /= sum;
            gv = 1 + (1 - gv);
        }
        double rep = gu + gv; // rep值
        bal = max_load - occupied[j];
        if (min_load != UINT64_MAX) {
            bal /= (epsilon + max_load - min_load);
        }
        // 计算结果应该有两部分组成，rep和bal
        // LOG(INFO) << "rep: " << rep << " bal: " << lambda * bal;
        double score_p = rep + lambda * bal;
        // LOG(INFO) << "score_p: " << score_p;
        CHECK_GE(score_p, 0) << "score_p: " << score_p;
        if (score_p > max_score) {
            max_score = score_p;
            max_p = j;
        }
    }
    return max_p;
}

