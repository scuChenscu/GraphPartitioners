#include "adwise.hpp"

using namespace std;

// HDRF是以边作为输入流的边分区算法，是为幂律分布图设计的图分区方法。其基本思路是优先对度数高的顶点进行切分，这样可以最小化镜像顶点的数量。
// HDRF是把边划分到不同的分区，即存在复制vertex
AdwisePartitioner::AdwisePartitioner(BaseGraph &baseGraph, const string &input, const string &algorithm,
                                 size_t num_partitions,
                                 int memory_size,
                                 double balance_ratio,
                                 double balance_lambda,
                                 bool shuffle) : EdgePartitioner(baseGraph, algorithm, num_partitions) {

    config_output_files();

    lambda = balance_lambda;

    if (shuffle) {
        fin.open(shuffled_binary_edgelist_name(input), ios::binary | ios::ate);
    } else {
        fin.open(binary_edgelist_name(input), ios::binary | ios::ate);
    }

    filesize = fin.tellg();
    fin.seekg(0, ios::beg);

    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));

    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t), filesize);

    max_partition_load = (uint64_t) balance_ratio * num_edges / num_partitions;
    // vertex度数
    degrees.resize(num_vertices, 0);
    // batch数，根据内存来划分
    num_batches = (filesize / ((size_t) memory_size * 1024 * 1024)) + 1;
    num_edges_per_batch = (num_edges / num_batches) + 1;
    // 记录每个分区的边负载
    edge_load.resize(num_partitions);
    // 应该是记录每个vertex的分区，用位图来记录，类似 0 0 0 0 1 0
    // is_mirrors.assign(num_vertices, dense_bitset(num_partitions));
    // true_vids.resize(num_vertices);
    // 分区度？不是很懂这个的作用
    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);

    partial_degrees.resize(num_vertices);

    max_degree = baseGraph.max_degree;
}

void AdwisePartitioner::batch_adwise(vector<edge_t> &edges) {
    // 计算vertex的度
    for (auto &e: edges) {
        ++degrees[e.first];
        ++degrees[e.second];

        int max_p = find_max_score_partition(e);
        // update_is_mirrors(e, max_p);
        update_min_max_load(max_p);
        // save_edge(e.first,e.second,max_p);
        ++part_degrees[e.first][max_p];
        ++part_degrees[e.second][max_p];
//        save_vertex(e.first,max_p);
//        save_vertex(e.second,max_p);
    }
}

// 选择得分最高的分区作为边的目标分区，得分主要由bal和rep两部分组成
int AdwisePartitioner::find_max_score_partition(edge_t &e) {
    auto degree_u = ++partial_degrees[e.first];
    auto degree_v = ++partial_degrees[e.second];

    uint32_t sum;
    double max_score = 0;
    // TODO 这里默认分区为0
    uint32_t max_p = 0;
    double bal, gv, gu;

    for (int j = 0; j < num_partitions; j++) {
        // 如果分区已经超过了最大分区负载，直接跳过
        if (edge_load[j] >= max_partition_load) {
            continue;
        }
        if (max_p == num_partitions) {
            max_p = j;
        }
        // 以下对应着adwise算法的核心实现
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
        bal = max_load - edge_load[j];
        if (min_load != UINT64_MAX) {
            bal /= (epsilon + max_load - min_load);
        }
        // 计算结果应该有两部分组成，rep和bal
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

void AdwisePartitioner::update_is_mirrors(edge_t &e, int max_p) {
    is_mirrors[e.first].set_bit_unsync(max_p);
    is_mirrors[e.second].set_bit_unsync(max_p);
    true_vids.set_bit_unsync(e.first);
    true_vids.set_bit_unsync(e.second);
}

void AdwisePartitioner::update_min_max_load(int max_p) {
    auto &load = ++edge_load[max_p];
    if (load > max_load) max_load = load;
    min_load = *min_element(edge_load.begin(), edge_load.end());
}

void AdwisePartitioner::batch_node_assign_neighbors(vector<edge_t> &edges) {
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
    }
}

void AdwisePartitioner::read_and_do(const string &opt_name) {
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), ios::beg);
    vector<edge_t> edges;
    auto num_edges_left = num_edges;
    for (uint32_t i = 0; i < num_batches; i++) {
        auto edges_per_batch = num_edges_per_batch < num_edges_left ? num_edges_per_batch : num_edges_left;
        edges.resize(edges_per_batch);
        fin.read((char *) &edges[0], sizeof(edge_t) * edges_per_batch);
        if (opt_name == "adwise") {
            batch_adwise(edges);
        } else if (opt_name == "node_assignment") {
            batch_node_assign_neighbors(edges);
        } else {
            LOG(ERROR) << "no valid opt function";
        }
        num_edges_left -= edges_per_batch;
    }
}

void AdwisePartitioner::split() {

    stringstream ss;
    ss << "Adwise: " << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    total_time.start();

    // 不断地移入移除边，直到满足w次操作，尝试调整窗口大小，并且满载w
    // 移入移除边需要去更新threshold
    while(index < num_edges) {
        if (W.size() < w) {
            double sc = calculate_score(index);
            // 计算加入C还是Q
            // 计算该边加入所有分区的最高得分，如果高于threshold，在C，否则在Q
            // TODO threshold始终为w中的平均值
            if (sc > threshold + 1) {
                C.insert(index++);
            } else {
                Q.insert(index++);
            }
            // 更新threshold
            threshold = (threshold * W.size() + sc) / (W.size() + 1) ;
            W.insert(index);
        }
        // 选择最好的e-p
        int e_id, p_id;
        get_best_assignment(e_id, p_id);
        // 分配
        edge_t &edge = edges[e_id];
        is_mirrors[edge.first].set_bit_unsync(p_id);
        is_mirrors[edge.second].set_bit_unsync(p_id);
        occupied[p_id]++;
        assigned_edges++;
        update_min_max_load(p_id);
        // 分配边
        // TODO 更新threshold
    }
    total_time.stop();
}

void AdwisePartitioner::get_best_assignment(int &e_id, int &p_id){
    // 不用计算，从C中选出最好的结果

    W.erase(e_id);
    C.erase(e_id);
    if (C.empty()) { // 计算Q

    }
    if (assigned_edges % w == 0) { // 可以用HDRF和NE算法来估计这个时间
        // C1：计算分配w条边的平均得分；
        // C2：剩余的边数和距离超过初始估计上限时间；
        if (true && true) {
            w = 2 * w;
            while(W.size() < w && index < num_edges) {
                W.insert(index++);
            }
        } else if (true) {
            // 向上取整
            int res = w % 2;
            w = w / 2 + res;
        }
    }
}

void AdwisePartitioner::calculate_set() {

}

double AdwisePartitioner::calculate_score(int e_id) {
    edge_t &e = edges[e_id];
    vid_t u = e.first;
    vid_t v = e.second;
    // 得分由三部分组成，平衡得分，度得分，簇得分
    double max_score = 0;
    for (int i = 0; i < num_partitions; i++) {
        // 平衡得分
        double bal, deg;
        bal = max_load - edge_load[i];
        if (min_load != UINT64_MAX) {
            bal /= (epsilon + max_load - min_load);
        }
        // 度得分
        double rep = 0;
        double du = 0;
        double dv = 0;
        if (is_mirrors[u].get(i)) {
            du = 2 - degrees[u] / (2 * max_degree);
        }
        if (is_mirrors[v].get(i)) {
            dv = 2 - degrees[v] / (2 * max_degree);
        }
        rep = du + dv;
        // 簇得分
        double cs = 0;
        int nc = 0;
        int n = 0;
        // 计算u,v 在当前分区的邻居数
        if (is_mirrors[u].get(i)) {
            rep (direction, 2) {
                adjlist_t &neighbors = direction ? adj_out[u] : adj_in[u];
                n += neighbors.size();
                // 遍历顶点的邻边
                for (size_t idx = 0; idx < neighbors.size(); idx++) {
                    vid_t &neighbor = direction ? edges[neighbors[idx].v].second : edges[neighbors[idx].v].first;
                    if (is_mirrors[neighbor].get(i)) {
                        nc++;
                    }
                }

            }
        }
        if (is_mirrors[v].get(i)) {
            rep (direction, 2) {
                adjlist_t &neighbors = direction ? adj_out[v] : adj_in[v];
                n += neighbors.size();
                // 遍历顶点的邻边
                for (size_t idx = 0; idx < neighbors.size(); idx++) {
                    vid_t &neighbor = direction ? edges[neighbors[idx].v].second : edges[neighbors[idx].v].first;
                    if (is_mirrors[neighbor].get(i)) {
                        nc++;
                    }
                }
            }
        }
        cs = (double)nc / (double)n;

        double sc = bal + rep + cs;
        if (sc > max_score) {
            max_score = sc;
        }
    }
    return max_score;
}