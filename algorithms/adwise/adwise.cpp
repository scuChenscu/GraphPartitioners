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
        if (get_window_size() < w) {
            // add_edge_into_window();
            C.insert(index++);
        }
        // 选择最好的e-p
        int e_id, p_id;
        double score = get_best_assignment(e_id, p_id);
        w_scores += score;
//        LOG(INFO) << "P id: " << p_id;
//        remove_edge_from_window(e_id, score);
        // 分配
        // edge_t &edge = edges[e_id];
        assign_edge(e_id, p_id);
        // 如果是新增副本，需要重新计算Q
        // TODO 直接不算了
//        if (!is_mirrors[edge.first].get(p_id)) {
//            // 需要计算的度得分和簇得分
//            // 遍历Q
//        }
//        if(!is_mirrors[edge.second].get(p_id)) {
//
//        }

    }
    LOG(INFO) << "Assign window";
    // 直接分配完
     if(get_window_size() > 0) {
         C.insert(Q.begin(), Q.end());
         while (!C.empty()) {
             // 分配边
             int e_id, p_id;
             get_best_assignment(e_id, p_id);
             assign_edge(e_id, p_id);
         }
     }

    total_time.stop();
}

void AdwisePartitioner::assign_edge(int e_id, int p_id) {
    edge_t &edge = edges[e_id];
    is_mirrors[edge.first].set_bit_unsync(p_id);
    is_mirrors[edge.second].set_bit_unsync(p_id);
    occupied[p_id]++;
    assigned_edges++;
    update_min_max_load(p_id);
}

double AdwisePartitioner::get_best_assignment(int &e_id, int &p_id){
    // TODO 这一步计算，类似HDRF，但是HDRF只计算一个，这个计算整个C
    double max_score = 0.0;
    for (const auto &ce_id : C) {
        int p;
        double score = calculate_score(ce_id, p);
        if (score > max_score) {
            max_score = score;
            e_id = ce_id;
            p_id = p;
        }
    }
    C.erase(e_id);

//    if (c % w == 0) { // 可以用HDRF和NE算法来估计这个时间
//        // C1：计算分配w条边的平均得分；
//        // C2：剩余的边数和距离超过初始估计上限时间；
//        if (w_scores / c > last_scores) {
//            w = 2 * w;
//            // c = 0;
//            while(get_window_size() < w && index < num_edges) {
//                // add_edge_into_window();
//                C.insert(index);
//                index++;
//            }
//        }
//        // 重置分数
//        last_scores = scores;
//        w_scores = 0.0;
//        // TODO 时间上限收缩窗口
////        else if (false) {
////            // 向上取整
////            int res = w % 2;
////            w = w / 2 + res;
////            c = 0;
////        }
//    }
    c++;
}

void AdwisePartitioner::remove_edge_from_window(int e_id, double sc) {
    w_scores += sc;
    scores = scores - sc;
    update_threshold();
    C.erase(e_id);
    // 更新threshold
    if (C.empty()) { // 计算Q
        // 从Q -> C，将大于threshold的边移入C
        for(const auto &qe_id : Q) {
            edge_t &e = edges[qe_id];
            int p;
            double score = calculate_score(qe_id, p);
            if (score > threshold + 1) {
                C.insert(qe_id);
                Q.erase(qe_id);
            }
        }
    }
}

void AdwisePartitioner::add_edge_into_window() {
    int pid;
    double sc = calculate_score(index, pid);
    // 计算加入C还是Q
    // 计算该边加入所有分区的最高得分，如果高于threshold，在C，否则在Q
    // TODO threshold始终为w中的平均值
    C.insert(index);
//    if (sc > threshold + 1) {
//        C.insert(index);
//    } else {
//        Q.insert(index);
//    }
//    // 更新threshold
//    scores += sc;
//    update_threshold();
    index++;
}

int AdwisePartitioner::get_window_size() {
    return C.size() + Q.size();
}

void AdwisePartitioner::update_threshold() {
    threshold = scores / (double)get_window_size();
}

double AdwisePartitioner::calculate_score(int e_id, int &pid) {
    edge_t &e = edges[e_id];
    vid_t u = e.first;
    vid_t v = e.second;
    // 得分由三部分组成，平衡得分，度得分，簇得分
    double max_score = -1;
    for (int i = 0; i < num_partitions; i++) {
        // 平衡得分
        double bal;
        bal = max_load - edge_load[i];
        if (min_load != UINT64_MAX) {
            bal /= (epsilon + max_load - min_load);
        }
        // 度得分
        double rep;
        double du = 0;
        double dv = 0;
        if (is_mirrors[u].get(i)) {
            du = 2 - (double)degrees[u] / (double)(2 * max_degree);
        }
        if (is_mirrors[v].get(i)) {
            dv = 2 - (double)degrees[v] / (double)(2 * max_degree);
        }
        rep = du + dv;
        // 簇得分
        double cs;
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

        cs = (double)nc / (double)(n + 0.1); // 避免除以0
        if (rep > 0) {
            LOG(INFO) << "bal: " << bal << " rep: " << rep << " cs: " << cs;

        }

        double sc = 2.8 * bal + rep + cs;
        if (sc > max_score) {
            max_score = sc;
            pid = i;
        }
    }
    // LOG(INFO) << "P id: " << pid;
    return max_score;
}