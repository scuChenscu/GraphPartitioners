#include "hdrf.hpp"

using namespace std;

// HDRF是以边作为输入流的边分区算法，是为幂律分布图设计的图分区方法。其基本思路是优先对度数高的顶点进行切分，这样可以最小化镜像顶点的数量。
// HDRF是把边划分到不同的分区，即存在复制vertex
HdrfPartitioner::HdrfPartitioner(const BaseGraph& baseGraph, const string& input, const string& algorithm,
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

    max_partition_load = (uint64_t)balance_ratio * num_edges / num_partitions;
    // vertex度数
    degrees.resize(num_vertices, 0);
    // batch数，根据内存来划分
    num_batches = (filesize / ((size_t) memory_size * 1024 * 1024)) + 1;
    num_edges_per_batch = (num_edges / num_batches) + 1;
    // 记录每个分区的边负载
    edge_load.resize(num_partitions);
    // 应该是记录每个vertex的分区，用位图来记录，类似 0 0 0 0 1 0
    is_mirrors.assign(num_vertices, dense_bitset(num_partitions));
    true_vids.resize(num_vertices);
    // 分区度？不是很懂这个的作用
    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);
}

void HdrfPartitioner::batch_hdrf(vector<edge_t> &edges) {
    // 计算vertex的度
    for (auto &e: edges) {
        ++degrees[e.first];
        ++degrees[e.second];

        int max_p = find_max_score_partition_hdrf(e);
        update_is_mirrors(e, max_p);
        update_min_max_load(max_p);
        // save_edge(e.first,e.second,max_p);
        ++part_degrees[e.first][max_p];
        ++part_degrees[e.second][max_p];
//        save_vertex(e.first,max_p);
//        save_vertex(e.second,max_p);
    }
}

// 选择得分最高的分区作为边的目标分区，得分主要由bal和rep两部分组成
int HdrfPartitioner::find_max_score_partition_hdrf(edge_t &e) {
    auto degree_u = degrees[e.first];
    auto degree_v = degrees[e.second];

    uint32_t sum;
    double max_score;
    uint32_t max_p = 0;
    double bal, gv, gu;

    for (int j = 0; j < num_partitions; j++) {
        // 如果分区已经超过了最大分区负载，直接跳过
        if (edge_load[j] >= max_partition_load) {
            continue;
        }
        // 以下对应着hdrf算法的核心实现
        gu = 0, gv = 0;
        sum = degree_u + degree_v;
        // 这一步是在干嘛？判断他们是不是在这个分区有？
        // 对边的两个顶点，如果在该分区有副本，计算g：先归一化，再计算，如果没有副本，g=0；
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
        if (min_load != UINT64_MAX && max_load != 0) {
            bal = (max_load - edge_load[j]) / (epsilon + max_load - min_load);
        } else {
            bal = max_load - edge_load[j];
        }
        // bal = (max_load - edge_load[p]) / (epsilon + max_load - min_load);
        // 计算结果应该有两部分组成，rep和bal
        double score_p = rep + lambda * bal;
        if (score_p < 0) {
            LOG(ERROR) << "ERROR: score_p < 0";
            LOG(ERROR) << "gu: " << gu;
            LOG(ERROR) << "gv: " << gv;
            LOG(ERROR) << "bal: " << bal;
            exit(-1);
        }
        if (score_p > max_score) {
            max_score = score_p;
            max_p = j;
        }
    }
    return max_p;
}

void HdrfPartitioner::update_is_mirrors(edge_t &e, int max_p) {
    is_mirrors[e.first].set_bit_unsync(max_p);
    is_mirrors[e.second].set_bit_unsync(max_p);
    true_vids.set_bit_unsync(e.first);
    true_vids.set_bit_unsync(e.second);
}

void HdrfPartitioner::update_min_max_load(int max_p) {
    auto &load = ++edge_load[max_p];
    if (load > max_load) max_load = load;
    min_load = *min_element(edge_load.begin(), edge_load.end());
}

void HdrfPartitioner::batch_node_assign_neighbors(vector<edge_t> &edges) {
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
    }
}

void HdrfPartitioner::read_and_do(const string& opt_name) {
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), ios::beg);
    vector<edge_t> edges;
    auto num_edges_left = num_edges;
    for (uint32_t i = 0; i < num_batches; i++) {
        auto edges_per_batch = num_edges_per_batch < num_edges_left ? num_edges_per_batch : num_edges_left;
        edges.resize(edges_per_batch);
        fin.read((char *) &edges[0], sizeof(edge_t) * edges_per_batch);
        if (opt_name == "hdrf") {
            batch_hdrf(edges);
        } else if (opt_name == "node_assignment") {
            batch_node_assign_neighbors(edges);
        } else {
            LOG(ERROR) << "no valid opt function";
        }
        num_edges_left -= edges_per_batch;
    }
}

void HdrfPartitioner::split() {

    stringstream ss;
    ss << "HDRF" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());


    Timer total_time;
    total_time.start();

    read_and_do("hdrf");

    //根据结点平衡性、随机分配的重叠度以及结点的度大小来判断
    size_t total_mirrors = 0;
    vector<vid_t> buckets(num_partitions);
    // 每个分区的容量
    double capacity = (double) true_vids.popcount() * 1.05 / num_partitions + 1;
    // 遍历vertex
    rep(i, num_vertices) {
        total_mirrors += is_mirrors[i].popcount(); // 这个是不是说明被复制了？
        double max_score = 0.0;
        vid_t which_p;
        bool unique = false;
        if (is_mirrors[i].popcount() == 1) { // 说明只有一个分区，没有复制
            unique = true;
        }
        // 遍历分区
        repv(j, num_partitions) {
            // 如果vertex i的j个分区为1
            if (is_mirrors[i].get(j)) {
//                double score=((i%p==j)?1:0)+(part_degrees[i][j]/(degrees[i]+1))+(buckets[j]< capacity?1:0);
                // 每个分区节点的度
                // TODO 该算法好像没有更新过part_degrees
                double score = (part_degrees[i][j] / (degrees[i] + 1)) + (buckets[j] < capacity ? 1 : 0);
                if (unique) {
                    which_p = j;
                } else if (max_score < score) {
                    max_score = score;
                    which_p = j;
                }
            }
        }
        ++buckets[which_p];
        save_vertex(i, which_p); // 这里为什么又是点分区
        balance_vertex_distribute[i] = which_p;
    }
    vertex_ofstream.close();
    repv(j, num_partitions) {
        LOG(INFO) << "each partition node count: " << buckets[j];
    }

    read_and_do("node_assignment");
    edge_ofstream.close();

    total_time.stop();
    // rep(i, p) LOG(INFO) << "edges in partition " << i << ": " << edge_load[i];
    // LOG(INFO) << "replication factor: " << (double)total_mirrors / true_vids.popcount();
    LOG(INFO) << "total partition time: " << total_time.get_time();

    calculate_replication_factor();
    stringstream result;
    result << "Cost Time: " << total_time.get_time()
           << " | Replication Factor: " << replication_factor
           << " | Epsilon: " << epsilon
           << " | Lambda: " << lambda
           << " | Max Load:" << max_load
           << " | Min Load: " << min_load
           << endl;
    appendToFile(result.str());
}