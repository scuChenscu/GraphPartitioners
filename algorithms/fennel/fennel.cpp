#include "fennel.hpp"


// FENNEL算法的主要思想：为了减少边切割，一个顶点应该被分配到有较多邻居的分区，同时还要加入一个惩罚因子，以防止一个分区变得过大
FennelPartitioner::FennelPartitioner(BaseGraph& baseGraph,const string& input, const string& algorithm, const size_t num_partitions, int memory_size, bool shuffle) :
        VertexPartitioner(baseGraph, algorithm, num_partitions) {
    config_output_files();
    LOG(INFO) << "begin init class";
    // TODO 为什么在这里就开始计算耗时
    //edge file
    if (shuffle) {
        fin.open(shuffled_binary_edgelist_name(input), std::ios::binary | std::ios::ate);
    } else {
        fin.open(binary_edgelist_name(input), std::ios::binary | std::ios::ate);
    }
    // 获取文件指针fin当前所处位置的文件偏移量
    filesize = fin.tellg();
    fin.seekg(0, std::ios::beg);

    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));
    // 文件大小 除以 内存大小 得到 批次数
    num_batches = (filesize / ((std::size_t) memory_size * 1024 * 1024)) + 1;
    // 每个batch的边数
    // TODO batch的作用是什么？内存不够，每次只能读取这么多？
    num_edges_per_batch = (num_edges / num_batches) + 1;

//    subg_vids.assign(p,dense_bitset(num_vertices));
//    true_vids.resize(num_vertices);
//    balance_vertex_distribute.resize(num_vertices);
//    node2neis.assign(num_vertices,dense_bitset(num_vertices));

    //unorderedset
    subg_vids.assign(p, unordered_set < vid_t > {});
//    true_vids.resize(num_vertices);

    balance_vertex_distribute.resize(num_vertices);

    node2neis.assign(num_vertices, unordered_set < vid_t > {});
    LOG(INFO) << "finish init";
}

// read_and_do主要有两种操作，分别是1. process neighbors 和2. node_assignment
void FennelPartitioner::read_and_do(string opt_name) {
    // seekg 将文件指针移动到指定位置的函数, 从 begin 开始移动到 offset位置，就是跳过num
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
    std::vector<edge_t> edges;
    auto num_edges_left = num_edges;
    // 遍历每一个batch：process neighbors就是（1）对每个vid建立邻接表，（2）unique vids
    // node_assignment
    for (uint32_t i = 0; i < num_batches; i++) {
        // 取较小值
        auto edges_per_batch = num_edges_per_batch < num_edges_left ? num_edges_per_batch : num_edges_left;
        edges.resize(edges_per_batch);
        // 把edge读入到edges，就是内存操作
        fin.read((char *) &edges[0], sizeof(edge_t) * edges_per_batch);
        if (opt_name == "node_assignment") {
            batch_node_assignment(edges);
        } else if (opt_name == "process neighbors") {
            for (auto &e: edges) {
                addNeighbors(e);
                true_vids.insert(e.first);
                true_vids.insert(e.second);
            }
            LOG(INFO) << "finish neis";
        } else {
            LOG(ERROR) << "no valid opt function";
        }
        num_edges_left -= edges_per_batch;
    }
}

// 根据edge的first和second，写入vertex的邻居信息，像是在构造邻接表
void FennelPartitioner::addNeighbors(edge_t &edge) {
    node2neis[edge.first].insert(edge.second);
    node2neis[edge.second].insert(edge.first);
}

size_t FennelPartitioner::intersection(vid_t vid, size_t partition) {
    if (!adjacency_list.count(vid)) return 0;
    int count = 0;
    for (auto neighbor : adjacency_list.find(vid)->second) {
        if (partition_vertices[partition].get(neighbor)) {
            count++;
        }
    }
    // LOG(INFO) << count;
    return count;
}

void FennelPartitioner::do_fennel() {
    // shuffle_vertices of streaming vertices
    vector<int> shuffle_vertices(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        shuffle_vertices[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(shuffle_vertices.begin(), shuffle_vertices.end(), g);

    // Initial partitions
    size_t i = 0;
    for (; i < num_partitions; i++) {
        assigned_vertex(shuffle_vertices[i], i);
    }

    const double gamma = 3 / 2.0;
    const double alpha =
            num_edges * pow(num_partitions, (gamma - 1)) / pow(num_vertices, gamma);
    const double load_limit = 1.1 * num_vertices / num_partitions;
    // TODO 感觉这里是不是会对前面的p个顶点进行冗余分区
    for (;i < num_vertices; i++) {
        vid_t v = shuffle_vertices[i];
        if (i % 10000 == 0)
            cout << i << "/" << num_vertices << endl;
            vector<double> from_scores(num_partitions, 0);
            // 计算每个分区的得分
            for (int id = 0; id < num_partitions; id++) {
                double partitionSize = num_vertices_each_partition[id];
                if (partitionSize <= load_limit) {
                    double firstVertexIntraCost;
                    double firstVertexInterCost = intersection(v, id);
                    firstVertexIntraCost = alpha * gamma * pow(partitionSize, gamma - 1);
                    from_scores[id] = firstVertexInterCost - firstVertexIntraCost;
                }
            }
            //最大值所在序列的位置
            int firstIndex = distance(from_scores.begin(),
                                      max_element(from_scores.begin(), from_scores.end()));
        assigned_vertex(v, firstIndex);
    }
}

// TODO 这个划分方法有点奇怪
// Fennel算法是把点划分到不同分分区
void FennelPartitioner::batch_node_assignment(vector<edge_t> &edges) {
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
        LOG(INFO) << sp << tp << endl;
        if (sp != tp) {
            edge_cut++;
        }
    }
    edge_cut_rate= (double) edge_cut / edges.size();
}

void FennelPartitioner::split() {
    stringstream ss;
    ss << "Fennel" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());
    total_time.start();
    do_fennel();
    total_time.stop();
}