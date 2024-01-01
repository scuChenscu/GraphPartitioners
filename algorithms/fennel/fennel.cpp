//
// Created by muzongshen on 2021/9/30.
//

#include "fennel.hpp"
#include <algorithm>
#include <cmath>
#include <random>

// FENNEL算法的主要思想：为了减少边切割，一个顶点应该被分配到有较多邻居的分区，同时还要加入一个惩罚因子，以防止一个分区变得过大
FennelPartitioner::FennelPartitioner(string input, string algorithm, int num_partition, int memsize, bool shuffle) {
    p = num_partition;
    config_output_files(input, algorithm, num_partition);
    LOG(INFO) << "begin init class";
    // TODO 为什么在这里就开始计算耗时
    total_time.start();
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
    num_batches = (filesize / ((std::size_t) memsize * 1024 * 1024)) + 1;
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

int FennelPartitioner::intersection(unordered_set<vid_t> &nums1, unordered_set<vid_t> &nums2) {
    // 建立unordered_set存储nums1数组(清除了重复的元素)
    unordered_set < vid_t > ans;
    for (auto num: nums2) {
        if (nums1.count(num) == 1)
            ans.insert(num);
    }
    return ans.size();
}

void FennelPartitioner::do_fennel() {
    // Ordering of streaming vertices
    // TODO 为什么不用true_vcount
    vector<int> ordering(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        ordering[i] = i;
    }
    // TODO 设置随机种子，乱序所有顶点
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(ordering.begin(), ordering.end(), g);

    // Initial partitions
    // TODO 这段代码的含义？为了保证每个分区初始顶点数大于0？
    for (int i = 0; i < p; ++i) {
        // subg_vids的size为p, 插入前p个顶点
        subg_vids[i].insert(ordering[i]);
        save_vertex(ordering[i], i);
    }

    int true_vcount = true_vids.size();
    const double gamma = 3 / 2.0;
    const double alpha =
            num_edges * pow(p, (gamma - 1)) / pow(true_vcount, gamma);
    const double load_limit = 1.1 * true_vcount / p;
    // TODO 感觉这里是不是会对前面的p个顶点进行冗余分区
    for (int i = 0; i < (vid_t) num_vertices; i++) {
        vid_t v = ordering[i];
        if (i % 10000 == 0)
            cout << i << "/" << num_vertices << endl;
        if (true_vids.find(v) != true_vids.end()) {
            vector<double> from_scores(p, 0);
            // 计算每个分区的得分
            for (int id = 0; id < p; id++) {
                double partitionSize = subg_vids[id].size();
                if (partitionSize <= load_limit) {
                    double firstVertexIntraCost;
                    double weightedGreedy =
                            (1 - (partitionSize / ((double) true_vcount / (double) p)));
                    double firstVertextInterCost = intersection(node2neis[v], subg_vids[id]);
                    firstVertexIntraCost = alpha * gamma * pow(partitionSize, gamma - 1);
                    from_scores[id] = firstVertextInterCost - firstVertexIntraCost;
                }
            }
            //最大值所在序列的位置
            int firstIndex = distance(from_scores.begin(),
                                      max_element(from_scores.begin(), from_scores.end()));
            balance_vertex_distribute[v] = firstIndex;
            subg_vids[firstIndex].insert(v);
            save_vertex(v, firstIndex);
        }
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
    LOG(INFO) << "begin process neighbors";

    stringstream ss;
    ss << "Fennel" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());


    read_and_do("process neighbors");
    LOG(INFO) << "begin write nodes";
    do_fennel();
    LOG(INFO) << "begin write edges";
    read_and_do("node_assignment");
    total_time.stop();
    edge_ofstream.close();
    LOG(INFO) << "total vertex count: " << true_vids.size();
    LOG(INFO) << "total partition time: " << total_time.get_time();

    stringstream result;
    result << "Cost Time: " << total_time.get_time()
           << "| Edge Cut: " << edge_cut
           << endl;
    appendToFile(result.str());
}

void FennelPartitioner::calculate_edge_cut() {

}