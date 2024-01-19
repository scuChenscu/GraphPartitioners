//
// Created by 陈键淞 on 2024/1/8.
//

#include "base_graph.hpp"
#include "../algorithms/ne/ne.hpp"
#include "../algorithms/dbh/dbh.hpp"
#include "../algorithms/hdrf/hdrf.hpp"
#include "../algorithms/ldg/ldg.hpp"
#include "../algorithms/fennel/fennel.hpp"
#include "../algorithms/greedy/greedy.hpp"
#include "../algorithms/rand/rand.hpp"
#include "../algorithms/model4/model4.hpp"
#include "../algorithms/model5/model5.hpp"
#include "../algorithms/model6/model6.hpp"

using namespace std;

BaseGraph::BaseGraph(const string& graph_name) {
    this->graph_name = graph_name;

    ifstream fin(binary_edgelist_name(graph_name),
                 ios::binary | ios::ate);
    // tellp 用于返回写入位置，
    // tellg 则用于返回读取位置也代表着输入流的大小
    auto filesize = fin.tellg();
    fin.seekg(0, ios::beg);

    //最大的下标+1
    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));
    LOG(INFO) << "File size: " << filesize
              << " | num_vertices: " << num_vertices
              << " | num_edges: " << num_edges
              << endl;

    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
             filesize);

    edges.resize(num_edges);
    degrees.resize(num_vertices, 0);
    true_vids.resize(num_vertices);
    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    // 初始化的时候构造图
     adj_out.build(edges);
     // 存储反向边
     adj_in.build_reverse(edges);
    gen.seed(DEFAULT_SEED);

    dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    visited = dense_bitset(num_vertices);
    indices.resize(num_vertices);
    reverse_indices.resize(num_vertices);

    construct_adjacency_list();

    // 重新索引
    if (REINDEX) {
        LOG(INFO) << "re_index" << endl;
        re_index();
    } else {

    }


}
// TODO 建立CSR
void BaseGraph::construct_adjacency_list() {
    LOG(INFO) << "Construct adjacency list..." << endl;
    // 遍历边集，建立每个顶点的邻居集合
    for (auto &edge: edges) {
        // 计算顶点度数
        degrees[edge.first]++;
        degrees[edge.second]++;
        true_vids.set_bit_unsync(edge.first);
        true_vids.set_bit_unsync(edge.second);

        if (adjacency_list.contains(edge.first)) {
            // LOG(INFO) << edge.first;
            adjacency_list.find(edge.first)->second.insert(edge.second);
        } else {
            set <vid_t> set;
            set.insert(edge.second);
            adjacency_list[edge.first] = set;
        }
        if (adjacency_list.contains(edge.second)) {
            // LOG(INFO) << edge.second;
            adjacency_list.find(edge.second)->second.insert(edge.first);
        } else {
            set <vid_t> set;
            set.insert(edge.first);
            adjacency_list[edge.second] = set;
        }
    }
}

void BaseGraph::partition() {
    // 所有的partitioner共享graph，即在一张图上面执行所有的partition算法，减少重复操作
    LOG(INFO) << "Start partition on " << graph_name << endl;
    stringstream ss;
    ss << "Graph Name: " << graph_name << endl;
    appendToFile(ss.str());
    ss.str("");  // 清空当前字符串内容
    ss.clear();   // 清空错误状态标志
    vector<Partitioner*> partitioners;
    for (auto num_partitions : partitions) {
        Partitioner* partitioner;
//        ss << "Number Partitions: " << num_partitions << endl;
//        appendToFile(ss.str());
//        ss.str("");  // 清空当前字符串内容
//        ss.clear();   // 清空错误状态标志
        for (auto& algorithm : algorithms) {
            // TODO this是一个指针
            if (algorithm == "ne") {
                partitioner = new NePartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "model4") {
                partitioner = new Model4Partitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "model5") {
                // TODO 三个核心参数：cores、balance_ratio、capacity_ratio
                if (SELF) {
                    partitioner = new Model5Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,1);
                    partitioners.push_back(partitioner);
                    continue;
                }
                // 三层循环cores、balance_ratio、capacity_ratio
                for(size_t cores = 1; cores <= MAX_CORES; cores++) {
                    for (auto ours_balance_ratio : OURS_BALANCE_RATIOS) {
                        for (auto ours_capacity_ratio : OURS_CAPACITY_RATIOS) {
                            partitioner = new Model5Partitioner(*this, graph_name, algorithm, num_partitions, ours_balance_ratio, ours_capacity_ratio, cores);
                            partitioners.push_back(partitioner);
                        }
                    }
                }
            } else if (algorithm == "model6") {
                partitioner = new Model6Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            } else if (algorithm == "dbh") {
                partitioner = new DbhPartitioner(*this, graph_name, algorithm, num_partitions, memory_size);
                partitioners.push_back(partitioner);
            } else if (algorithm == "hdrf") {
                partitioner = new HdrfPartitioner(*this, graph_name, algorithm, num_partitions, memory_size,
                                                  balance_ratio, lambda, isShuffle);
                partitioners.push_back(partitioner);
            } else if (algorithm == "ldg") {
                partitioner = new LdgPartitioner(*this, graph_name, algorithm, num_partitions, memory_size, isShuffle);
                partitioners.push_back(partitioner);
            } else if (algorithm == "fennel") {
                partitioner = new FennelPartitioner(*this, graph_name, algorithm, num_partitions, memory_size,
                                                    isShuffle);
                partitioners.push_back(partitioner);
            } else if (algorithm == "greedy") {
                partitioner = new GreedyPartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "rand") {
                partitioner = new RandPartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            }
            else {
                LOG(ERROR) << "Unknown algorithm: " << algorithm;
                continue;
            }
        }
    }
    for (auto partitioner : partitioners) {
        // LOG(INFO) << "Execute " << algorithm << " on: " << graph_name << endl;
        partitioner->split();
        // LOG(INFO) << "Finish " << algorithm << " on: " << graph_name << endl;
        partitioner->calculate_indices();
        delete partitioner;
        // 将指针设置为 nullptr，以防止悬垂指针
    }
}

void BaseGraph::re_index() {
    // LOG(INFO) << adjacency_list.size();
    queue<vid_t> v_queue;
    auto start = std::chrono::high_resolution_clock::now(); // 记录开始时间
    // 随机选择顶点，进行广度遍历，重新索引
    vid_t index = 0;
    vid_t vid = dis(gen);
    // 基于该顶点进行深度遍历，对每个顶点重新索引
    // TODO 该顶点可能没有邻居
    v_queue.push(vid);
    while (!v_queue.empty()) {
        // LOG(INFO) << index;
        vid_t v = v_queue.front();
        v_queue.pop();
        if (visited.get(v)) {
            continue;
        }
        visited.set_bit_unsync(v);
        // 将v加入到indices,重新索引
        reverse_indices[v] = index; // vid所在indices的下标为index
        indices[index++] = v;

        // 获取v的邻居顶点
        if (!adjacency_list.contains(v)) continue;
        // LOG(INFO) << v;
        set <vid_t> neighbor_set = adjacency_list.find(v)->second;
        // 将neighbor_set加入v_queue和v_set中
        for (auto &i: neighbor_set) {
            v_queue.push(i);
        }
    }
    auto end = std::chrono::high_resolution_clock::now(); // 记录结束时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // 计算时间差
    LOG(INFO) << "re_index time: " << duration.count() << "ms" << endl;
}