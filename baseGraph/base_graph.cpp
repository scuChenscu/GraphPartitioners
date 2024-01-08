//
// Created by 陈键淞 on 2024/1/8.
//

#include "base_graph.hpp"
#include "../algorithms/ne/ne.hpp"
#include "../algorithms/dbh/dbh.hpp"
#include "../algorithms/hdrf/hdrf.hpp"
#include "../algorithms/ldg/ldg.hpp"
#include "../algorithms/fennel/fennel.hpp"
#include "../algorithms/model4/model4.hpp"

BaseGraph::BaseGraph(const string &graph_name) {
    this->graph_name = graph_name;

    ifstream fin(binary_edgelist_name(graph_name),
                 std::ios::binary | std::ios::ate);
    // tellp 用于返回写入位置，
    // tellg 则用于返回读取位置也代表着输入流的大小
    auto filesize = fin.tellg();
    fin.seekg(0, std::ios::beg);

    //最大的下标+1
    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));
    LOG(INFO) << "file size: " << filesize
              << " | num_vertices: " << num_vertices
              << " | num_edges: " << num_edges
              << endl;

    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
             filesize);

    edges.resize(num_edges);
    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    // 初始化的时候构造图
    adj_out.build(edges);
    // 存储反向边
    adj_in.build_reverse(edges);
    construct_adj_list();
    fin.close();
}

// 实现析构函数
//BaseGraph::~BaseGraph() {
//    edges.clear();
//    vertices.clear();
//    adjacency_list.clear();
//     delete adj_out;
//     delete adj_in;
//}

void BaseGraph::construct_adj_list() {
    LOG(INFO) << "Construct adjacency list..." << endl;
    // 遍历边集，建立每个顶点的邻居集合
    for (auto &edge: edges) {
        if (adjacency_list.contains(edge.first)) {
            adjacency_list.find(edge.first)->second.insert(edge.second);
        } else {
            std::set < vid_t > set;
            set.insert(edge.second);
            adjacency_list[edge.first] = set;
        }
        if (adjacency_list.contains(edge.second)) {
            adjacency_list.find(edge.second)->second.insert(edge.first);
        } else {
            std::set < vid_t > set;
            set.insert(edge.first);
            adjacency_list[edge.second] = set;
        }
    }
}

// TODO 因为每次执行图分区算法会依赖很多的数据结构，只能通过重新创建partitioner的方式来实现
void BaseGraph::partition() {
    // 所有的partitioner共享graph，即在一张图上面执行所有的partition算法，减少重复操作
    LOG(INFO) << "Execute algorithms on: " << graph_name << endl;
    stringstream graph;
    graph << "Graph Name: " << graph_name << endl;
    appendToFile(graph.str());
    int size = sizeof(algorithms) / sizeof(algorithms[0]);
    int partition = sizeof(partitions) / sizeof(partitions[0]);
    // TODO 这里就频繁地创建和销毁partitioner，是否可以优化
    for (int j = 0; j < partition; j++) {
        size_t num_partition = partitions[j];
        stringstream part;
        part << "Partition: " << num_partition << endl;
        LOG(INFO) << part.str();
        appendToFile(part.str());
        for (int i = 0; i < size; i++) {
            string algorithm = algorithms[i];
            LOG(INFO) << "Initialize " << algorithm << " partitioner: " << endl;
            Partitioner* partitioner;
            if (algorithm == "ne") {
                partitioner = new NePartitioner(*this, graph_name, algorithm, num_partition);
            } else if (algorithm == "model4") {
                partitioner = new Model4Partitioner(*this, graph_name, algorithm, num_partition);
            } else if (algorithm == "dbh") {
                partitioner = new DbhPartitioner(*this, graph_name, algorithm, num_partition, memory_size);
            } else if (algorithm == "hdrf") {
                partitioner = new HdrfPartitioner(*this, graph_name, algorithm, num_partition, memory_size,
                                                  balance_ratio, lambda, isShuffle);
            } else if (algorithm == "ldg") {
                partitioner = new LdgPartitioner(*this, graph_name, algorithm, num_partition, memory_size, isShuffle);
            } else if (algorithm == "fennel") {
                partitioner = new FennelPartitioner(*this, graph_name, algorithm, num_partition, memory_size,
                                                    isShuffle);
            } else {
                LOG(ERROR) << "Unknown algorithm: " << algorithm;
                continue;
            }
            LOG(INFO) << "Execute " << algorithm << " on: " << graph_name << endl;
            partitioner->split();
            LOG(INFO) << "Finish " << algorithm << " on: " << graph_name << endl;
            partitioner->calculate_indices();
            // LOG(INFO) << "Try to delete " << algorithm << " partitioner" << endl;
            // TODO 这里显式调用partitioner会报错
            // delete partitioner;
            // 日志打印是异步
            // LOG(INFO) << "Delete " << algorithm << " partitioner successfully" << endl;
        }
    }
}