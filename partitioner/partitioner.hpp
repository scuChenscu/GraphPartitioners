#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "../baseGraph/base_graph.hpp"
#include "../utils/util.hpp"
#include "../utils/dense_bitset.hpp"

using namespace std;

/**
 * Partitioner，分区器，所有图分区算法都继承Partitioner，实现split方法
 */

class Partitioner {
protected:
    // 分区算法执行耗时
    Timer total_time;
    // 输出文件流对象
    ofstream edge_ofstream;
    ofstream vertex_ofstream;
    // 持有引用
    BaseGraph& graph;
    size_t num_partitions;
    // TODO 每次分区会修改edges，已确认是拷贝，不会影响graph
    vector<edge_t> edges;
    vector<vid_t> vertices;
    // 度数
    vector<size_t> degrees;
    // 真实顶点数
    dense_bitset true_vids;
    string algorithm;
    graph_t& adj_out;
    graph_t& adj_in;
    // TODO 这个不能改，因为二进制文件是用vid_t类型写入
    vid_t num_vertices;
    size_t num_edges;
    string graph_name;

public:
    explicit Partitioner(BaseGraph& baseGraph, const string& algorithm, size_t num_partitions);

    virtual void split() {
        LOG(ERROR) << "Derived class has not override split function" << endl;
    };

    virtual void calculate_indices() {
        LOG(ERROR) << "Derived class has not override calculate_indices function" << endl;
    }

    virtual void print_indices() {
        LOG(ERROR) << "Derived class has not override print_indices function" << endl;
    }

    // 配置边分区、顶点分区结果输出文件名
    void config_output_files();

    // 输出顶点分区结果：vertex_id partition_id
    void save_vertex(vid_t vertex_id, int partition_id);

    // 输出边分区结果：from to partition_id
    void save_edge(vid_t from, vid_t to, int partition_id);

    ~Partitioner() {
        edge_ofstream.close();
        vertex_ofstream.close();
    }
};



