#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "../baseGraph/base_graph.hpp"
#include "../utils/util.hpp"

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
    // 需要执行分区算法的图
    BaseGraph graph;
public:
    explicit Partitioner(const BaseGraph& baseGraph) : graph(nullptr){
        graph = baseGraph;
    };

    // 分区算法
    void split() {};

    // 配置边分区、顶点分区结果输出文件名
    void config_output_files(const string &graph_name, const string &algorithm, int num_partition) {
        edge_ofstream.open(edge_partition_filename(graph_name, algorithm, num_partition));
        vertex_ofstream.open(vertex_partition_filename(graph_name, algorithm, num_partition));
    }

    // 输出顶点分区结果：vertex_id partition_id
    void save_vertex(vid_t vertex_id, int partition_id) {
        vertex_ofstream << vertex_id << " " << partition_id << endl;
    }

    // 输出边分区结果：from to partition_id
    void save_edge(vid_t from, vid_t to, int partition_id) {
        edge_ofstream << from << " " << to << " " << partition_id << endl;
    }

    ~Partitioner() {
        vertex_ofstream.close();
        edge_ofstream.close();
    }
};
