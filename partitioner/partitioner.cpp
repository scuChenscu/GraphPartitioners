#pragma once

#include "partitioner.hpp"

Partitioner::Partitioner(BaseGraph& baseGraph, const string& algorithm, const size_t num_partitions): graph(baseGraph), num_partitions(num_partitions), algorithm(algorithm){
    this->graph = baseGraph;
    // 浅拷贝
    // BaseGraph g = baseGraph;
    // BaseGraph *b = baseGraph;
    // TODO 改成引用传递
    // const BaseGraph& g2 = baseGraph;

    // baseGraph.adjacency_list = nullptr;

    this->graph_name = graph.graph_name;
    this-> num_partitions = num_partitions;
    this->edges = graph.edges;
    this->vertices = graph.vertices;
    this->num_vertices = graph.num_vertices;
    this->num_edges = graph.num_edges;
    this->algorithm = algorithm;

    this->adj_out = graph.adj_out;

    this->adj_in = graph.adj_in;
}
//
//Partitioner::~Partitioner() {
//    // LOG(INFO) << "Partitioner::~Partitioner()" << endl;
//     vertex_ofstream.close();
//     edge_ofstream.close();
//}

// 配置边分区、顶点分区结果输出文件名
void Partitioner::config_output_files() {
    edge_ofstream.open(edge_partition_filename(graph_name, algorithm, num_partitions));
    vertex_ofstream.open(vertex_partition_filename(graph_name, algorithm, num_partitions));
}

// 输出顶点分区结果：vertex_id partition_id
void Partitioner::save_vertex(vid_t vertex_id, int partition_id) {
    vertex_ofstream << vertex_id << " " << partition_id << endl;
}

// 输出边分区结果：from to partition_id
void Partitioner::save_edge(vid_t from, vid_t to, int partition_id) {
    edge_ofstream << from << " " << to << " " << partition_id << endl;
}

Partitioner::~Partitioner() {
    LOG(INFO) << "Partitioner::~Partitioner()" << endl;
}