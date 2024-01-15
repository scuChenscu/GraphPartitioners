#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "partitioner.hpp"
#include "../baseGraph/base_graph.hpp"
#include "../utils/util.hpp"


class VertexPartitioner : public Partitioner {
public:
    explicit VertexPartitioner(BaseGraph &baseGraph,const string& algorithm, size_t num_partitions);
    void calculate_indices() override;
    void print_indices() override;
protected:
    // 边割数
    size_t edge_cut;
    // 边割率
    double edge_cut_rate;
    // 负载均衡
    double load_balance;
    // 已经分配的顶点数
    atomic<size_t> assigned_vertices;
    // 每个分区顶点的数目
    vector<size_t> num_vertices_each_partition;
    // 每个顶点所在的分区
    vector<size_t>  vertex_partition;
    // 每个分区的顶点
    vector<dense_bitset> partition_vertices;
    // 每个分区的理论顶点数
    size_t capacity;
    size_t max_vertex;
    size_t min_vertex;

    void assigned_vertex(vid_t vid, size_t partition);
    void calculate_edge_cut();
    void calculate_load_balance();
    void print_partition_vertices();
};