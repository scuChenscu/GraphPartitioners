#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <atomic>
#include "partitioner.hpp"
#include "../baseGraph/base_graph.hpp"
#include "../utils/dense_bitset.hpp"
#include "../utils/util.hpp"

class EdgePartitioner : public Partitioner {
public:
     EdgePartitioner(BaseGraph& baseGraph, const string& algorithm, size_t num_partitions);
     ~EdgePartitioner();
protected:
    // 复制因子
    double replication_factor;
    // 负载均衡因子
    double alpha;
    // 顶点总副本数
    size_t replicas;

    size_t max_edge = numeric_limits<size_t>::min();

    size_t min_edge = numeric_limits<size_t>::max();

    // 分区边负载上限调节因子
    double balance_ratio;

    size_t max_degree = 0;
    size_t min_degree = UINT16_MAX;
    size_t avg_degree;

    // 已经分配的边数
    atomic<size_t> assigned_edges;
    // 每个分区顶点的数目
    vector<size_t> num_vertices_each_partition;

    // 每个分区的边的数量
    vector<size_t> occupied;
    // 每个顶点所属的分区，每个顶点可能属于多个分区
    vector<dense_bitset> is_mirrors;
    // TODO 需要一个结构来记录每个边所属的分区

    uint64_t min_load = UINT64_MAX;
    uint64_t max_load = 0;

    dense_bitset visited;
    vector<size_t> indices;
    vector<size_t> reverse_indices;

    // 计算副本
    virtual void calculate_replication_factor();
    // 计算负载均衡因子
    virtual void calculate_alpha();

    void calculate_indices() override;

    void print_indices() override;

    void print_partition_edges();

    void update_min_max_load(int max_p) {
        auto &load = occupied[max_p];
        if (load > max_load) max_load = load;
        min_load = *min_element(occupied.begin(), occupied.end());
    }
};