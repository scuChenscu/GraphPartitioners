#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
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


    // 已经分配的边数
    atomic<size_t> assigned_edges;
    // 每个分区顶点的数目
    vector<size_t> num_vertices_each_partition;


    // TODO 邻接表，每个顶点的邻居顶点；该数据应该来自Graph
    unordered_map<vid_t, set<vid_t>> adjacency_list;
    // 每个分区的边的数量
    vector<size_t> occupied;
    // 每个顶点所属的分区，每个顶点可能属于多个分区
    vector<dense_bitset> is_mirrors;
    // TODO 需要一个结构来记录每个边所属的分区

    // 计算副本
    void calculate_replication_factor();

    virtual // 计算负载均衡因子
    void calculate_alpha();

    void calculate_indices() override;

    void print_indices() override;
};