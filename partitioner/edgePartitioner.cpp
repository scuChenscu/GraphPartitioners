#pragma once

#include "edgePartitioner.hpp"

EdgePartitioner::EdgePartitioner(BaseGraph &baseGraph, const std::string &algorithm,
                                 const size_t num_partitions) : Partitioner(baseGraph, algorithm, num_partitions) {
    this->adjacency_list = baseGraph.adjacency_list;
    this->occupied.assign(num_partitions, 0);
    this->is_mirrors.assign(num_vertices, dense_bitset(num_partitions));
    replication_factor = 0.0;
    alpha = 0.0;
    replicas = 0;
    balance_ratio = 0.0;
    assigned_edges = 0;
}

EdgePartitioner::~EdgePartitioner() = default;

// 计算副本
void EdgePartitioner::calculate_replication_factor() {
    for (auto &is_mirror: is_mirrors) {
        replicas += is_mirror.popcount();
    }
    replication_factor = (double) replicas / (double) graph.num_vertices;
}

// 计算负载均衡因子
void EdgePartitioner::calculate_alpha() {
    max_edge = *max_element(occupied.begin(), occupied.end()); // 获取最大值
    min_edge = *min_element(occupied.begin(), occupied.end());

    alpha = (double) max_edge * (double) num_partitions / (double) graph.num_edges;
}

void EdgePartitioner::calculate_indices() {
    calculate_replication_factor();
    calculate_alpha();
    print_indices();
}

void EdgePartitioner::print_indices() {
    stringstream result;
    result << "Cost Time: " << total_time.get_time()
           << " | Replication Factor: " << replication_factor
           << " | Alpha: " << alpha
           << " | Replicas: " << replicas
           << " | Partition: " << num_partitions
           << " | Max Edge: " << max_edge
           << " | Min Edge: " << min_edge
           << " | Balance Ratio: " << balance_ratio
           << " | Avg Edge: " << num_edges / num_partitions
           << " | Edges: " << num_edges
           << " | Vertices: " << num_vertices
           << endl;
    appendToFile(result.str());

    LOG(INFO) << result.str() << endl;
}
