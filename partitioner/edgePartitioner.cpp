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
    LOG(INFO) << "Calculating replication factor..." << endl;
    for (auto &is_mirror: is_mirrors) {
        replicas += is_mirror.popcount();
    }
    replication_factor = (double) replicas / (double) true_vids.popcount();
}

// 计算负载均衡因子
void EdgePartitioner::calculate_alpha() {
    LOG(INFO) << "Calculating alpha..." << endl;
    max_edge = *max_element(occupied.begin(), occupied.end()); // 获取最大值
    min_edge = *min_element(occupied.begin(), occupied.end());

    alpha = (double) max_edge * (double) num_partitions / (double) graph.num_edges;
}

void EdgePartitioner::calculate_indices() {
    LOG(INFO) << "Calculating indices..." << endl;
    CHECK_EQ(assigned_edges, num_edges);
    calculate_replication_factor();
    calculate_alpha();
    print_indices();
    print_partition_edges();
}

void EdgePartitioner::print_partition_edges() {
    repv(j, num_partitions) {
        LOG(INFO) << "Partition " << j << " Edge Count: " << occupied[j];
    }
}

void EdgePartitioner::print_indices() {
    stringstream result;
    result << fixed << setprecision(5);
    result << "Cost Time: " << total_time.get_time()
           << " | Replication Factor: " << replication_factor
           << " | Alpha: " << alpha
           << " | Replicas: " << replicas
           << " | Partition: " << num_partitions
           << " | Max Edge: " << max_edge
           << " | Min Edge: " << min_edge
           // << " | Balance Ratio: " << balance_ratio
           << " | Avg Edge: " << num_edges / num_partitions
           << " | Edges: " << num_edges
           << " | Vertices: " << num_vertices
           << " | True Vertices: " << true_vids.popcount()
           << endl;
    LOG(INFO) << result.str();
    appendToFile(result.str());
}
