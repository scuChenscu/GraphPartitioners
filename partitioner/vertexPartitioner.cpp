#pragma once
#include "vertexPartitioner.hpp"

VertexPartitioner::VertexPartitioner(BaseGraph &baseGraph, const string& algorithm, const size_t num_partitions) : Partitioner(baseGraph, algorithm, num_partitions) {
    assigned_vertices = 0;
    num_vertices_each_partition.resize(num_partitions, 0);
    vertex_partition.resize(num_vertices);
    partition_vertices.resize(num_partitions, dense_bitset(num_vertices));
    capacity = num_vertices / num_partitions + 1;
}

void VertexPartitioner::assigned_vertex(vid_t vid, size_t partition) {
    assigned_vertices++;
    num_vertices_each_partition[partition]++;
    vertex_partition[vid] = partition;
    partition_vertices[partition].set_bit_unsync(vid);
}

void VertexPartitioner::calculate_indices() {
    LOG(INFO) << "Calculating indices..." << endl;
    calculate_edge_cut();
    calculate_load_balance();
    print_indices();
    print_partition_vertices();
}

void VertexPartitioner::print_partition_vertices() {
    repv(j, num_partitions) {
        LOG(INFO) << " Partition " << j << " Vertex Count: " << num_vertices_each_partition[j];
    }
}

void VertexPartitioner::calculate_edge_cut() {
    LOG(INFO) << "Calculating edge cut..." << endl;
    // 边的两个顶点不在同一个分区
    for (auto edge : edges) {
        if (vertex_partition[edge.first] != vertex_partition[edge.second]) {
            edge_cut++;
        }
    }
    edge_cut_rate = double(edge_cut) / (double)num_edges;
}

void VertexPartitioner::calculate_load_balance() {
    LOG(INFO) << "Calculating load balance..." << endl;
    max_vertex = *max_element(num_vertices_each_partition.begin(), num_vertices_each_partition.end()); // 获取最大值
    min_vertex = *min_element(num_vertices_each_partition.begin(), num_vertices_each_partition.end());

    load_balance = (double) max_vertex * (double) num_partitions / (double) num_vertices;
}


void VertexPartitioner::print_indices() {
    stringstream result;
    result << fixed << setprecision(5);
    result << "Cost Time: " << total_time.get_time()
           << " | Edge Cut: " << edge_cut
           << " | Edge Cut Rate: " << edge_cut_rate
           << " | Load balance: " << load_balance
           << " | Max vertex: " << max_vertex
           << " | Min vertex: " << min_vertex
           << " | Edges: " << num_edges
           << " | Vertices: " << num_vertices
           << endl;
    LOG(INFO) << result.str();
    appendToFile(result.str());
}
