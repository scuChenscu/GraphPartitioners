#pragma once

#include "dbh.hpp"

using namespace std;
// 以边为基本的划分单位
// Degree Based Hashing (DBH)[6]分割通过判断结点的度信息来切分结点分配边。对于幂律图来说低度结点的局部性很容易保持，高度结点因为关联太多结点如果将边全部分配在一个子图上不太可能，因此该算法尽最大可能保持低度结点的局部性。
// DBH也是把边划分到不同的分区，计算负载因子
DbhPartitioner::DbhPartitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                               size_t num_partitions, int memory_size) : EdgePartitioner(baseGraph, algorithm, num_partitions){
    config_output_files();
    ifstream fin(binary_edgelist_name(input),
                 std::ios::binary | std::ios::ate);
    //degree file
    filesize = fin.tellg();
    fin.seekg(0, std::ios::beg);
    this->memory_size = memory_size;
    degrees.resize(num_vertices);
    ifstream degree_file(degree_name(input), std::ios::binary);
    degree_file.read((char *) &degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();

    num_batches = (filesize / ((std::size_t) memory_size * 1024 * 1024)) + 1;
    num_edges_per_batch = (num_edges / num_batches) + 1;

    true_vids.resize(num_vertices);
    is_mirrors.assign(num_vertices, dense_bitset(num_partitions));
    counter.assign(num_partitions, 0);
    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);
}

void DbhPartitioner::read_and_do(string opt_name) {
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
    std::vector<edge_t> edges;
    auto num_edges_left = num_edges;
    for (uint32_t i = 0; i < num_batches; i++) {
        auto edges_per_batch = num_edges_per_batch < num_edges_left ? num_edges_per_batch : num_edges_left;
        edges.resize(edges_per_batch);
        fin.read((char *) &edges[0], sizeof(edge_t) * edges_per_batch);
        if (opt_name == "dbh") {
            batch_dbh(edges);
        } else if (opt_name == "node_assignment") {
            batch_node_assignment(edges);
        } else {
            LOG(ERROR) << "no valid opt function";
        }
        num_edges_left -= edges_per_batch;
    }
}

void DbhPartitioner::batch_dbh(vector<edge_t> &edges) {
    for (auto &e: edges) {
        vid_t w = degrees[e.first] <= degrees[e.second] ? e.first : e.second;
        int bucket = w % num_partitions;
        counter[bucket]++;
        // save_edge(e->first,e->second,bucket);
        is_mirrors[e.first].set_bit_unsync(bucket);
        is_mirrors[e.second].set_bit_unsync(bucket);
        ++part_degrees[e.first][bucket];
        ++part_degrees[e.second][bucket];
        true_vids.set_bit_unsync(e.first);
        true_vids.set_bit_unsync(e.second);
    }
}

void DbhPartitioner::batch_node_assignment(vector<edge_t> &edges) {
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
    }
}

void DbhPartitioner::split() {
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "DBH" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    total_time.start();
    for (auto &e: edges) {
        vid_t vid = degrees[e.first] <= degrees[e.second] ? e.first : e.second;
        size_t partition = vid % num_partitions;
        occupied[partition]++;
        is_mirrors[e.first].set_bit_unsync(partition);
        is_mirrors[e.second].set_bit_unsync(partition);
        ++part_degrees[e.first][partition];
        ++part_degrees[e.second][partition];
        true_vids.set_bit_unsync(e.first);
        true_vids.set_bit_unsync(e.second);
    }
    total_time.stop();
}
