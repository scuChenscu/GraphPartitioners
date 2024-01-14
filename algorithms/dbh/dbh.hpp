#pragma once

#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include "../../utils/dense_bitset.hpp"
#include "../../utils/graph.hpp"
#include "../../utils/min_heap.hpp"
#include "../../partitioner/edgePartitioner.hpp"
#include "../../utils/util.hpp"
#include "../../baseGraph/base_graph.hpp"

using namespace std;

class DbhPartitioner : public EdgePartitioner {
private:
    size_t filesize;
    int memory_size;
    vector<vid_t> degrees;
    ifstream fin;
    uint32_t num_batches;
    uint32_t num_edges_per_batch;
    dense_bitset true_vids;
    // 记录每个分区的边数
    vector<size_t> counter;
    vector<vector<vid_t> > part_degrees;
    vector<int> balance_vertex_distribute;



protected:
    void read_and_do(string opt_name);

    void batch_dbh(vector<edge_t> &edges);

    void batch_node_assignment(vector<edge_t> &edges);

public:
    DbhPartitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                   size_t num_partitions, int memory_size);

    void split() override;
    // void calculate_replication_factor() override;
};
