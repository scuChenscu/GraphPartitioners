#pragma once


#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_set>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include "../../partitioner/vertexPartitioner.hpp"
#include "../../utils/dense_bitset.hpp"
#include "../../utils/util.hpp"

using namespace std;

class FennelPartitioner : public VertexPartitioner {
private:
    vid_t num_vertices;
    size_t num_edges;
    size_t filesize;

    int p;
    ifstream fin;

    uint32_t num_batches;
    uint32_t num_edges_per_batch;

    vector<unordered_set<vid_t>> node2neis;
    // 用于过滤重复vid, 得到unique vids
    unordered_set<vid_t> true_vids;
    vector<unordered_set<vid_t>> subg_vids;
    vector<int> balance_vertex_distribute;

protected:
    void read_and_do(string opt_name);

    void do_fennel();

    void batch_node_assignment(vector<edge_t> &edges);

    void addNeighbors(edge_t &edge);

    int intersection(unordered_set<vid_t> &nums1, unordered_set<vid_t> &nums2);

public:
    FennelPartitioner(BaseGraph& baseGraph, const string& input, const string& algorithm, const size_t num_partitions, int memory_size, bool shuffle);

    void split();

    void calculate_edge_cut();
};
