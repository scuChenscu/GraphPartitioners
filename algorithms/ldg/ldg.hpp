#pragma once
#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include "../../partitioner/partitioner.hpp"
#include "../../utils/dense_bitset.hpp"
#include "../../utils/util.hpp"
#include "../../partitioner/vertexPartitioner.hpp"

using namespace std;

class LdgPartitioner : public VertexPartitioner {
private:
    vid_t num_vertices;
    size_t num_edges;
    size_t filesize;

    int p;
    ifstream fin;

    uint32_t num_batches;
    uint32_t num_edges_per_batch;

    vector<unordered_set<vid_t> > node2neis;
    unordered_set<vid_t> true_vids;
    vector<unordered_set<vid_t> > subg_vids;
    vector<int> balance_vertex_distribute;

protected:
    void read_and_do(string opt_name);

    void do_ldg();

    void batch_node_assignment(vector<edge_t> &edges);

    void addNeighbors(edge_t &edge);

    int intersection(unordered_set<vid_t> &nums1, unordered_set<vid_t> &nums2);

public:
    LdgPartitioner(const BaseGraph &baseGraph, const string& input, const string& algorithm, const size_t num_partitions, int memory_size, bool shuffle);

    void split();

    void calculate_edge_cut();
};
