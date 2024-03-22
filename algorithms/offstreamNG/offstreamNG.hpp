#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <unordered_map>
#include <set>
#include <queue>
#include "../../utils/dense_bitset.hpp"
#include "../../utils/graph.hpp"
#include "../../utils/min_heap.hpp"
#include "../../partitioner/edgePartitioner.hpp"
#include "../../utils/util.hpp"
#include "../../baseGraph/base_graph.hpp"

using namespace std;
/* Neighbor Expansion (NE) */
class OffstreamNGPartitioner : public EdgePartitioner {
private:
    const double BALANCE_RATIO = 1.00;

    unordered_map<vid_t, vid_t> indices; // new_vid, old_vid
    // TODO 用dense_bitset
    set<vid_t> v_set; // 重新索引时已经被处理的顶点
    
    queue<vid_t> v_queue;

    string input;

    size_t current_partition;
    double average_degree;
    size_t capacity;


    vector<vector<vid_t> > part_degrees;
    vector<int> balance_vertex_distribute;
    MinHeap<vid_t, vid_t> d; // 顶点的度

    MinHeap<vid_t, vid_t> min_heap;

    vector<vid_t> degrees;

    vector<int8_t> master;
    vector<dense_bitset> is_cores, is_boundaries;
    dense_bitset true_vids;


    //随机数生成器
    //std::random_device rd;
    mt19937 gen;
    //均匀分布区间
    uniform_int_distribution<vid_t> dis;


    vector<edge_t> off_part;
    vector<edge_t> stream_part;

    vector<vector<int>> vertex_partitions;
    vector<dense_bitset> vp_set;

    double ratio = 0.5;

    size_t check_edge(const edge_t *e);

    void assign_edge(size_t partition, vid_t from, vid_t to);

    void add_boundary(vid_t vid);

    // 根据算法定义，把顶点加入的核心集时，需要把它的所有边都加入到边集合中
    void occupy_vertex(vid_t vid, vid_t d);


    bool get_free_vertex(vid_t &vid);

    bool get_target_vertex(vid_t &vid);

    void assign_remaining();

    void assign_master();

    size_t count_mirrors();

    size_t leastLoad(const dense_bitset& bitmap);

public:
    OffstreamNGPartitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                  size_t num_partitions);
    void split() override;
    // 广度遍历，重新索引，用于将顶点分块
    void reindex();

};

