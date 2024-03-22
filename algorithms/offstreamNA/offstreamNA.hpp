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
class OffstreamNAPartitioner : public EdgePartitioner {
private:
    const double BALANCE_RATIO = 1.00;

    unordered_map<vid_t, vid_t> indices; // new_vid, old_vid
    
    queue<vid_t> v_queue;

    string input;

    size_t current_partition;
    size_t capacity;


    vector<int> balance_vertex_distribute;

    MinHeap<vid_t, vid_t> min_heap;
    vector<vid_t> degrees;

    vector<int8_t> master;
    vector<dense_bitset> is_cores, is_boundaries;
    // dense_bitset true_vids;

    vector<edge_t> off_part;
    vector<edge_t> stream_part;

    // 需要维护的状态
    vector<edge_t> window;

    unordered_set<vid_t> edges_in_window;
    // 第一次访问边的时候更新
    vector<size_t> partial_degree;
    // 每个顶点在分区的邻居数目
    vector<vector<size_t>> vertex_partitions;
    // 每个顶点在分区是否有副本
    vector<dense_bitset> vp_set;
    // is_mirrors 每个分区的顶点副本
//    size_t max_load;
//    size_t min_load;

    const double epsilon = 1;

    //随机数生成器
    //std::random_device rd;
    mt19937 gen;
    //均匀分布区间
    uniform_int_distribution<vid_t> dis;

    uint64_t max_partition_load;

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
    double calculate_lb_score(size_t partition_id);
    double calculate_rf_score(vid_t v, vid_t u, size_t p);
    double calculate_cs_score(vid_t v, vid_t u, size_t p);
    void remove_edge_from_window();
    void add_edge_to_window(edge_t& edge);
    int find_max_score_partition(edge_t &e);

public:
    OffstreamNAPartitioner(BaseGraph& baseGraph, string  input, const string& algorithm,
                  size_t num_partitions);
    void split() override;
    // 广度遍历，重新索引，用于将顶点分块
    void reindex();

};

