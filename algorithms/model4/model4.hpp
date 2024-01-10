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
#include <thread>
#include "../../utils/dense_bitset.hpp"
#include "../../partitioner/edgePartitioner.hpp"
#include "../../utils/graph.hpp"
#include "../../utils/min_heap.hpp"
#include "../../utils/util.hpp"

using namespace std;

/* Neighbor Expansion (NE) */
class Model4Partitioner : public EdgePartitioner {
private:
    const double BALANCE_RATIO = 1.00;
    const double CAPACITY_RATIO = 0.50;

    vector<vid_t> indices; // new_vid, old_vid
    // TODO 记录每个原始顶点在indices的下表
    vector<size_t> reverse_indices;

    string input;
    
    double average_degree;
    size_t capacity;
    size_t num_vertices_each_cores;
    vector<vector<vid_t> > part_degrees;
    vector<size_t> balance_vertex_distribute;
    size_t cores;
    size_t avg_vertices_each_partition;
    // MinHeap<vid_t, vid_t> d; // 顶点的度
    // 存储边
    vector<edge_t> edges;
    // 图结构
    graph_t adj_out;
    graph_t adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    // 为每个分区维护一个min_heap
    vector<MinHeap<vid_t, vid_t> > min_heaps;
    size_t current_partition;
    //每个分区边的数量
    vector<size_t> occupied;
    // TODO 需要一个结构存储每个分区顶点的数目
    vector<size_t> num_vertices_in_partition;

    vector<vid_t> degrees;

    vector<int8_t> master;

    // is_cores和is_boundaries是每个分区独立的dense_bitset
    vector<dense_bitset> is_cores, is_boundaries;
    dense_bitset true_vids;
    vector<dense_bitset> is_mirrors;
    dense_bitset visited;

    //随机数生成器
    //std::random_device rd;
    mt19937 gen;
    //均匀分布区间
    uniform_int_distribution<vid_t> dis;
    uniform_int_distribution<vid_t> sub_dis;

    size_t check_edge(const edge_t *e);

    void assign_edge(size_t index, vid_t from, vid_t to);

    void add_boundary(vid_t vid);

    void sub_add_boundary(vid_t vid, size_t partition);


    // 根据算法定义，把顶点加入的核心集时，需要把它的所有边都加入到边集合中
    void occupy_vertex(vid_t vid, vid_t d);

    void sub_occupy_vertex(vid_t vid, vid_t d,  size_t index);


    bool get_free_vertex(vid_t &vid);

    bool sub_get_free_vertex(vid_t &vid, vid_t index);
    void assign_remaining();

    void assign_master();

    size_t count_mirrors();

    void sub_split(size_t partition_index);

public:
    Model4Partitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                      size_t num_partitions);

    void split();

    // 广度遍历，重新索引，用于将顶点分块
    void re_index();
};

