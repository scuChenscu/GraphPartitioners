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
#include <algorithm>
#include "../../utils/dense_bitset.hpp"
#include "../../utils/graph.hpp"
#include "../../utils/min_heap.hpp"
#include "../../partitioner/edgePartitioner.hpp"
#include "../../utils/util.hpp"
#include "../../baseGraph/base_graph.hpp"

using namespace std;
/* Neighbor Expansion (NE) */
class DpqnePartitioner : public EdgePartitioner {
private:
    const double BALANCE_RATIO = 1.00;

    unordered_map<vid_t, vid_t> indices; // new_vid, old_vid
    
    queue<vid_t> v_queue;

    string input;

    size_t current_partition;
    size_t capacity;


    vector<vector<vid_t> > part_degrees;
    vector<int> balance_vertex_distribute;

    vector<unordered_set<vid_t>> degree2vertices;
    int number_of_vertices = 0;
    // int minIndex = 0;
    // 对应在degree2vertices中的度数
    unordered_map<vid_t, vid_t> vertex2degree;

    int threshold;
    int candidates; // 要计算的个数

    MinHeap<vid_t, vid_t> high_min_heap;

    size_t front_partition;
    double avg_factor;
    double front_factor;

    // MinHeap<vid_t, vid_t> min_heap;


    vector<vid_t> degrees;

    vector<int8_t> master;
    vector<dense_bitset> is_cores, is_boundaries;
    dense_bitset true_vids;

    //随机数生成器
    //std::random_device rd;
    mt19937 gen;
    //均匀分布区间
    uniform_int_distribution<vid_t> dis;

    size_t check_edge(const edge_t *e);

    void assign_edge(size_t partition, vid_t from, vid_t to);

    void add_boundary(vid_t vid);

    // 根据算法定义，把顶点加入的核心集时，需要把它的所有边都加入到边集合中
    void occupy_vertex(vid_t vid, vid_t d);
    void decrease_key(vid_t vid);


    bool get_free_vertex(vid_t &vid);

    bool get_target_vertex(vid_t &vid);

    void assign_remaining();

    void assign_master();

    size_t count_mirrors();
    bool get_min_value(vid_t & degree, vid_t& vid, vid_t& position);
    void remove(vid_t vid, vid_t position);
public:
    DpqnePartitioner(BaseGraph& baseGraph, string  input, const string& algorithm,
                  size_t num_partitions);
    void split() override;
    // 广度遍历，重新索引，用于将顶点分块
    void reindex();

};

