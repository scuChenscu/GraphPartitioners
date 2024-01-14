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
#include <bitset>
#include <set>
#include <queue>
#include <thread>
#include "../../utils/dense_bitset.hpp"
#include "../../partitioner/edgePartitioner.hpp"
#include "../../utils/graph.hpp"
#include "../../utils/min_heap.hpp"
#include "../../utils/util.hpp"

using namespace std;

class Model5Partitioner : public EdgePartitioner {
private:
    size_t cores;
    double capacity_ratio;
    // new_vid, old_vid
    vector<vid_t> indices;
    // TODO 记录每个原始顶点在indices的下标
    vector<size_t> reverse_indices;
    vector<size_t> v_lock;
    string input;
    double average_degree;
    dense_bitset dirty_vertices;
    size_t capacity;
    size_t num_vertices_each_cores;
    vector<vector<vid_t>> part_degrees;
    vector<size_t> balance_vertex_distribute;

    graph_t adj_directed;

    MinHeap<vid_t, vid_t> min_heap;
    // 为每个分区维护一个min_heap
    vector<MinHeap<vid_t, vid_t>> min_heaps;
    size_t current_partition;

    // 每条边的分区
    vector<size_t> edge_partition;

    vector<vid_t> degrees;

    vector<int8_t> master;

    // is_cores和is_boundaries是每个分区独立的dense_bitset
    vector<dense_bitset> is_cores, is_boundaries;
    dense_bitset true_vids;

    vector<dense_bitset> reverse_is_mirrors;

    dense_bitset visited;
    // 记录已经被使用过的边
    dense_bitset assigned;
    // 顶点锁，保证只能有一个线程持有某个顶点
    dense_bitset vertex_lock;
    //随机数生成器
    //std::random_device rd;
    mt19937 gen;
    //均匀分布区间
    uniform_int_distribution<vid_t> dis;
    uniform_int_distribution<vid_t> sub_dis;

    vector<unordered_map<size_t, edge_t*>> vertex_adjacent_edges;

    size_t check_edge(const edge_t *e);

    void assign_edge(size_t index, vid_t from, vid_t to);

    void add_boundary(vid_t vid);

    void sub_add_boundary(vid_t vid, size_t partition);


    bool get_free_vertex(vid_t &vid);

    bool sub_get_free_vertex(vid_t &vid, vid_t partition);

    void assign_remaining();

    void assign_master();

    // size_t count_mirrors();

    void sub_split(size_t partition);

public:
    Model5Partitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                      size_t num_partitions, double balance_ratio, double capacity_ratio, size_t cores);

    void build_vertex_adjacent_edges();

    void split() override;

    // 广度遍历，重新索引，用于将顶点分块
    void re_index();

    bool acquire_vertex(size_t vid) {
        bool success = __sync_bool_compare_and_swap(&v_lock[vid], 0, 1);
        std::thread::id currentThreadId = std::this_thread::get_id();
        if (success) {
            // LOG(INFO) << "Thread ID: " << currentThreadId << " acquires vid: " << vid << std::endl;
        } else {
            // LOG(INFO) << "Thread ID: " << currentThreadId << " fails to acquire vid: " << vid << std::endl;
        }
        return success;
    }

    bool release_vertex(size_t vid) {
        bool success = __sync_bool_compare_and_swap(&v_lock[vid], 1, 0);
        if (success) {
            // LOG(INFO) << "Thread ID: " << std::this_thread::get_id() << " releases vid: " << vid << std::endl;
        } else {
            // LOG(ERROR) << "Thread ID: " << std::this_thread::get_id() << " fails to release vid: " << vid << std::endl;
        }
        return success;
    }

    void calculate_replication_factor() override;
};

