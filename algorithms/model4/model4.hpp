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
    // TODO 不应该用set，应该用vector
    // set<vid_t> v_set; // 重新索引时已经被处理的顶点

    unordered_map<vid_t, set<vid_t>> adj_list; // 邻接表
    string input;
    vid_t num_p_v;
    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p, bucket;
    double average_degree;
    size_t capacity;

    vector<vector<vid_t> > part_degrees;
    vector<int> balance_vertex_distribute;
    // MinHeap<vid_t, vid_t> d; // 顶点的度
    // 存储边
    vector<edge_t> edges;
    // 图结构
    graph_t adj_out;
    graph_t adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    // 为每个分区维护一个min_heap
    vector<MinHeap<vid_t, vid_t> > min_heaps;

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

    int check_edge(const edge_t *e) {
        rep (i, bucket) {
            auto &is_boundary = is_boundaries[i];
            if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
                occupied[i] < capacity) {
                return i;
            }
        }

        rep (i, bucket) {
            auto &is_core = is_cores[i], &is_boundary = is_boundaries[i];
            if ((is_core.get(e->first) || is_core.get(e->second)) &&
                occupied[i] < capacity) {
                if (is_core.get(e->first) && degrees[e->second] > average_degree)
                    continue;
                if (is_core.get(e->second) && degrees[e->first] > average_degree)
                    continue;
                is_boundary.set_bit(e->first);
                is_boundary.set_bit(e->second);
                return i;
            }
        }

        return p;
    }

    void assign_edge(int partition, vid_t from, vid_t to) {
        // save_edge(from, to, bucket);
        true_vids.set_bit_unsync(from);
        true_vids.set_bit_unsync(to);
        is_mirrors[from].set_bit_unsync(partition);
        is_mirrors[to].set_bit_unsync(partition);
        assigned_edges++;
        occupied[partition]++;
        degrees[from]--;
        degrees[to]--;
    }

    void add_boundary(vid_t vid) {
        // 获取到当前分区的核心集和边界集
        auto &is_core = is_cores[bucket];
        auto &is_boundary = is_boundaries[bucket];

        // 如果已经被加入到边界集，直接返回
        if (is_boundary.get(vid))
            return;
        // 把顶点加入到边界集
        is_boundary.set_bit_unsync(vid);

        // 如果顶点没有在核心集中，直接把顶点的度数据加入到最小堆
        // 因为在核心集中的顶点已经不需要处理了
        if (!is_core.get(vid)) {
            min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
        }

        //仅支持无向图，在计算neighbor的时候有向和无向会导致邻居的差别从而影响分割

        rep (direction, 2) {
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            // 遍历顶点的邻边
            for (size_t i = 0; i < neighbors.size();) {
                if (edges[neighbors[i].v].valid()) {
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    if (is_core.get(u)) { // 如果顶点在核心集中
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid);
                        // TODO 这个要修改
                        min_heap.decrease_key(vid); // 默认移除一条边
                        edges[neighbors[i].v].remove();
                        //TODO 交换到最后位置，然后长度减1
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else if (is_boundary.get(u) &&
                               occupied[bucket] < capacity) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                        // 将边加入边集
                        assign_edge(bucket, direction ? vid : u, direction ? u : vid);
                        min_heap.decrease_key(vid);
                        min_heap.decrease_key(u);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else
                        i++;
                } else {
                    //swap是pop的前提，先交换到最后位置然后把长度减1
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }

    void sub_add_boundary(vid_t vid, int cur_p) {
        // 获取到当前分区的核心集和边界集
        auto &is_core = is_cores[cur_p];
        auto &is_boundary = is_boundaries[cur_p];
        // 当从S\C引入顶点时，它的邻居顶点可能已经加入到边界集
        if (is_boundary.get(vid)) return;
        // 2. 如果顶点不在边界集，把顶点加入到边界集
        is_boundary.set_bit_unsync(vid);
        // 当我们引入新的顶点到边界集的时候，需要把它的度信息加入到最小堆
        if (!is_core.get(vid)) {
            min_heaps[cur_p].insert(adj_out[vid].size() + adj_in[vid].size(), vid);

        }

        // 3. 把顶点的在边界集的邻居顶点的边加入到当前边集
        rep (direction, 2) {
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            // 遍历顶点的邻边
            for (size_t i = 0; i < neighbors.size();) {
                if (edges[neighbors[i].v].valid()) {
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    if (is_core.get(u)) { // 如果顶点在核心集中
                        assign_edge(cur_p, direction ? vid : u,
                                    direction ? u : vid);
                        min_heaps[cur_p].decrease_key(vid); // 默认移除一条边
                        edges[neighbors[i].v].remove();
                        //TODO 交换到最后位置，然后长度减1
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else if (is_boundary.get(u) &&
                               occupied[cur_p] < capacity * CAPACITY_RATIO) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                        // 将边加入边集
                        assign_edge(cur_p, direction ? vid : u, direction ? u : vid);
                        min_heaps[cur_p].decrease_key(vid);
                        min_heaps[cur_p].decrease_key(u);
                        edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else
                        i++;
                } else {
                    //swap是pop的前提，先交换到最后位置然后把长度减1
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }


    // 根据算法定义，把顶点加入的核心集时，需要把它的所有边都加入到边集合中
    void occupy_vertex(vid_t vid, vid_t d) {
        CHECK(!is_cores[bucket].get(vid)) << "add " << vid << " to core again";
        // 核心集是vector<dense_bitset>，dense_bitset是一个稠密位图
        // 对位图的vid位置置1，表示vid被分配到bucket分区
        is_cores[bucket].set_bit_unsync(vid);
        // 如果顶点的度为0，不需要处理
        if (d == 0) return;

        // 这一步的意义是什么？
        add_boundary(vid);

        // 将顶点的邻居加入到边界集，这里需要过滤掉核心集中的顶点
        for (auto &i: adj_out[vid]) {
            if (edges[i.v].valid())
                add_boundary(edges[i.v].second);
        }
        adj_out[vid].clear();

        for (auto &i: adj_in[vid]) {
            if (edges[i.v].valid())
                add_boundary(edges[i.v].first);
        }
        adj_in[vid].clear();
    }

    void sub_occupy_vertex(vid_t vid, vid_t d, int cur_p) {
        CHECK(!is_cores[cur_p].get(vid)) << "add " << vid << " to core again";
        // 核心集是vector<dense_bitset>，dense_bitset是一个稠密位图
        // 对位图的vid位置置1，表示vid被分配到bucket分区

        // 1. 把顶点加入核心集和边界集，不需要考虑重复设置的场景
        is_cores[cur_p].set_bit_unsync(vid);
        // is_boundaries[cur_p].set_bit_unsync(vid);
        // 顶点的来源有两种情况，一是从V，二是从S\C
        // 如果顶点的度为0，不需要处理
        if (d == 0) return;
        sub_add_boundary(vid, cur_p);

        // 2. 遍历vid的邻居顶点，把不属于边界集的顶点加入到边界集
        for (auto &i: adj_out[vid]) {
            if (edges[i.v].valid()) { // 因为是边分区算法，所以我们不能引入重复的边，只需要考虑未分配的边
                sub_add_boundary(edges[i.v].second, cur_p);
            }
        }

        adj_out[vid].clear();

        for (auto &i: adj_in[vid]) {
            if (edges[i.v].valid())
                sub_add_boundary(edges[i.v].first, cur_p);
        }
        adj_in[vid].clear();
    }


    bool get_free_vertex(vid_t &vid) {
        //随机选择一个节点
        vid = dis(gen);
        vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
        //TODO 什么叫已经超出平衡范围
        //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
        while (count < num_vertices &&
               (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                2 * average_degree ||
                is_cores[bucket].get(vid))) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    bool sub_get_free_vertex(vid_t &vid, vid_t cur_p) {
        //随机选择一个节点
        vid_t min = cur_p * num_p_v;
        vid_t index = sub_dis(gen); // 生成 1000 / 10 = 100；0-99
        vid_t max = (cur_p + 1) * num_p_v;
        // 选择数据范围在[min, max]

        vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点

        vid = indices[min + index];

        //TODO 什么叫已经超出平衡范围
        //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
        while (count < num_p_v &&
               (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                2 * average_degree ||
                is_cores[cur_p].get(vid))) { // 此时候选集跟核心集是一致的，只需要判断一个即可
            // 应该是一个链式寻找的过程
            index = (index + ++count) % num_p_v;
            vid = indices[min + index];
        }
        if (count == num_p_v)
            return false;
        return true;
    }
    void assign_remaining();

    void assign_master();

    size_t count_mirrors();

    void sub_split(int i);

public:
    Model4Partitioner(const BaseGraph& baseGraph, const string& input, const string& algorithm,
                      size_t num_partitions);

    void split();

    // 广度遍历，重新索引，用于将顶点分块
    void re_index() {
        queue<vid_t> v_queue;
        auto start = std::chrono::high_resolution_clock::now(); // 记录开始时间
        // 随机选择顶点，进行广度遍历，重新索引
        vid_t index = 0;
        vid_t vid = dis(gen);
        // 基于该顶点进行深度遍历，对每个顶点重新索引
        v_queue.push(vid);
        while (!v_queue.empty()) {
            // LOG(INFO) << index;
            vid_t v = v_queue.front();
            v_queue.pop();
            if (visited.get(v)) {
                continue;
            }
            visited.set_bit_unsync(v);
            // 将v加入到indices,重新索引
            indices[index++] = v;

            // 获取v的邻居顶点
            set < vid_t > neighbor_set = adj_list.find(v)->second;
            // 将neighbor_set加入v_queue和v_set中
            for (auto &i: neighbor_set) {
                v_queue.push(i);
            }
        }
        auto end = std::chrono::high_resolution_clock::now(); // 记录结束时间
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // 计算时间差
        LOG(INFO) << "re_index time: " << duration.count() << "ms" << endl;
    }


};

