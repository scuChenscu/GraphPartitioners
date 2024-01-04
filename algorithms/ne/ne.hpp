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
#include "../../partitioner/partitioner.hpp"
#include "../../utils/graph.hpp"
#include "../../utils/min_heap.hpp"
#include "../../partitioner/partitioner.hpp"
#include "../../utils/util.hpp"


using namespace std;

/* Neighbor Expansion (NE) */
class NePartitioner : public Partitioner {
private:
    const double BALANCE_RATIO = 1.00;

    unordered_map<vid_t, vid_t> indices; // new_vid, old_vid
    set<vid_t> v_set; // 重新索引时已经被处理的顶点

    unordered_map<vid_t, set<vid_t>> adj_list; // 邻接表

    queue<vid_t> v_queue;

    std::string input;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p, bucket;
    double average_degree;
    size_t capacity;

    vector<vector<vid_t> > part_degrees;
    vector<int> balance_vertex_distribute;
    MinHeap<vid_t, vid_t> d; // 顶点的度
    // 存储边
    std::vector<edge_t> edges;
    // 图结构
    graph_t adj_out, adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    //每个分区边的数量
    std::vector<size_t> occupied;
    // TODO 需要一个结构存储每个分区顶点的数目
    vector<size_t> num_vertices_in_partition;

    std::vector<vid_t> degrees;

    std::vector<int8_t> master;
    std::vector<dense_bitset> is_cores, is_boundarys;
    dense_bitset true_vids;
    vector<dense_bitset> is_mirrors;

    //随机数生成器
    //std::random_device rd;
    std::mt19937 gen;
    //均匀分布区间
    std::uniform_int_distribution<vid_t> dis;

    int check_edge(const edge_t *e) {
        rep (i, bucket) {
            auto &is_boundary = is_boundarys[i];
            if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
                occupied[i] < capacity) {
                return i;
            }
        }

        rep (i, bucket) {
            auto &is_core = is_cores[i], &is_boundary = is_boundarys[i];
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
        auto &is_boundary = is_boundarys[bucket];

        // 如果已经被加入到边界集，直接返回
        if (is_boundary.get(vid))
            return;
        // 把顶点加入到边界集
        is_boundary.set_bit_unsync(vid);

        // 如果顶点没有在核心集中，直接把顶点的度数据加入到最小堆
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

    bool get_target_vertex(vid_t &vid) {
        // TODO 将随机选择顶点改成选择度最小的顶点，或者是距离当前分区所有节点距离最近的顶点
        // TODO 以上这个计算不太现实
        // 选择度最小的顶点，因为这样跨分区的边从一定概率来说是最小的
        if (d.size() == 0) return false;
        vid_t degree;
        d.get_min(degree, vid);
        d.remove(vid);
        return true;
    }

    void assign_remaining();

    void assign_master();

    size_t count_mirrors();

public:
    NePartitioner(std::string input, std::string algorithm, int num_partition);

    void split() override;


    void calculate_replication_factor() {
        // 每个边集的顶点数求和除以总的顶点数
        for (auto &is_mirror: is_mirrors)
            repv(j, p) {
                if (is_mirror.get(j)) {
                    replicas++;
                }
            }
        // TODO 这个没法计算replication_factor_1
        replicas_1 = replicas - is_mirrors.back().popcount();
        avg_vertex_1 = replicas_1 / (p - 1);

        replication_factor = (double) replicas / num_vertices;
        avg_vertex = replicas / p;
    }

    void calculate_alpha() {
        max_edge = *max_element(occupied.begin(), occupied.end()); // 获取最大值
        min_edge = *min_element(occupied.begin(), occupied.end());

        alpha = (double) max_edge * p / (double) num_edges;
    }

    // TODO 计算方式不对
    void calculate_rho() {
        int variance = 0;
        int variance_1 = 0;
        // 每个分区减去平均
        for (int i = 0; i < num_vertices_in_partition.size() - 1; i++) {
            variance_1 += (num_vertices_in_partition[i] - avg_vertex_1) * (num_vertices_in_partition[i] - avg_vertex_1);
            variance += (num_vertices_in_partition[i] - avg_vertex) * (num_vertices_in_partition[i] - avg_vertex);
        }
        rho_1 = variance_1 / (p - 1);
        rho = variance + (num_vertices_in_partition[num_vertices_in_partition.size() - 1] - avg_vertex) *
                         (num_vertices_in_partition[num_vertices_in_partition.size() - 1] - avg_vertex);
        rho /= p;// 最后一个分区的方差 =]
    }
    // 广度遍历，重新索引，用于将顶点分块
    void re_index() {
        // 随机选择顶点，进行广度遍历，重新索引
        vid_t index = 0;
        vid_t vid = dis(gen);
        // 基于该顶点进行深度遍历，对每个顶点重新索引
        v_queue.push(vid);
        while (!v_queue.empty()) {
            vid_t v = v_queue.front();
            v_queue.pop();
            // 将v加入到indices,重新索引
            indices.insert(std::pair<vid_t, vid_t>(index++, v));

            // 获取v的邻居顶点
            set < vid_t > neighbor_set = adj_list.find(v)->second;
            // 将顶点v的邻居加入到队列中，注意去重
            std::set<int> differenceSet;

            // 使用 std::set_difference 求差值
            std::set_difference(neighbor_set.begin(), neighbor_set.end(),
                                v_set.begin(), v_set.end(),
                                std::inserter(differenceSet, differenceSet.begin()));
            // 将neighbor_set加入v_queue和v_set中
            for (auto &i: differenceSet) {
                v_queue.push(i);
                v_set.insert(i);
            }

        }
    }

    void construct_adj_list(vector<edge_t> &edges) {
        // 遍历边集，建立每个顶点的邻居集合
        for (auto &i: edges) {
            if (adj_list.contains(i.first)) {
                adj_list.find(i.first)->second.insert(i.second);
            } else {
                std::set < vid_t > set;
                set.insert(i.second);
                adj_list[i.first] = set;
            }
            if (adj_list.contains(i.second)) {
                adj_list.find(i.second)->second.insert(i.first);
            } else {
                std::set < vid_t > set;
                set.insert(i.first);
                adj_list[i.second] = set;
            }
        }
    }
};

