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

#include "../../utils/dense_bitset.hpp"
#include "../../partitioner/partitioner.hpp"
#include "../ne/graph.hpp"
#include "../ne/min_heap.hpp"
#include "../../partitioner/partitioner.hpp"
#include "../../utils/util.hpp"

/* Neighbor Expansion (NE) */
class OursPartitioner : public Partitioner {
private:
    // TODO 因为NE算法在设定负载均衡的时候，是直接给个参数，会导致最后一个分区可能出现顶点数较少的情况，导致负载均衡做的不够好
    // 如果我们能够在图遍历的过程中，动态地计算到这个负载均衡的值，可能实现边割率和负载均衡的同时优化
    // 比如根据当前已分配顶点跟外界的边数，当前分区负载距离标准负载的距离，来动态地调整负载均衡
    const double BALANCE_RATIO = 1.00;

    std::string input;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    int p, bucket;
    double average_degree;
    size_t capacity;

    vector<vector<vid_t> > part_degrees;
    vector<int> balance_vertex_distribute;
    // 这玩意存储啥？像是存储所有的边，edge_t 是一个first second的结构图，所有edges就是存储边的
    std::vector<edge_t> edges;
    // 图结构
    graph_t adj_out, adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    //每个分区边的数量
    std::vector<size_t> occupied;
    // 顶点的度
    std::vector<vid_t> degrees;
    MinHeap<vid_t, vid_t> d;
    std::vector<int8_t> master;
    std::vector<dense_bitset> is_cores, is_boundarys;
    dense_bitset true_vids;
    vector<dense_bitset> is_mirrors;

    vector<int> connected_components;

    int components;

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
//        d.decrease_key(from);
//        d.decrease_key(to);
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
    OursPartitioner(std::string input, std::string algorithm, int num_partition);

    void split() override;

    virtual void calculate_replication_factor();
};

