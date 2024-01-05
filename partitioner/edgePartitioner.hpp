#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "partitioner.hpp"
#include "../utils/util.hpp"

using namespace std;

/**
 * 以边为基本划分单位的图分区算法
 */
class EdgePartitioner : public Partitioner {

protected:
    // 复制因子
    double replication_factor;
    // 负载均衡因子
    double alpha;
    // 顶点总副本数
    size_t replicas;
    // 分区数
    size_t partition;

    size_t max_edge = numeric_limits<size_t>::min();

    size_t min_edge = numeric_limits<size_t>::max();

    // 分区边负载上限调节因子
    double balance_ratio;

    // TODO 邻接表，每个顶点的邻居顶点；该数据应该来自Graph
    unordered_map<vid_t, set<vid_t>> adjacency_list;
    // 每个分区的边的数量
    vector<size_t> occupied;
    // 每个顶点所属的分区，每个顶点可能属于多个分区
    vector<dense_bitset> is_mirrors;

    // 计算副本
    void calculate_replication_factor()  {
        for(auto &is_mirror : is_mirrors) {
            replicas += is_mirror.popcount();
        }
        replication_factor = (double) replicas / (double)graph.num_vertices;
    }

    // 计算负载均衡因子
    void calculate_alpha()  {
        max_edge = *max_element(occupied.begin(), occupied.end()); // 获取最大值
        min_edge = *min_element(occupied.begin(), occupied.end());

        alpha = (double) max_edge * (double)partition / (double) graph.num_edges;
    }
};