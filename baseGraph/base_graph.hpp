#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
// #include "../partitioner/partitioner.hpp"
#include "../utils/util.hpp"
#include "../utils/graph.hpp"


using namespace std;
class Partitioner; // 解决循环引用的问题
class BaseGraph {
public:
    // 邻接表，每个顶点的邻居顶点
    unordered_map<vid_t, set<vid_t>> adjacency_list;
    // 总边数
    size_t num_edges{};
    // 总顶点数
    size_t num_vertices{};
    // 存储边
    vector<edge_t> edges;
    // 存储点
    vector<vid_t> vertices;
    // 存储每个顶点的邻边，区分出边和入边
    graph_t adj_out;
    graph_t adj_in;

    // 构造函数
    explicit BaseGraph(const string &input) {
        ifstream fin(binary_edgelist_name(input),
                     std::ios::binary | std::ios::ate);
        // tellp 用于返回写入位置，
        // tellg 则用于返回读取位置也代表着输入流的大小
        auto filesize = fin.tellg();
        fin.seekg(0, std::ios::beg);

        //最大的下标+1
        fin.read((char *) &num_vertices, sizeof(num_vertices));
        fin.read((char *) &num_edges, sizeof(num_edges));

        CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
                 filesize);

        edges.resize(num_edges);
        fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
        // 初始化的时候构造图
        adj_out.build(edges);
        // 存储反向边
        adj_in.build_reverse(edges);

        construct_adj_list();
    };

    // 析构函数
    ~BaseGraph() {
        edges.clear();
        vertices.clear();
        // delete adj_out;
        // delete adj_in;
    }

    void construct_adj_list() {
        // 遍历边集，建立每个顶点的邻居集合
        for (auto &edge: edges) {
            if (adjacency_list.contains(edge.first)) {
                adjacency_list.find(edge.first)->second.insert(edge.second);
            } else {
                std::set < vid_t > set;
                set.insert(edge.second);
                adjacency_list[edge.first] = set;
            }
            if (adjacency_list.contains(edge.second)) {
                adjacency_list.find(edge.second)->second.insert(edge.first);
            } else {
                std::set < vid_t > set;
                set.insert(edge.first);
                adjacency_list[edge.second] = set;
            }
        }
    }
    // TODO 因为每次执行图分区算法会依赖很多的数据结构，只能通过重新创建partitioner的方式来实现
    static void partition() {
        // 在这里去创建partitioner
        Partitioner *partitioner;

    }
};