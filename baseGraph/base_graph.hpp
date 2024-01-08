#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
#include "../utils/util.hpp"
#include "../utils/graph.hpp"


using namespace std;

// int partitions[] = {2, 4, 8, 16, 32, 64};
static const size_t partitions[] = { 8, 16, 32, 64};
// int partitions[] = {64};
static const int memory_size = 4096;
static const double lambda = 1.1;
static const double balance_ratio = 1.05;
// const string algorithms[] = {"ne", "dbh", "hdrf", "ldg", "fennel"};
static const string algorithms[] = { "ne"};
static const bool isShuffle = false;

class BaseGraph {
public:
    // 邻接表，每个顶点的邻居顶点
    unordered_map<vid_t, set<vid_t>> adjacency_list;
    // 总边数
    size_t num_edges;
    // 总顶点数
    vid_t num_vertices;
    // 存储边
    vector<edge_t> edges;
    // 存储点
    vector<vid_t> vertices;
    // 存储每个顶点的邻边，区分出边和入边
    graph_t adj_out;
    graph_t adj_in;
    
    string graph_name;

    // 构造函数
    explicit BaseGraph(const string &input);

    // 析构函数
//    ~BaseGraph();

    void construct_adj_list();
    void partition();
};