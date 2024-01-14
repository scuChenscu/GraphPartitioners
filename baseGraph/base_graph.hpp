#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
#include <thread>
#include "../utils/util.hpp"
#include "../utils/graph.hpp"
#include "../utils//dense_bitset.hpp"
using namespace std;

// int partitions[] = {2, 4, 8, 16, 32, 64};
static const size_t partitions[] = { 4, 8};
// int partitions[] = {64};
static const int memory_size = 4096;
static const double lambda = 1.1;
static const double balance_ratio = 1.05;
// const string algorithms[] = {"ne", "dbh", "hdrf", "ldg", "fennel"};
static const string algorithms[] = {  "model5"};
static  const string graph_suffix = "4elt.graph";
static const bool isShuffle = false;
const static string input = "../graphs/input";


// Ours参数
static const bool SELF = true;
static const double OURS_BALANCE_RATIO = 1.00;
static const double OURS_CAPACITY_RATIO = 0.00;
//static const size_t CORES = 1;
//static const size_t MAX_CORES = thread::hardware_concurrency();
//static const double OURS_CAPACITY_RATIOS[] = { 0.15, 0.30, 0.45, 0.60, 0.75,0.90};
//static const double OURS_BALANCE_RATIOS[] = {1.00, 1.05, 1.10};
// 简化
static const size_t CORES = 1;
static const size_t MAX_CORES = 4;
static const double OURS_CAPACITY_RATIOS[] = { 0.30, 0.45, 0.60};
static const double OURS_BALANCE_RATIOS[] = {1.00};


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
    // 存储每个顶点的度数
    vector<size_t> degrees;
    dense_bitset true_vids;


    string graph_name;

    // 构造函数
    explicit BaseGraph(const string &graph_name);

    void construct_adjacency_list();
    void partition();
};