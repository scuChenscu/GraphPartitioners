#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
#include <thread>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include "../utils/util.hpp"
#include "../utils/graph.hpp"
#include "../utils//dense_bitset.hpp"
using namespace std;

// int partitions[] = {2, 4, 8, 16, 32, 64};
static const size_t partitions[] = { 32 };
// int partitions[] = {64};
static const int memory_size = 4096;
static const double lambda = 1.1;
static const double balance_ratio = 1.05;
// const string algorithms[] = {"ne", "dbh", "hdrf", "ldg", "fennel"};
static const string algorithms[] = {      "ne" };
// com-amazon.graph不是强连通图，废弃
static  const string graph_suffix = "com-dblp.graph";
static const bool isShuffle = false;
const static string input = "../graphs/medium-scale";
const static bool REINDEX = false;
// Ours参数
static const bool SELF = true;
static const double OURS_BALANCE_RATIO = 1.00;
static const double OURS_CAPACITY_RATIO = 1.05;
static const size_t CORES = thread::hardware_concurrency();
//static const size_t CORES = 1;
//static const size_t MAX_CORES = thread::hardware_concurrency();
//static const double OURS_CAPACITY_RATIOS[] = { 0.15, 0.30, 0.45, 0.60, 0.75,0.90};
//static const double OURS_BALANCE_RATIOS[] = {1.00, 1.05, 1.10};
// 简化
static const size_t MAX_CORES = thread::hardware_concurrency();
static const double OURS_CAPACITY_RATIOS[] = {  0.9};
static const double OURS_BALANCE_RATIOS[] = {1.00};
static  const unsigned DEFAULT_SEED = 985;

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
    vector<unsigned long> degrees;
    size_t max_degree = 0;
    size_t min_degree = UINT64_MAX;
    double avg_degree;
    long long total_degree = 0;
    dense_bitset true_vids;


    string graph_name;

    vector<size_t> reverse_indices;
    vector<size_t> indices;
    dense_bitset visited;
    mt19937 gen;
    //均匀分布区间
    uniform_int_distribution<vid_t> dis;
    // 构造函数
    explicit BaseGraph(const string &graph_name);

    void construct_adjacency_list();
    void partition();
    void re_index();
};