#pragma once
#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include "../../partitioner/partitioner.hpp"
#include "../../utils/dense_bitset.hpp"
#include "../../utils/util.hpp"
#include "../../partitioner/edgePartitioner.hpp"

using namespace std;
//贪心策略 Greedy[11],定义了针对不同类型的边的划分
//规则:1) 两个端点都是新顶点,将边分配到负载最小的分区;
// 2) 一个端点是新顶点,另一个已经在分区中,将边分配到该分区中;
// 3) 两个端点都在分区中,且存在于同一个分区中,将边分配到共同分区中;
// 4) 两个端点都在分区中,但不在同一个分区中,将边分配到这些分区中负载最小的分区
class GreedyPartitioner : public EdgePartitioner {

private:
    // 集合
    vector<set<size_t>> vertex_partitions;
    size_t leastLoad(set<size_t> set);
    void assign_edge(size_t index, vid_t from, vid_t to);

public:
    GreedyPartitioner(BaseGraph& baseGraph,
                      const string& input,
                      const string &algorithm,
                      size_t num_partitions);
    void split() override;
};