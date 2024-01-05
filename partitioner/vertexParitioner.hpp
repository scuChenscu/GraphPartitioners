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
class VertexPartitioner : public Partitioner {
protected:
    // 边割率
    size_t edge_cut;
    double edge_cut_rate;
    // 负载均衡
    double load_balance;
};