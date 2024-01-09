#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "partitioner.hpp"
#include "../baseGraph/base_graph.hpp"
#include "../utils/util.hpp"


class VertexPartitioner : public Partitioner {
public:
    explicit VertexPartitioner(BaseGraph &baseGraph,const string& algorithm, size_t num_partitions);
protected:
    // 边割数
    size_t edge_cut;
    // 边割率
    double edge_cut_rate;
    // 负载均衡
    double load_balance;

};