#pragma once
#include "vertexPartitioner.hpp"

VertexPartitioner::VertexPartitioner(BaseGraph &baseGraph, const string& algorithm, const size_t num_partitions) : Partitioner(baseGraph, algorithm, num_partitions) {

}