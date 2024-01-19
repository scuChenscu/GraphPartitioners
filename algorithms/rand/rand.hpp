#include <random>
#include "../../partitioner/edgePartitioner.hpp"

class RandPartitioner : public EdgePartitioner {
public:
    RandPartitioner(BaseGraph& baseGraph, const string &input, const string& algorithm, size_t num_partitions);
    void split() override;
private:

    void assign_edge(size_t index, vid_t from, vid_t to);
};
