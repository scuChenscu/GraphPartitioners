#include "rand.hpp"

RandPartitioner::RandPartitioner(BaseGraph& baseGraph, const string &input, const string& algorithm, size_t num_partitions) : EdgePartitioner(baseGraph, algorithm, num_partitions) {};
void RandPartitioner::split() {
    stringstream ss;
    ss << "rand" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    total_time.start();

    for (int i = 0; i < edges.size(); i++) {
        auto &edge = edges[i];

        size_t partition = i % num_partitions;
        assign_edge(partition, edge.first, edge.second);
    }
    total_time.stop();
}

void RandPartitioner::assign_edge(size_t index, vid_t from, vid_t to) {
    true_vids.set_bit_unsync(from);
    true_vids.set_bit_unsync(to);
    is_mirrors[from].set_bit_unsync(index);
    is_mirrors[to].set_bit_unsync(index);
    assigned_edges++;
    occupied[index]++;
}