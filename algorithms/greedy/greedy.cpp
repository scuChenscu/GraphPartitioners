#include <climits>
#include "greedy.hpp"

GreedyPartitioner::GreedyPartitioner(BaseGraph& baseGraph,
                                     const string& input,
                                     const string &algorithm,
                                     size_t num_partitions) : EdgePartitioner(baseGraph, algorithm, num_partitions){
    vertex_partitions.assign(num_vertices, vector<int>(num_partitions));
    vp_set.assign(num_vertices, set<size_t>());
}

//贪心策略 Greedy,定义了针对不同类型的边的划分
//规则:1) 两个端点都是新顶点,将边分配到负载最小的分区;
// 2) 一个端点是新顶点,另一个已经在分区中,将边分配到该分区中;
// 3) 两个端点都在分区中,且存在于同一个分区中,将边分配到共同分区中;
// 4) 两个端点都在分区中,但不在同一个分区中,将边分配到这些分区中负载最小的分区
void GreedyPartitioner::split() {
    stringstream ss;
    ss << "Greedy" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    total_time.start();

    for (auto &edge: edges) {
        vid_t first = edge.first;
        vid_t second = edge.second;

        vector<int> first_partition = vertex_partitions[first];
        vector<int> second_partition = vertex_partitions[second];

        set<size_t>& fp = vp_set[first];
        set<size_t>& sp = vp_set[second];


        set<size_t> intersections;
        set<size_t> unions;
        for (int i = 0; i < num_partitions; i++) {
            if (first_partition[i] == 1 && second_partition[i] == 1) {
                intersections.insert(i);
            }
            if (first_partition[i] == 1 || second_partition[i] == 1) {
                unions.insert(i);
            }
        }

        size_t partition = 0;
        if (!intersections.empty()) {
            partition = leastLoad(intersections);
        } else if (intersections.empty() && !unions.empty()) {
            partition = leastLoad(unions);
        } else if (fp.empty() && !sp.empty()) {
            partition = leastLoad(fp);
        } else if (sp.empty() && !fp.empty()) {
            partition = leastLoad(fp);
        } else if (fp.empty() && sp.empty()){
            // 找出occupied最小值所在的下标
            int min = occupied[0];
            for (size_t i = 1; i < num_partitions; i++) {
                if (occupied[i] < min) {
                    min = occupied[i];
                    partition = i;
                }
            }
        }

        assign_edge(partition, first, second);
        vertex_partitions[first][partition] = 1;
        vertex_partitions[second][partition] = 1;
        fp.insert(partition);
        sp.insert(partition);
    }
    total_time.stop();
}
size_t GreedyPartitioner::leastLoad(set<size_t> set) {
    // 遍历集合元素，找出最小occupied负载
    int min = INT_MAX;
    size_t partition_id;
    for (auto &partition: set) {
        if (occupied[partition] < min) {
            min = occupied[partition];
            partition_id = partition;
        }
    }
    return partition_id;
}

void GreedyPartitioner::assign_edge(size_t index, vid_t from, vid_t to) {
    true_vids.set_bit_unsync(from);
    true_vids.set_bit_unsync(to);
    is_mirrors[from].set_bit_unsync(index);
    is_mirrors[to].set_bit_unsync(index);
    assigned_edges++;
    occupied[index]++;
}
