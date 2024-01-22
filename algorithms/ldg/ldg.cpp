#include "ldg.hpp"

// LDG是把点分到不同的分区，需要计算边割率
// LDG计算：遍历所有的边，判断边的两个顶点是否在同一个分区
// Linear Deterministic Greedy(LDG)
// LDG是一种贪心算法，它以顶点作为输入流，是一种点分区算法。
// 它希望能把顶点分配到邻居最多的分区，以减小跨分区边的数量。
// LDG需要保存之前的顶点信息，因此不适用于无边界流。
// LDG和Fennel算法的区别在于partition_score的计算方式不太一样
LdgPartitioner::LdgPartitioner(BaseGraph &baseGraph, const string &input, const string &algorithm,
                               const size_t num_partitions, int memory_size, bool shuffle) :
        VertexPartitioner(baseGraph, algorithm, num_partitions) {
    config_output_files();

    //edge file
    if (shuffle) {
        fin.open(shuffled_binary_edgelist_name(input), std::ios::binary | std::ios::ate);
    } else {
        fin.open(binary_edgelist_name(input), std::ios::binary | std::ios::ate);
    }
    filesize = fin.tellg();
    fin.seekg(0, std::ios::beg);

    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));

    num_batches = (filesize / ((std::size_t) memory_size * 1024 * 1024)) + 1;
    num_edges_per_batch = (num_edges / num_batches) + 1;

    subg_vids.assign(num_partitions, unordered_set < vid_t > {});
//    true_vids.resize(num_vertices);
    balance_vertex_distribute.resize(num_vertices);
    node2neis.assign(num_vertices, unordered_set < vid_t > {});
    LOG(INFO) << "finish init";
}


void LdgPartitioner::read_and_do(string opt_name) {
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
    std::vector<edge_t> edges;
    auto num_edges_left = num_edges;
    for (uint32_t i = 0; i < num_batches; i++) {
        auto edges_per_batch = num_edges_per_batch < num_edges_left ? num_edges_per_batch : num_edges_left;
        edges.resize(edges_per_batch);
        fin.read((char *) &edges[0], sizeof(edge_t) * edges_per_batch);
        if (opt_name == "node_assignment") {
            batch_node_assignment(edges);
            LOG(INFO) << "finish edge";
        } else if (opt_name == "process neighbors") {
            for (auto &e: edges) {
                addNeighbors(e);
                true_vids.insert(e.first);
                true_vids.insert(e.second);
            }
            LOG(INFO) << "finish neis";
        } else {
            LOG(ERROR) << "no valid opt function";
        }
        num_edges_left -= edges_per_batch;
    }
}

void LdgPartitioner::addNeighbors(edge_t &edge) {
    node2neis[edge.first].insert(edge.second);
    node2neis[edge.second].insert(edge.first);
}

size_t LdgPartitioner::intersection(vid_t vid, size_t partition) {
    if (!adjacency_list.count(vid)) return 0;
    int count = 0;
    for (auto neighbor : adjacency_list.find(vid)->second) {
        if (partition_vertices[partition].get(neighbor)) {
            count++;
        }
    }
    // LOG(INFO) << count;
    return count;
}

void LdgPartitioner::do_ldg() {
    // shuffle_vertices of streaming vertices
    vector<int> shuffle_vertices(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        shuffle_vertices[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(shuffle_vertices.begin(), shuffle_vertices.end(), g);

    // Initial partitions
    size_t i = 0;
    for (; i < num_partitions; i++) {
        assigned_vertex(shuffle_vertices[i], i);
    }

    for (; i < num_vertices; i++) {
        if (i % 10000 == 0)
            LOG(INFO) << i << "/" << num_vertices << endl;
        vid_t v = shuffle_vertices[i];
        vector<double> scores(num_partitions, 0);
        for (int id = 0; id < num_partitions; id++) {
            double weightedGreedy = 1 - (double)num_vertices_each_partition[id] / (double)capacity;
            // TODO 如果每个顶点的邻居都很少，可能会导致intersection为0，出现大量0的情况
            double neighbors_in_partition = intersection(v, id) + 1e-6;
            scores[id] = neighbors_in_partition * weightedGreedy;
        }
        //最大值所在序列的位置
        int firstIndex = distance(scores.begin(),
                                  max_element(scores.begin(), scores.end()));
        assigned_vertex(v, firstIndex);
    }
    LOG(INFO) << "finish ldg";
}

void LdgPartitioner::batch_node_assignment(vector<edge_t> &edges) {
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
        if (sp != tp) {
            edge_cut++;
        }
    }

    edge_cut_rate = double(edge_cut) / edges.size();
}

void LdgPartitioner::split() {
    // read_and_do("process neighbors");

    stringstream ss;
    ss << "LDG" << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());
    total_time.start();

    do_ldg();
    // read_and_do("node_assignment");
    total_time.stop();
}