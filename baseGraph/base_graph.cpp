//
// Created by 陈键淞 on 2024/1/8.
//

#include "base_graph.hpp"
#include "../algorithms/ne/ne.hpp"
#include "../algorithms/dbh/dbh.hpp"
#include "../algorithms/hdrf/hdrf.hpp"
#include "../algorithms/adwise/adwise.hpp"
#include "../algorithms/ldg/ldg.hpp"
#include "../algorithms/fennel/fennel.hpp"
#include "../algorithms/greedy/greedy.hpp"
#include "../algorithms/rand/rand.hpp"
#include "../algorithms/model4/model4.hpp"
#include "../algorithms/model5/model5.hpp"
#include "../algorithms/model6/model6.hpp"
#include "../algorithms/model7/model7.hpp"
#include "../algorithms/model8/model8.hpp"
#include "../algorithms/model9/model9.hpp"
#include "../algorithms/model10/model10.hpp"
#include "../algorithms/model11/model11.hpp"
#include "../algorithms/model12/model12.hpp"
#include "../algorithms/dne/dne.hpp"
#include "../algorithms/offstreamNG/offstreamNG.hpp"
#include "../algorithms/offstreamNH/offstreamNH.hpp"
#include "../algorithms/offstreamNWG/offstreamNWG.hpp"
#include "../algorithms/timene/timene.hpp"
#include "../algorithms/dcne/dcne.hpp"
#include "../algorithms/bne/bne.hpp"
#include "../algorithms/dpqne/dpqne.hpp"
using namespace std;

BaseGraph::BaseGraph(const string& graph_name) {
    this->graph_name = graph_name;
    ifstream fin;
    if (need_to_shuffle) {
        fin.open(shuffled_binary_edgelist_name(graph_name),
                     ios::binary | ios::ate);
    } else {
        fin.open(binary_edgelist_name(graph_name),
                     ios::binary | ios::ate);
    }

    // tellp 用于返回写入位置，
    // tellg 则用于返回读取位置也代表着输入流的大小
    auto filesize = fin.tellg();
    fin.seekg(0, ios::beg);

    //最大的下标+1
    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));
    LOG(INFO) << "File size: " << filesize
              << " | num_vertices: " << num_vertices
              << " | num_edges: " << num_edges
              << endl;

    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
             filesize);

    edges.resize(num_edges);
    degrees.resize(num_vertices, 0);
    true_vids.resize(num_vertices);

    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    if (algrithm_type == "offstream") {
        // 将edges分为两个集合
        size_t size = edges.size();
        size_t halfSize = size * EDGE_RATIO; // 整数除法得到一半大小
        if (size % 2 != 0) {
            ++halfSize; // 如果总数是奇数，则第二个向量多一个元素
        }

        off_part.reserve(halfSize);
        stream_part.reserve(size - halfSize);

        std::copy(edges.begin(), edges.begin() + halfSize, std::back_inserter(off_part));
        std::copy(edges.begin() + halfSize, edges.end(), std::back_inserter(stream_part));

        CHECK_EQ(edges.size(), off_part.size() + stream_part.size()) << "edges no equals!";

        adj_out.build(off_part);
        adj_in.build_reverse(off_part);

    } else {
        // 初始化的时候构造图
        adj_out.build(edges);
        // 存储反向边
        adj_in.build_reverse(edges);
    }
    gen.seed(DEFAULT_SEED);

    dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    visited = dense_bitset(num_vertices);
    indices.resize(num_vertices);
    reverse_indices.resize(num_vertices);

    construct_adjacency_list();

    // 重新索引
    if (REINDEX) {
        LOG(INFO) << "re_index" << endl;
        re_index();
    } else {
        for (int i = 0; i < num_vertices;i++) {
            indices[i] = i;
        }
    }


}
// TODO 建立CSR
void BaseGraph::construct_adjacency_list() {
    LOG(INFO) << "Construct adjacency list..." << endl;
    // 遍历边集，建立每个顶点的邻居集合
    for (auto &edge: edges) {
        // 计算顶点度数
        degrees[edge.first]++;
        degrees[edge.second]++;
        true_vids.set_bit_unsync(edge.first);
        true_vids.set_bit_unsync(edge.second);
        continue;
        if (adjacency_list.count(edge.first) > 0 ) {
            // LOG(INFO) << edge.first;
            adjacency_list.find(edge.first)->second.insert(edge.second);
        } else {
            set <vid_t> set;
            set.insert(edge.second);
            adjacency_list[edge.first] = set;
        }
        if (adjacency_list.count(edge.second)) {
            // LOG(INFO) << edge.second;
            adjacency_list.find(edge.second)->second.insert(edge.first);
        } else {
            set <vid_t> set;
            set.insert(edge.first);
            adjacency_list[edge.second] = set;
        }
    }
    // 计算顶点度数
    for (int i = 0; i < num_vertices; i++) {
        if (degrees[i] > max_degree) {
            max_degree = degrees[i];
        } else if (degrees[i] < min_degree) {
            min_degree = degrees[i];
        }
        total_degree += degrees[i];
    }
    LOG(INFO) << "total degree: " << total_degree << " , num_vertices: " << num_vertices << endl;
    avg_degree = (double)total_degree / (double )num_vertices;

    // 计算大于平均度数的有多少

    int over_ad = 0;
    for (int i = 0; i < num_vertices; i++) {
        if ((double)degrees[i] > avg_degree) {
            over_ad++;
        }
    }

    LOG(INFO) << "avg_degree: " <<  avg_degree << " over_ad: " << over_ad  << " ratio: " << (double)over_ad / num_vertices << endl;

}

void BaseGraph::partition() {
    // 所有的partitioner共享graph，即在一张图上面执行所有的partition算法，减少重复操作
    LOG(INFO) << "Start partition on " << graph_name << endl;
    stringstream ss;
    ss << "Graph Name: " << graph_name << endl;
    appendToFile(ss.str());
    ss.str("");  // 清空当前字符串内容
    ss.clear();   // 清空错误状态标志
    vector<Partitioner*> partitioners;
    for (auto num_partitions : partitions) {
        Partitioner* partitioner;
//        ss << "Number Partitions: " << num_partitions << endl;
//        appendToFile(ss.str());
//        ss.str("");  // 清空当前字符串内容
//        ss.clear();   // 清空错误状态标志
        for (auto& algorithm : algorithms) {
            // TODO this是一个指针
            if (algorithm == "ne") {
                partitioner = new NePartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "model4") {
                partitioner = new Model4Partitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "model5") {
                // TODO 三个核心参数：cores、balance_ratio、capacity_ratio
                if (SELF) {
                    partitioner = new Model5Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,1);
                    partitioners.push_back(partitioner);
                    continue;
                }
                // 三层循环cores、balance_ratio、capacity_ratio
                for(size_t cores = 1; cores <= MAX_CORES; cores++) {
                    for (auto ours_balance_ratio : OURS_BALANCE_RATIOS) {
                        for (auto ours_capacity_ratio : OURS_CAPACITY_RATIOS) {
                            partitioner = new Model5Partitioner(*this, graph_name, algorithm, num_partitions, ours_balance_ratio, ours_capacity_ratio, cores);
                            partitioners.push_back(partitioner);
                        }
                    }
                }
            } else if (algorithm == "model6") {
                partitioner = new Model6Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "model7") {
                partitioner = new Model7Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "model8") {
                partitioner = new Model8Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "model9") {
                partitioner = new Model9Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "model10") {
                partitioner = new Model10Partitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "model11") {
                partitioner = new Model11Partitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "model12") {
                partitioner = new Model12Partitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "dne") {
                partitioner = new DnePartitioner(*this, graph_name, algorithm, num_partitions,OURS_BALANCE_RATIO, OURS_CAPACITY_RATIO,CORES);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "dbh") {
                partitioner = new DbhPartitioner(*this, graph_name, algorithm, num_partitions, memory_size);
                partitioners.push_back(partitioner);
            } else if (algorithm == "hdrf") {
                partitioner = new HdrfPartitioner(*this, graph_name, algorithm, num_partitions, memory_size,
                                                  balance_ratio, lambda, isShuffle);
                partitioners.push_back(partitioner);
            } else if (algorithm == "adwise") {
                partitioner = new AdwisePartitioner(*this, graph_name, algorithm, num_partitions, memory_size,
                                                  balance_ratio, lambda, isShuffle);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "ldg") {
                partitioner = new LdgPartitioner(*this, graph_name, algorithm, num_partitions, memory_size, isShuffle);
                partitioners.push_back(partitioner);
            } else if (algorithm == "fennel") {
                partitioner = new FennelPartitioner(*this, graph_name, algorithm, num_partitions, memory_size,
                                                    isShuffle);
                partitioners.push_back(partitioner);
            } else if (algorithm == "greedy") {
                partitioner = new GreedyPartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "rand") {
                partitioner = new RandPartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "offstreamNG") {
                partitioner = new OffstreamNGPartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            } else if (algorithm == "offstreamNH") {
            partitioner = new OffstreamNHPartitioner(*this, graph_name, algorithm, num_partitions);
            partitioners.push_back(partitioner);
        }   else if (algorithm == "offstreamNWG") {
            partitioner = new OffstreamNWGPartitioner(*this, graph_name, algorithm, num_partitions);
            partitioners.push_back(partitioner);
        } else if (algorithm == "timene") {
                partitioner = new TimernePartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
        } else if (algorithm == "dcne") {
            partitioner = new DcnePartitioner(*this, graph_name, algorithm, num_partitions);
            partitioners.push_back(partitioner);
        } else if (algorithm == "bne") {
                partitioner = new BnePartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            }
            else if (algorithm == "dpqne") {
                partitioner = new DpqnePartitioner(*this, graph_name, algorithm, num_partitions);
                partitioners.push_back(partitioner);
            }

            else {
                LOG(ERROR) << "Unknown algorithm: " << algorithm;
                continue;
            }
        }
    }
    for (auto partitioner : partitioners) {
        // LOG(INFO) << "Execute " << algorithm << " on: " << graph_name << endl;
        partitioner->split();
        // LOG(INFO) << "Finish " << algorithm << " on: " << graph_name << endl;
        partitioner->calculate_indices();
        delete partitioner;
        // 将指针设置为 nullptr，以防止悬垂指针
    }
}

void BaseGraph::re_index() {
    // LOG(INFO) << adjacency_list.size();
    queue<vid_t> v_queue;
    auto start = std::chrono::high_resolution_clock::now(); // 记录开始时间
    // 随机选择顶点，进行广度遍历，重新索引
    vid_t index = 0;
    vid_t vid = dis(gen);
    // 基于该顶点进行深度遍历，对每个顶点重新索引
    // TODO 该顶点可能没有邻居
    v_queue.push(vid);
    while (!v_queue.empty()) {
        // LOG(INFO) << index;
        vid_t v = v_queue.front();
        v_queue.pop();
        if (visited.get(v)) {
            continue;
        }
        visited.set_bit_unsync(v);
        // 将v加入到indices,重新索引
        reverse_indices[v] = index; // vid所在indices的下标为index
        indices[index++] = v;

        // 获取v的邻居顶点
        if (!adjacency_list.count(v)) continue;
        // LOG(INFO) << v;
        set <vid_t> neighbor_set = adjacency_list.find(v)->second;
        // 将neighbor_set加入v_queue和v_set中
        for (auto &i: neighbor_set) {
            v_queue.push(i);
        }
    }
    auto end = std::chrono::high_resolution_clock::now(); // 记录结束时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // 计算时间差
    LOG(INFO) << "re_index time: " << duration.count() << "ms" << endl;
}