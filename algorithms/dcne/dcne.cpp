#include "dcne.hpp"

#include <utility>

using namespace std;

//固定随机数
// 构造函数
DcnePartitioner::DcnePartitioner(BaseGraph& baseGraph, string input, const string &algorithm,
                             size_t num_partitions)
        : EdgePartitioner(baseGraph, algorithm, num_partitions), input(std::move(input)), gen(985) {
    config_output_files();
    current_partition = 0;
    assigned_edges = 0;
    capacity = num_edges * BALANCE_RATIO / num_partitions + 1;
    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundaries.assign(num_partitions, dense_bitset(num_vertices));
    dis.param(
            uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    front_factor = 0.5;
    avg_factor = 2.0;
    front_partition = num_partitions * front_factor;

}

//最后一个子图就是剩下边组合而成
void DcnePartitioner::assign_remaining() {
    auto &is_boundary = is_boundaries[num_partitions - 1], &is_core = is_cores[num_partitions - 1];
    repv(u, num_vertices) for (auto &i: adj_out[u])
            if (edges[i.v].valid()) {
                assign_edge(num_partitions - 1, u, edges[i.v].second);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(edges[i.v].second);
            }
    repv(i, num_vertices) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            //在其他子图中是核心集的不予理睬，不能设置为本子图的核心集
            rep(j, num_partitions - 1) if (is_cores[j].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
        }
    }
}

size_t DcnePartitioner::count_mirrors() {
    size_t result = 0;
    rep(i, num_partitions) result += is_boundaries[i].popcount();
    return result;
}

void DcnePartitioner::split() {
    total_time.start();
    // 初始化最小堆，用于存储S\C的顶点信息
    high_min_heap.reserve(num_vertices);
    low_min_heap.reserve(num_vertices);
    LOG(INFO) << "Start DC-NE partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "DC-NE" << endl
       << "BALANCE RATIO:" << BALANCE_RATIO
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());


    // 前p-1个分区
    for (current_partition = 0; current_partition < front_partition; current_partition++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid, position;
            if (!get_min_value(degree, vid, position)) { // 当S\C为空时，从V\C中随机选择顶点
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                    break;
                }
                degree = adj_out[vid].size() + adj_in[vid].size();
            } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
                if (position == 0) {
                    low_min_heap.remove(vid);
                } else {
                    high_min_heap.remove(vid);
                }
            }
            // 将顶点加入到核心集C
            occupy_vertex(vid, degree);
        }
        low_min_heap.clear();
        high_min_heap.clear();

        rep(direction, 2) repv(vid, num_vertices) {
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        i++;
                    } else {
                        swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
    }

    for (; current_partition < num_partitions - 1; current_partition++) {
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid;
            if (!low_min_heap.get_min(degree,vid)) { // 当S\C为空时，从V\C中随机选择顶点
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                    break;
                }
                degree = adj_out[vid].size() + adj_in[vid].size();
            } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
                low_min_heap.remove(vid);
            }
            // 将顶点加入到核心集C
            occupy_vertex(vid, degree);
        }
        low_min_heap.clear();

        rep(direction, 2) repv(vid, num_vertices) {
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        i++;
                    } else {
                        swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
    }


    current_partition = num_partitions - 1;
    // 把剩余的边放入最后一个分区
    assign_remaining();
    CHECK_EQ(assigned_edges, num_edges);
    total_time.stop();
}

bool DcnePartitioner::get_free_vertex(vid_t &vid) {
    //随机选择一个节点
    vid = dis(gen);
    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices &&
           (adj_out[vid].size() + adj_in[vid].size() == 0 ||
            adj_out[vid].size() + adj_in[vid].size() >
            2 * avg_degree ||
            is_cores[current_partition].get(vid))) {
        vid = (vid + ++count) % num_vertices;
    }
    if (count == num_vertices)
        return false;
    return true;
}

void DcnePartitioner::occupy_vertex(vid_t vid, vid_t degree) {
    CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
    is_cores[current_partition].set_bit_unsync(vid);
    if (degree == 0) return;
    add_boundary(vid);

    for (auto &i: adj_out[vid]) {
        if (edges[i.v].valid())
            add_boundary(edges[i.v].second);
    }
    adj_out[vid].clear();

    for (auto &i: adj_in[vid]) {
        if (edges[i.v].valid())
            add_boundary(edges[i.v].first);
    }
    adj_in[vid].clear();
}

void DcnePartitioner::assign_edge(size_t partition, vid_t from, vid_t to) {
    is_mirrors[from].set_bit_unsync(partition);
    is_mirrors[to].set_bit_unsync(partition);
    assigned_edges++;
    occupied[partition]++;
}

void DcnePartitioner::add_boundary(vid_t vid) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[current_partition];
    auto &is_boundary = is_boundaries[current_partition];

    if (is_boundary.get(vid)) {
        return;
    }
    is_boundary.set_bit_unsync(vid);

    if (!is_core.get(vid)) {
        vid_t degree = adj_out[vid].size() + adj_in[vid].size();
        if (current_partition < front_partition) {
            if (degree < avg_factor * avg_degree) {
                low_min_heap.insert(degree, vid);
            } else {
                high_min_heap.insert(degree, vid);
            }
        } else {
            low_min_heap.insert(degree, vid);
        }

    }

    rep (direction, 2) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        // 遍历顶点的邻边
        for (size_t i = 0; i < neighbors.size();) {
            // 判断邻居edges[neighbors[i].v]是否已经分配
            if (edges[neighbors[i].v].valid()) {
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                if (is_core.get(u)) { // 如果顶点在核心集中，不考虑负载
                    assign_edge(current_partition, direction ? vid : u,
                                direction ? u : vid);
                    if (current_partition < front_partition) {
                        if (high_min_heap.contains(vid)) {
                            high_min_heap.decrease_key(vid);
                        } else {
                            low_min_heap.decrease_key(vid);
                        }
                    } else {
                        low_min_heap.decrease_key(vid);
                    }

                    edges[neighbors[i].v].remove();
                    swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else if (is_boundary.get(u) &&
                           occupied[current_partition] < capacity) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                    // 将边加入边集
                    assign_edge(current_partition, direction ? vid : u, direction ? u : vid);
                    if (current_partition < front_partition) {
                        if (high_min_heap.contains(vid)) {
                            high_min_heap.decrease_key(vid);
                        } else {
                            low_min_heap.decrease_key(vid);
                        }

                        if (high_min_heap.contains(u)) {
                            high_min_heap.decrease_key(u);
                        } else {
                            low_min_heap.decrease_key(u);
                        }
                    } else {
                        low_min_heap.decrease_key(u);
                        low_min_heap.decrease_key(vid);
                    }
                    edges[neighbors[i].v].remove();
                    swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else
                    i++;
            } else {
                swap(neighbors[i], neighbors.back());
                neighbors.pop_back();
            }
        }
    }
}

bool DcnePartitioner::get_min_value(vid_t & degree, vid_t& vid, vid_t& position) {
    if (low_min_heap.size() == 0 && high_min_heap.size() == 0) return false;
    if (high_min_heap.size() == 0) {
        low_min_heap.get_min(degree, vid);
        position = 0;
        return true;
    }
    if (low_min_heap.size() == 0) {
        high_min_heap.get_min(degree, vid);
        position = 1;
        return true;
    }
    vid_t degree2;
    vid_t vid2;
    if (high_min_heap.get_min(degree, vid) && low_min_heap.get_min(degree2, vid2)) {
        if (degree2 < degree) {
            degree = degree2;
            vid = vid2;
            position = 0;
        } else {
            position = 1;
        }
    }
    return true;
}

void DcnePartitioner::remove(vid_t vid, vid_t position) {
    if (position == 0) {
        low_min_heap.remove(vid);
    } else {
        high_min_heap.remove(vid);
    }
}