#include "model9.hpp"

using namespace std;

Model9Partitioner::Model9Partitioner(BaseGraph& baseGraph,
                                     const string &input,
                                     const string &algorithm,
                                     size_t num_partitions,
                                     double ours_balance_ratio,
                                     double ours_capacity_ratio,
                                     size_t cores) : EdgePartitioner(baseGraph, algorithm, num_partitions),
                                                              input(input), gen(989) {
    config_output_files();

    this->cores = cores;
    this->balance_ratio = ours_balance_ratio;
    this->capacity_ratio = ours_capacity_ratio;


    current_partition = 0;
    average_degree = (double) num_edges * 2 / num_vertices;
    capacity = (double) num_edges * balance_ratio / (double)num_partitions + 1;
    num_vertices_each_cores = num_vertices / cores;

//    adj_directed.resize(num_vertices);
    assigned = dense_bitset(num_edges);
    edge_partition.resize(num_edges);
    vertex_lock = dense_bitset(num_vertices);

    v_lock.resize(num_vertices, 0);
    // 每个都分配num_vertices大小，保证不会有问题
    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundaries.assign(num_partitions, dense_bitset(num_vertices));
    true_vids.resize(num_vertices);
    is_mirrors.assign(num_vertices, dense_bitset(num_partitions));

    reverse_is_mirrors.assign(num_partitions, dense_bitset(num_vertices));

    master.assign(num_vertices, -1);
    dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    sub_dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices_each_cores - 1));

    edges.resize(num_edges);

    vertex_adjacent_edges.resize(num_vertices, unordered_map<size_t, edge_t *>());

    // degrees.resize(num_vertices);
//    std::ifstream degree_file(degree_name(input), std::ios::binary);
//    degree_file.read((char *) &degrees[0], num_vertices * sizeof(vid_t));
//    degree_file.close();

    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);
    // acquired_partition.resize(num_partitions, -1);


}

//最后一个子图就是剩下边组合而成
void Model9Partitioner::assign_remaining() {
    auto &is_boundary = is_boundaries[num_partitions - 1];
    auto &is_core = is_cores[num_partitions - 1];

    int rep = 0;
    // 遍历边集合
    for (auto &e: edges) {
        if (e.valid()) {
            assign_edge(num_partitions - 1, e.first, e.second);
            if (!is_boundary.get(e.first)) {
                rep++;
            }
            if (!is_boundary.get(e.second)) {
                rep++;
            }

            // TODO 这里要判断核心集合吗？一个顶点只能在一个核心集
            is_boundary.set_bit_unsync(e.first);
            is_boundary.set_bit_unsync(e.second);
            is_core.set_bit_unsync(e.first);
            is_core.set_bit_unsync(e.second);
            rep(partition, num_partitions - 1) {
                if (is_cores[partition].get(e.first)) {
                    is_core.set_unsync(e.first, false);
                }
                if (is_cores[partition].get(e.second)) {
                    is_core.set_unsync(e.second, false);
                }
            }
        }
    }
    LOG(INFO) << "Rep: " << rep;
}

void Model9Partitioner::assign_master() {
    vector<vid_t> count_master(num_partitions, 0);
    vector<vid_t> quota(num_partitions, num_vertices);
    long long sum = num_partitions * num_vertices;
    uniform_real_distribution<double> distribution(0.0, 1.0);
    vector<dense_bitset::iterator> pos(num_partitions);
    rep(b, num_partitions) pos[b] = is_boundaries[b].begin();
    vid_t count = 0;
    while (count < num_vertices) {
        long long r = distribution(gen) * sum;
        //随机选择哪个子图先赋值
        int k;
        for (k = 0; k < num_partitions; k++) {
            if (r < quota[k])
                break;
            r -= quota[k];
        }
        //选出当前位置还未被赋值的结点
        while (pos[k] != is_boundaries[k].end() && master[*pos[k]] != -1)
            pos[k]++;
        if (pos[k] != is_boundaries[k].end()) {
            count++;
            master[*pos[k]] = k;
            count_master[k]++;
            quota[k]--;
            sum--;
        }
    }
}

void Model9Partitioner::split() {
    LOG(INFO) << "Build vertex adjacent edges" << endl;
    // 记录每个顶点的边数
    build_vertex_adjacent_edges();

    // 根据边的两端顶点所属的cores建立有向边
//    LOG(INFO) << "Build direction edge" << endl;
//    adj_directed.build_directed(edges, reverse_indices);

    min_heaps.resize(num_partitions, MinHeap<vid_t, vid_t>(num_vertices));
    capacity_ratios.resize(num_partitions, 0);
    expansions.resize(num_partitions,0);
    LOG(INFO) << "Start Model9 partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "Model9"
       << endl
       << "去除BFS分块，使用K-way多线程来同时进行k个分区划分，最后将剩余的边使用HDRF来分配"
       << "目标是提高NE算法的运行速度，同时保证RF尽可能低"
       << "对不同的线程设置不同的扩展速度和负载上限"
       << endl
       << "Balance Ratio: " << balance_ratio
       << " | Capacity Ratio: " << capacity_ratio
       << " | Cores: " << cores
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    // 启动线程
    if (cores > 1) {
        if (cores > num_partitions) {
            cores = num_partitions;
        }
        // cores = 1;
        // 扩展速度和负载上限
        for (int i = 0; i < cores; i++) {
            expansions[i] = cores - i;
            // expansions[i] = 1;
            // 主要是容量上限带来的影响
            capacity_ratios[i] = capacity_ratio * (1 - i * 0.05);
            // capacity_ratios[i] = 0;
        }

        thread threads[cores];
        for (int i = 0; i < cores; ++i) {
            threads[i] = thread(&Model9Partitioner::sub_split, this, i);
        }
        total_time.start();
        Timer t;
        t.start();
        for (int i = 0; i < cores; ++i) {
            threads[i].join();
        }
        t.stop();
        // LOG(INFO) << "Model9 multi-partitioning time: " << t.get_time() << endl;
        // LOG(INFO) << "Model9 multi-partitioning finished!" << endl;
    } else {
        total_time.start();
        LOG(INFO) << "Cores: " << cores << " , less than 1, use single thread" << endl;
    }

    for (int p = 0; p < num_partitions; p++) {
        if (occupied[p] > max_load) max_load = occupied[p];
        if (occupied[p] < min_load) min_load = occupied[p];
    }

    // LOG(INFO) << "Start building the first num_partitions - 1 full partition" << endl;
    // 使用HDRF算法分配剩余的边
    Timer t;
    t.start();
    for (auto &edge: edges) {
        if (edge.valid()) {
             // TODO 一共有三种情况：1. 该边两个顶点都没有被抢占，2.该边两个顶点都被抢占，3.该边只有一个顶点被抢占
            int partition = find_max_score_partition(edge);
//            is_mirrors[edge.first].set_bit_unsync(partition);
//            is_mirrors[edge.second].set_bit_unsync(partition);
            assign_edge(partition, edge.first, edge.second);
            update_min_max_load(partition);
        }
    }
    // 换一种方式
//    for (auto &edge: edges) {
//        if (edge.valid()) { // 用贪心的方式
//
//        }
//    }
    t.stop();
    // LOG(INFO) << "HDRF time: " << t.get_time() << endl;
//    for (current_partition = 0; current_partition < num_partitions - 1; current_partition++) {
//        Timer t;
//        t.start();
//        // min_heap = min_heaps[current_partition];
//        // LOG(INFO) << "Min_heap " << current_partition << "  size: " << min_heaps[current_partition].size() << endl;
//        while (occupied[current_partition] < capacity * balance_ratio) {
//            vid_t degree, vid;
//            bool has_min = false;
//            while (min_heaps[current_partition].size() > 0) {
//                min_heaps[current_partition].get_min(degree, vid);
//                min_heaps[current_partition].remove(vid);
//                if (!vertex_adjacent_edges[vid].empty()) {
//                    // 更新degree
//                    // CHECK_EQ(degree, vertex_adjacent_edges[vid].size());
//                    degree = vertex_adjacent_edges[vid].size();
//                    has_min = true;
//                    break;
//                }
//            }
//            if (!has_min) { // 当S\C为空时，从V\C中随机选择顶点
//                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
//                    break;
//                }
//                // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
//                degree = vertex_adjacent_edges[vid].size();
//            }
//
//            CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
//            // 把顶点x加入到S和C
//            is_cores[current_partition].set_bit_unsync(vid);
//            is_boundaries[current_partition].set_bit_unsync(vid);
//
//
//            auto iterator = vertex_adjacent_edges[vid].begin();
//            while(iterator != vertex_adjacent_edges[vid].end()){
//                edge_t &edge = *iterator->second;
//                if (edge.valid()) {
//                    vid_t& neighbor = edge.first == vid ? edge.second : edge.first;
//                        if (edge.valid()) {
//                            add_boundary(neighbor);
//                            if (!edge.valid()) {
//                                vertex_adjacent_edges[vid].erase(iterator++);
//                            } else iterator++;
//                        } else {
//                            iterator++;
//                        }
//                } else iterator++;
//            }
//        }
//        t.stop();
//        LOG(INFO) << "Partition " << current_partition << " time: " << t.get_time() << endl;
//        LOG(INFO) << "End partition " << current_partition  << " edge count: " << occupied[current_partition] << endl;
//
//    }
//    LOG(INFO) << "Assigned edges: " << assigned_edges;
//    LOG(INFO) << "Remaining edges: " << num_edges - assigned_edges;
//    // 能否修改为将剩余的边手动加入到最合适的分区？
//    LOG(INFO) << "Start building the last full partition" << endl;
//    current_partition = num_partitions - 1;
//    // 把剩余的边放入最后一个分区
//    Timer t;
//    t.start();
//    assign_remaining();
//    t.stop();
//    LOG(INFO) << "Partition " << current_partition << " time: " << t.get_time() << endl;
//    LOG(INFO) << "End partition " << current_partition  << " edge count: " << occupied[current_partition] << endl;
//
//    CHECK_EQ(assigned_edges, num_edges);
    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time() << endl;
}

// TODO 这里要注意不能用current_partition，current_partition是全局变量，只能用在后续建立完整分区
void Model9Partitioner::sub_split(size_t index) {
    // LOG(INFO) << "Start sub_split " << index << endl;
    Timer t;
    t.start();
    // LOG(INFO) << capacity * capacity_ratio;
    while (occupied[index] <= capacity * capacity_ratios[index]) {
        // 尝试从当前分区所属的最小堆中获取顶点
        vector<vid_t> vids;
            vid_t degree, v;
            if (min_heaps[index].size() > 0) {
                while (min_heaps[index].size() > 0 && vids.size() < expansions[index]) {
                    min_heaps[index].get_min(degree, v);
                    // TODO，这里不能移除太早
                    min_heaps[index].remove(v);
                    vids.push_back(v);
                }
            } else {
                while(vids.size() < expansions[index]) {
                    if (sub_get_free_vertex(v, index)) {
                        degree = vertex_adjacent_edges[v].size();
                        vids.push_back(v);
                    } else {
                        break;
                    }
                }
            }
        if (vids.empty()) break;
        // 一次性选择多个顶点加入到核心集
        for (vid_t vid : vids) {
            CHECK(!is_cores[index].get(vid)) << "add " << vid << " to core again";
            is_cores[index].set_bit_unsync(vid);
            // 可能是冗余操作，但是不重要
            is_boundaries[index].set_bit_unsync(vid);
            // 把他们的邻居全加入到S\C中

        }

        for (vid_t vid : vids) {
//            CHECK(!is_cores[index].get(vid)) << "add " << vid << " to core again";
//            is_cores[index].set_bit_unsync(vid);
//            is_boundaries[index].set_bit_unsync(vid);
            auto iterator = vertex_adjacent_edges[vid].begin();
            while(iterator != vertex_adjacent_edges[vid].end()){
                edge_t &edge = *iterator->second;
                if (edge.valid()) {
                    vid_t neighbor = edge.first == vid ? edge.second : edge.first;
                    if (acquire_vertex(neighbor)) {
                        if (edge.valid()) {
                            sub_add_boundary(neighbor, index);
                            if (!edge.valid()) {
                                vertex_adjacent_edges[vid].erase(iterator++);
                            } else iterator++;
                        } else {
                            release_vertex(neighbor);
                            iterator++;
                        }
                    } else iterator++;
                } else iterator++;
            }
            release_vertex(vid);
        }
        vids.clear();
    }

    t.stop();
    // LOG(INFO) << "End sub_split " << index << " time: " << t.get_time() << endl;
    // 确保前面不会有误
    // LOG(INFO) << "Occupy Core " << index << " finished: " << occupied[index] << endl;
}




bool Model9Partitioner::sub_get_free_vertex(vid_t &vid, vid_t index) {
    vid_t min = index * num_vertices_each_cores;
    vid_t idx = sub_dis(gen);

    for (size_t count = 0; count < num_vertices_each_cores; count++) {
        idx = (idx + count) % num_vertices_each_cores;
        vid = indices[min + idx];
        if (!vertex_adjacent_edges[vid].empty() && vertex_adjacent_edges[vid].size() < 2 * average_degree &&
            !is_cores[index].get(vid)) {
            if (acquire_vertex(vid)) {
                if (vertex_adjacent_edges[vid].empty() || vertex_adjacent_edges[vid].size() > 2 * average_degree ||
                    is_cores[index].get(vid)) { // 不符合条件
                    // TODO 探测邻居
                    release_vertex(vid);
                } else {
                    // acquired_partition[vid] = index;
                    return true;
                }
            }
        }
    }
    return false;
}

bool Model9Partitioner::get_free_vertex(vid_t &vid) {
    //随机选择一个节点
    vid = dis(gen);
    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    //TODO 什么叫已经超出平衡范围
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices &&
           (vertex_adjacent_edges[vid].size() == 0 || vertex_adjacent_edges[vid].size() > 2 * average_degree || is_cores[current_partition].get(vid))) {
        vid = (vid + ++count) % num_vertices;
    }
    if (count == num_vertices)
        return false;
    return true;
}

void Model9Partitioner::sub_add_boundary(vid_t vid, size_t index) {
    auto &is_core = is_cores[index];
    auto &is_boundary = is_boundaries[index];
    // CHECK(!is_boundary.get(vid)) << vid << " is in boundary";
    if (is_boundary.get(vid)) return;
    is_boundary.set_bit_unsync(vid);

    CHECK(!is_core.get(vid)) << vid << " is in core";

    size_t degree = vertex_adjacent_edges[vid].size();
    CHECK_GE(degree, 0) << "Degree is 0";
    min_heaps[index].insert(degree, vid);

    auto iterator = vertex_adjacent_edges[vid].begin();
    while (iterator != vertex_adjacent_edges[vid].end()) {
        edge_t &edge = *iterator->second;
        size_t e_id = iterator->first;
        if (edge.valid()) {
            vid_t neighbor = edge.first == vid ? edge.second : edge.first;
            if (is_core.get(neighbor) && occupied[index] < capacity * capacity_ratios[index]) {
                assign_edge(index, edge.first, edge.second);
                min_heaps[index].decrease_key(vid);
                vertex_adjacent_edges[vid].erase(iterator++);
                edge.remove();
                edge.set_partition(index);
                // vertex_adjacent_edges[neighbor].erase(e_id);
            } else if (is_boundary.get(neighbor) && occupied[index] < capacity * capacity_ratios[index]) {
                assign_edge(index, edge.first, edge.second);
                min_heaps[index].decrease_key(vid);
                min_heaps[index].decrease_key(neighbor);
                vertex_adjacent_edges[vid].erase(iterator++);
                vertex_adjacent_edges[neighbor].erase(e_id);
                edge.remove();
                edge.set_partition(index);
            } else {
                iterator++;
            }
        } else {
            iterator++;
        }
    }
}


void Model9Partitioner::add_boundary(vid_t vid) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[current_partition];
    auto &is_boundary = is_boundaries[current_partition];

    // 如果已经被加入到边界集，直接返回
    if (is_boundary.get(vid)) return;
    // 把顶点加入到边界集
    is_boundary.set_bit_unsync(vid);

    CHECK(!is_core.get(vid)) << vid << " is in core";

    size_t degree = vertex_adjacent_edges[vid].size();
    CHECK_GE(degree, 0) << "Degree is 0";
    // TODO 感觉是这个的问题？
    min_heaps[current_partition].insert(degree, vid);

    auto iterator = vertex_adjacent_edges[vid].begin();
    while (iterator != vertex_adjacent_edges[vid].end()) {
        edge_t &edge = *iterator->second;
        size_t e_id = iterator->first;
        if (edge.valid()) {
            vid_t neighbor = edge.first == vid ? edge.second : edge.first;
            if (is_core.get(neighbor) && occupied[current_partition] < capacity * balance_ratio) {
                assign_edge(current_partition, edge.first, edge.second);
                min_heaps[current_partition].decrease_key(vid);
                vertex_adjacent_edges[vid].erase(iterator++);
                edge.remove();
            } else if (is_boundary.get(neighbor) && occupied[current_partition] < capacity * balance_ratio) {
                assign_edge(current_partition, edge.first, edge.second);
                min_heaps[current_partition].decrease_key(vid);
                min_heaps[current_partition].decrease_key(neighbor);
                vertex_adjacent_edges[vid].erase(iterator++);
                vertex_adjacent_edges[neighbor].erase(e_id);
                edge.remove();
            } else {
                iterator++;
            }
        } else {
            iterator++;
        }
    }
}

void Model9Partitioner::assign_edge(size_t index, vid_t from, vid_t to) {

    // edge_partition[e_id] = index;
    // save_edge(from, to, num_partitions);
    true_vids.set_bit_unsync(from);
    true_vids.set_bit_unsync(to);
    // TODO 并发问题，这里为什么会有并发问题？每个顶点不是只被一个线程持有吗？
    is_mirrors[from].set_bit_unsync(index);
    is_mirrors[to].set_bit_unsync(index);

    reverse_is_mirrors[index].set_bit_unsync(from);
    reverse_is_mirrors[index].set_bit_unsync(to);

    assigned_edges.fetch_add(1);
    // LOG(INFO) << "Assign edge " << from << " " << to << " to partition " << index;
    // LOG(INFO) << "Assigned edges: " << assigned_edges << endl;
    occupied[index]++;
//    // TODO 这两个不太重要
//    degrees[from]--;
//    degrees[to]--;
}

size_t Model9Partitioner::check_edge(const edge_t *e) {
    rep (i, num_partitions) {
        auto &is_boundary = is_boundaries[i];
        if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
            occupied[i] < capacity) {
            return i;
        }
    }

    rep (i, num_partitions) {
        auto &is_core = is_cores[i], &is_boundary = is_boundaries[i];
        if ((is_core.get(e->first) || is_core.get(e->second)) &&
            occupied[i] < capacity) {
            if (is_core.get(e->first) && degrees[e->second] > average_degree)
                continue;
            if (is_core.get(e->second) && degrees[e->first] > average_degree)
                continue;
            is_boundary.set_bit(e->first);
            is_boundary.set_bit(e->second);
            return i;
        }
    }
    return num_partitions;
}

void Model9Partitioner::build_vertex_adjacent_edges() {
    // 遍历edges
    for (size_t i = 0; i < edges.size(); i++) {
        if (edges[i].valid()) {
            vid_t from = edges[i].first;
            vid_t to = edges[i].second;
            vertex_adjacent_edges[from].emplace(i, &edges[i]);
            vertex_adjacent_edges[to].emplace(i, &edges[i]);
        }
    }
}

//void Model9Partitioner::calculate_replication_factor() {
//    LOG(INFO) << "Calculating replication factor..." << endl;
//    // 对每个分区的边界集求和
//    for (size_t i = 0; i < reverse_is_mirrors.size(); i++) {
//        replicas += is_boundaries[i].popcount();
//    }
//    replication_factor = replication_factor = (double) replicas / (double) true_vids.popcount();
//}