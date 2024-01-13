#include "model5.hpp"

//固定随机数
// 构造函数
Model5Partitioner::Model5Partitioner(BaseGraph &baseGraph, const string &input, const string &algorithm,
                                     size_t num_partitions) : EdgePartitioner(baseGraph, algorithm, num_partitions),
                                                              input(input), gen(985) {
    config_output_files();
    std::ifstream fin(binary_edgelist_name(input),
                      std::ios::binary | std::ios::ate);
    // tellp 用于返回写入位置，
    // tellg 则用于返回读取位置也代表着输入流的大小
    auto filesize = fin.tellg();
//  LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    //最大的下标+1
    fin.read((char *) &num_vertices, sizeof(num_vertices));
    fin.read((char *) &num_edges, sizeof(num_edges));

    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
             filesize);
    avg_vertices_each_partition = num_partitions / cores;
    current_partition = 0;
    average_degree = (double) num_edges * 2 / num_vertices;
    assigned_edges = 0;
    capacity = (double) num_edges * BALANCE_RATIO / num_partitions + 1;
    occupied.assign(num_partitions, 0);
    num_vertices_in_partition.assign(num_partitions, 0);
    num_vertices_each_cores = num_vertices / cores;
    avg_vertices_each_partition = num_vertices / num_partitions;
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    adj_directed.resize(num_vertices);
    visited = dense_bitset(num_vertices);
    assigned = dense_bitset(num_edges);
    edge_partition.resize(num_edges);
    vertex_lock = dense_bitset(num_vertices);
    indices.resize(num_vertices);
    reverse_indices.resize(num_vertices);
    v_lock.resize(num_vertices, 0);
    // 每个都分配num_vertices大小，保证不会有问题
    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundaries.assign(num_partitions, dense_bitset(num_vertices));
    true_vids.resize(num_vertices);
    is_mirrors.assign(num_vertices, dense_bitset(num_partitions));
    master.assign(num_vertices, -1);
    dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    sub_dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices_each_cores - 1));

    edges.resize(num_edges);
    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);

    vertex_adjacent_edges.resize(num_vertices, unordered_map<size_t, edge_t *>());

    // 初始化的时候构造图
    // adj_out.build(edges);
    // 存储反向边
    // adj_in.build_reverse(edges);

    adj_out = graph.adj_out;
    adj_in = graph.adj_in;

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(input), std::ios::binary);
    degree_file.read((char *) &degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();

    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);
    fin.close();
}

//最后一个子图就是剩下边组合而成
void Model5Partitioner::assign_remaining() {
    auto &is_boundary = is_boundaries[num_partitions - 1];
    auto &is_core = is_cores[num_partitions - 1];

    // 遍历边集合
    for (auto &e: edges) {
        if (e.valid()) {
            assign_edge(num_partitions - 1, e.first, e.second);
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
}

void Model5Partitioner::assign_master() {
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

size_t Model5Partitioner::count_mirrors() {
    size_t result = 0;
    rep(i, num_partitions) result += is_boundaries[i].popcount();
    return result;
}

void Model5Partitioner::split() {
    LOG(INFO) << "build_vertex_adjacent_edges" << endl;

    build_vertex_adjacent_edges();
    assigned.set_bit(10);
    // 构造顶点的邻接表
    LOG(INFO) << "construct_adj_list" << endl;
    // TODO，感觉这个不一定需要，因为adj_out和adj_in已经存储了所有边的信息
    // 这里的adj_list是直接根据邻接表拿到邻居顶点
    // 重新索引
    LOG(INFO) << "re_index" << endl;
    re_index();
    // 根据边的两端顶点所属的cores建立有向边
    LOG(INFO) << "Build direction edge" << endl;
    adj_directed.build_directed(edges, reverse_indices);
    // 每个分区的理论顶点数
    // 0, num_vertices/p, 2*num_vertices/p, 3*num_vertices/p,...
    // p个分区，每个分区的顶点范围是[num_partitions * i, p * (i + 1) - 1]

    // 初始化最小堆，用于存储S\C的顶点信息
    min_heap.reserve(num_vertices);

    min_heaps.resize(num_partitions, MinHeap<vid_t, vid_t>(num_vertices));
    // TODO 存储每个顶点的度，暂时不需要
//    d.reserve(num_vertices);
//    repv(vid, num_vertices) {
//        d.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
//    }
    LOG(INFO) << "Start Model5 partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "Model5"
       << endl
       << "基于BFS对顶点进行分块，使用K-way多线程来同时进行k个分区划分，最后将未分配的顶点和边统一收集处理"
       << "目标是提高NE算法的运行速度，同时保证RF尽可能低"
       << endl
       << "BALANCE RATIO:" << BALANCE_RATIO
       << " | Cores: " << cores
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());
    // TODO 这个不是根据分区数，而是根据机器的核心数来确定线程数
    // const size_t numThreads = c;
    thread threads[cores];
    // 启动线程
    for (int i = 0; i < cores; ++i) {
        threads[i] = thread(&Model5Partitioner::sub_split, this, i);
    }
    total_time.start();
    for (int i = 0; i < cores; ++i) {
        threads[i].join();
    }


    // TODO 先保证单个分区不出现错误
    // LOG(INFO) << "Start sub_split" << endl;
    // sub_split(0);
    LOG(INFO) << "Model5 multi-partitioning finished!" << endl;
    // total_time.stop();
    // LOG(INFO) << "Total time: " << total_time.get_time() << endl;
    // return;
    // 前p-1个分区
    LOG(INFO) << "Start building the first num_partitions - 1 full partition" << endl;
    // 用位图来维护每个分区新引入的顶点
    dirty_vertices.resize(num_vertices);
    for (current_partition = 0; current_partition < num_partitions - 1; current_partition++) {
        // LOG(INFO) << "Start partition " << current_partition  << " edge count: " << occupied[current_partition] << endl;
        min_heap = min_heaps[current_partition];
        LOG(INFO) << "Min_heap " << current_partition << "  size: " << min_heap.size() << endl;
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid;

            bool has_min = false;
            while (min_heap.size() > 0) {
                min_heap.get_min(degree, vid);
                min_heap.remove(vid);
                if (!vertex_adjacent_edges[vid].empty()) {
                    // 更新degree
                    degree = vertex_adjacent_edges[vid].size();
                    has_min = true;
                    break;
                }
            }
            if (!has_min) { // 当S\C为空时，从V\C中随机选择顶点
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                    break;
                }
                // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
                degree = adj_directed[vid].size();
            }
            // 把顶点加入到C，即核心集
            CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
            is_cores[current_partition].set_bit_unsync(vid);
            is_boundaries[current_partition].set_bit_unsync(vid);

            for (auto &pair: vertex_adjacent_edges[vid]) {
                edge_t &edge = *pair.second;
                if (edge.valid()) {
                    vid_t neighbor = edge.first == vid ? edge.second : edge.first;
                    if (edge.valid()) {
                        add_boundary(neighbor);
                    }
                }
            }
        }
        min_heap.clear();
        LOG(INFO) << "End partition " << current_partition  << " edge count: " << occupied[current_partition] << endl;
    }
    LOG(INFO) << "Start building the last full partition" << endl;
    current_partition = num_partitions - 1;
    // 把剩余的边放入最后一个分区
    assign_remaining();
    LOG(INFO) << "End partition " << current_partition  << " edge count: " << occupied[current_partition] << endl;

    CHECK_EQ(assigned_edges, num_edges);




    //根据结点平衡性、随机分配的重叠度以及结点的度大小来判断
    vector<vid_t> current_partitions(num_partitions);
    capacity = (double) true_vids.popcount() * 1.05 / num_partitions + 1;

    // 复制因子
    rep(i, num_vertices) {
        double max_score = 0.0;
        vid_t which_p;
        bool unique = false;
        // 判断顶点是否只属于一个分区
        if (is_mirrors[i].popcount() == 1) {
            unique = true;
        }
        // 计算顶点的分区最高得分
        repv(j, num_partitions) {
            if (is_mirrors[i].get(j)) {
                num_vertices_in_partition[j]++; // 每个分区的顶点数
                double score = (part_degrees[i][j] / (degrees[i] + 1)) + (current_partitions[j] < capacity ? 1 : 0);
                if (unique) {
                    which_p = j;
                    break;
                } else if (max_score < score) {
                    max_score = score;
                    which_p = j;
                }
            }
        }
        // 用最高得分所在的分区作为顶点的最终分区
        ++current_partitions[which_p];
        save_vertex(i, which_p);
        balance_vertex_distribute[i] = which_p;
    }
    // vertex_ofstream.close();
    repv(j, num_partitions) {
        LOG(INFO) << "Partition " << j << " Vertex Count: " << current_partitions[j];
        LOG(INFO) << "Partition " << j << " Edge Count: " << occupied[j];
    }

    ifstream fin(binary_edgelist_name(input), std::ios::binary | std::ios::ate);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
    edges.resize(num_edges);
    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
    size_t idx = 0;
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];

        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
    }
    // edge_ofstream.close();

    total_time.stop();

    LOG(INFO) << "total partition time: " << total_time.get_time() << endl;


}

// TODO 这里要注意不能用current_partition，current_partition是全局变量，只能用在后续建立完整分区
void Model5Partitioner::sub_split(size_t index) {
    LOG(INFO) << "Start sub_split " << index << endl;
    while (occupied[index] <= capacity * CAPACITY_RATIO) {
        vid_t degree, vid;
        // 尝试从当前分区所属的最小堆中获取顶点
        if (!min_heaps[index].get_min(degree, vid)) { // 当S\C为空时，从V\C中随机选择顶点
            // LOG(INFO) << "sub_get_free_vertex" << endl;
            if (!sub_get_free_vertex(vid, index)) { // 当V\C已经没有顶点，结束算法
                break;
            }
            degree = vertex_adjacent_edges[vid].size();
        } else {
            min_heaps[index].remove(vid);
            // LOG(INFO) << "Min_heaps " << index << " remove: " << vid << " degree: " << degree << endl;
        }
        // sub_occupy_vertex(vid, degree, index);
        CHECK(!is_cores[index].get(vid)) << "add " << vid << " to core again";
        is_cores[index].set_bit_unsync(vid);
        is_boundaries[index].set_bit_unsync(vid);

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
        // LOG(INFO) << "Min_heap " << index << "  size: " << min_heaps[index].size() << endl;
        release_vertex(vid);
    }
    // 确保前面不会有误
    LOG(INFO) << "Occupy Core " << index << " finished: " << occupied[index] << endl;
    // LOG(INFO) << "Min_heap " << index << "  size: " << min_heaps[index].size() << endl;
}


void Model5Partitioner::re_index() {
    queue<vid_t> v_queue;
    auto start = std::chrono::high_resolution_clock::now(); // 记录开始时间
    // 随机选择顶点，进行广度遍历，重新索引
    vid_t index = 0;
    vid_t vid = dis(gen);
    // 基于该顶点进行深度遍历，对每个顶点重新索引
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
        set < vid_t > neighbor_set = adjacency_list.find(v)->second;
        // 将neighbor_set加入v_queue和v_set中
        for (auto &i: neighbor_set) {
            v_queue.push(i);
        }
    }
    auto end = std::chrono::high_resolution_clock::now(); // 记录结束时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // 计算时间差
    LOG(INFO) << "re_index time: " << duration.count() << "ms" << endl;
}

bool Model5Partitioner::sub_get_free_vertex(vid_t &vid, vid_t index) {
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
                    release_vertex(vid);
                } else {
                    return true;
                }
            }
        }
    }
    return false;
}

bool Model5Partitioner::get_free_vertex(vid_t &vid) {
    //随机选择一个节点
    vid = dis(gen);
    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    //TODO 什么叫已经超出平衡范围
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices &&
           (adj_directed[vid].size() == 0 || adj_directed[vid].size() > 2 * average_degree || is_cores[current_partition].get(vid))) {
        vid = (vid + ++count) % num_vertices;
    }
    if (count == num_vertices)
        return false;
    return true;
}

void Model5Partitioner::sub_occupy_vertex(vid_t vid, vid_t d, size_t index) {
    CHECK(!is_cores[index].get(vid)) << "add " << vid << " to core again";
    is_cores[index].set_bit_unsync(vid);
    if (d == 0) return;
    sub_add_boundary(vid, index);

    for (auto &i: adj_directed[vid]) {
        if (edges[i.v].valid()) {
            vid_t neighbor = edges[i.v].second == vid ? edges[i.v].first : edges[i.v].second;
            if (acquire_vertex(neighbor)) {
                if (edges[i.v].valid()) {
                    // TODO 这里记录了指向顶点的边，直接分配，因为neighbor没有持有该边
                    assign_edge(index, vid, neighbor);
                    edges[i.v].remove();
                    sub_add_boundary(neighbor, index);
                } else {
                    release_vertex(neighbor);
                }
            }
        }
    }
}


void Model5Partitioner::occupy_vertex(vid_t vid, vid_t d) {
    CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
    is_cores[current_partition].set_bit_unsync(vid);
    if (d == 0) return;

    add_boundary(vid);

    for (auto &i: adj_directed[vid]) {
        if (edges[i.v].valid()) {
            vid_t neighbor = edges[i.v].second == vid ? edges[i.v].first : edges[i.v].second;
                if (edges[i.v].valid()) {
                    // TODO 这里记录了指向顶点的边，直接分配，因为neighbor没有持有该边
                    assign_edge(current_partition, vid, neighbor);
                    edges[i.v].remove();
                    sub_add_boundary(neighbor, current_partition);
                }
        }
    }
}

void Model5Partitioner::sub_add_boundary(vid_t vid, size_t index) {
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
            if (is_core.get(neighbor) && occupied[index] < capacity * CAPACITY_RATIO) {
                assign_edge(index, edge.first, edge.second);
                min_heaps[index].decrease_key(vid);
                vertex_adjacent_edges[vid].erase(iterator++);
                edge.remove();
                // vertex_adjacent_edges[neighbor].erase(e_id);
            } else if (is_boundary.get(neighbor) && occupied[index] < capacity * CAPACITY_RATIO) {
                assign_edge(index, edge.first, edge.second);
                min_heaps[index].decrease_key(vid);
                min_heaps[index].decrease_key(neighbor);
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


void Model5Partitioner::add_boundary(vid_t vid) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[current_partition];
    auto &is_boundary = is_boundaries[current_partition];

    // 如果已经被加入到边界集，直接返回
    if (is_boundary.get(vid))
        return;
    // 把顶点加入到边界集
    is_boundary.set_bit_unsync(vid);

    CHECK(!is_cores[current_partition].get(vid)) << vid << " is in core";

    CHECK(!is_cores[current_partition].get(vid)) << vid << " is in core";

    size_t degree = vertex_adjacent_edges[vid].size();
    CHECK_GE(degree, 0) << "Degree is 0";
    min_heaps[current_partition].insert(degree, vid);

    for (auto &pair: vertex_adjacent_edges[vid]) {
        edge_t &edge = *pair.second;
        size_t e_id = pair.first;
        if (edge.valid()) {
            vid_t neighbor = edge.first == vid ? edge.second : edge.first;
            if (is_core.get(neighbor) && occupied[current_partition] < capacity) {
                assign_edge(current_partition, edge.first, edge.second);
                min_heaps[current_partition].decrease_key(vid);
                dirty_vertices.set_bit_unsync(vid);
                edge.remove();
            } else if (is_boundary.get(neighbor) && occupied[current_partition] < capacity) {
                assign_edge(current_partition, edge.first, edge.second);
                min_heaps[current_partition].decrease_key(vid);
                dirty_vertices.set_bit_unsync(vid);
                min_heaps[current_partition].decrease_key(neighbor);
                dirty_vertices.set_bit_unsync(neighbor);
                edge.remove();
            }
        }
    }
}

void Model5Partitioner::assign_edge(size_t index, vid_t from, vid_t to) {
    // edge_partition[e_id] = index;
    // save_edge(from, to, num_partitions);
    true_vids.set_bit_unsync(from);
    true_vids.set_bit_unsync(to);
    // TODO 并发问题
    is_mirrors[from].set_bit_unsync(index);
    is_mirrors[to].set_bit_unsync(index);
//    is_mirrors[from].set_bit(index);
//    is_mirrors[to].set_bit(index);
    // TODO 这里需要原子操作，将assigned_edges改成原子类
    assigned_edges.fetch_add(1);
    // LOG(INFO) << "Assign edge " << from << " " << to << " to partition " << index;
    // LOG(INFO) << "Assigned edges: " << assigned_edges << endl;
    occupied[index]++;
    // TODO 这两个不太重要
    degrees[from]--;
    degrees[to]--;
}

size_t Model5Partitioner::check_edge(const edge_t *e) {
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

void Model5Partitioner::build_vertex_adjacent_edges() {
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

void Model5Partitioner::calculate_alpha() {
    // 对每个分区的边界集求和
    for (size_t i = 0; i < num_partitions; i++) {
        replicas += is_boundaries[i].popcount();
    }
    replication_factor = replication_factor = (double) replicas / (double) num_vertices;
}