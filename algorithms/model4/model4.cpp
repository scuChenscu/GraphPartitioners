#include "model4.hpp"
//固定随机数
// 构造函数
Model4Partitioner::Model4Partitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                                     size_t num_partitions) : EdgePartitioner(baseGraph, algorithm, num_partitions), input(input), gen(985) {
    config_output_files();
    cores = thread::hardware_concurrency();
    total_time.start();
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
    visited = dense_bitset(num_vertices);
    indices.resize(num_vertices);
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
    // 初始化的时候构造图
    adj_out.build(edges);
    // 存储反向边
    adj_in.build_reverse(edges);

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(input), std::ios::binary);
    degree_file.read((char *) &degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();

    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);
    fin.close();
}

//最后一个子图就是剩下边组合而成
void Model4Partitioner::assign_remaining() {
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

void Model4Partitioner::assign_master() {
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

size_t Model4Partitioner::count_mirrors() {
    size_t result = 0;
    rep(i, num_partitions) result += is_boundaries[i].popcount();
    return result;
}

void Model4Partitioner::split() {
    // 构造顶点的邻接表
    LOG(INFO) << "construct_adj_list" << endl;
    // TODO，感觉这个不一定需要，因为adj_out和adj_in已经存储了所有边的信息
    // 这里的adj_list是直接根据邻接表拿到邻居顶点
    // 重新索引
    LOG(INFO) << "re_index" << endl;
    re_index();
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
    LOG(INFO) << "Start Model4 partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "Model4"
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
        threads[i] = thread(&Model4Partitioner::sub_split,this, i);
    }
    for (int i = 0; i < cores; ++i) {
        threads[i].join();
    }


    // TODO 先保证单个分区不出现错误
    // LOG(INFO) << "Start sub_split" << endl;
    // sub_split(0);
    LOG(INFO) << "Model4 multi-partitioning finished!" << endl;
    return;
    // TODO 资源清理
//    min_heaps.clear();
//    vertex_ofstream.close();
//    edge_ofstream.close();
//    adjacency_list.clear();
//    visited.clear();
//
//    part_degrees.clear();
//    is_cores.clear();
//    is_boundaries.clear();
//    is_mirrors.clear();
//    true_vids.clear();
//    min_heap.clear();
//    edges.clear();
//    degrees.clear();
//    delete &adj_out;
//    // delete &adj_in;
//    // v_set.clear();

    // return;
    min_heap = min_heaps[0];
    // 前p-1个分区
    LOG(INFO) << "Start building the first num_partitions-1 full partition" << endl;
    for (current_partition = 0; current_partition < num_partitions - 1; current_partition++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid;
            if (!min_heap.get_min(degree, vid)) { // 当S\C为空时，从V\C中随机选择顶点
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                    break;
                }
                // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
                degree = adj_out[vid].size() + adj_in[vid].size();
            } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
                min_heap.remove(vid);
                // d.remove(vid);
            }
            // 把顶点加入到C，即核心集
            occupy_vertex(vid, degree);
        }
        //TODO 清空最小堆
        // LOG(INFO) << "Clear min heap" << endl;
        min_heap.clear();

        //TODO 感觉不应该在这里操作，应该在分配顶点的时候就更新
        rep(direction, 2) repv(vid, num_vertices) {
                // 获取所有顶点的邻接表
                adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size();) {
                    if (edges[neighbors[i].v].valid()) {
                        i++;
                    } else {
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    }
                }
            }
    }
    LOG(INFO) << "Start building the last full partition" << endl;
    current_partition = num_partitions - 1;
    // 把剩余的边放入最后一个分区
    assign_remaining();

    // CHECK_EQ(assigned_edges, num_edges);

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
    }

    ifstream fin(binary_edgelist_name(input), std::ios::binary | std::ios::ate);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
    edges.resize(num_edges);
    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
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
void Model4Partitioner::sub_split(size_t index) {
    LOG(INFO) << "Start sub_split " << index << endl;
    // 根据i的值从v_set中获取指定范围的顶点集合，选择不同的起始点
    // 然后在该顶点集合上执行NE算法
    // 当前分区的边数小于负载上限时，添加顶点到核心集C
    // MinHeap<vid_t, vid_t> cur_min_heap = min_heaps[num_partitions_i];
    // 设定一个上限
    // occupied存储的是每个分区边的数量
    while (occupied[index] <= capacity * CAPACITY_RATIO) {
        // LOG(INFO) << "Core " << index << " occupied edge count: " << occupied[index] << endl;
        vid_t degree, vid;
        // 尝试从当前分区所属的最小堆中获取顶点
        if (!min_heaps[index].get_min(degree, vid)) { // 当S\C为空时，从V\C中随机选择顶点
            // LOG(INFO) << "sub_get_free_vertex" << endl;
            if (!sub_get_free_vertex(vid, index)) { // 当V\C已经没有顶点，结束算法
                break;
            }
            // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
            degree = adj_out[vid].size() + adj_in[vid].size();
        } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
            min_heaps[index].remove(vid);
        }
        // 把顶点加入到C，即核心集
        sub_occupy_vertex(vid, degree, index);
    }
    // 确保前面不会有误
    LOG(INFO) << "Occupy Core "<< index << " finished: " << occupied[index] << endl;

    //TODO 因为每个分区独立使用min_heap,不需要清空最小堆
//    min_heaps[current_partition].clear();
//    is_cores[current_partition].clear();
//    is_boundaries[current_partition].clear();
    //TODO 为什么要全局扫描，移除已经分配的邻边
    rep(direction, 2) repv(vid, num_vertices) {
            // 获取所有顶点的邻接表
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            for (size_t i = 0; i < neighbors.size();) {
                if (edges[neighbors[i].v].valid()) {
                    i++;
                } else {
                    // 将边移动到最后
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
}

void Model4Partitioner::re_index()  {
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

bool Model4Partitioner::sub_get_free_vertex(vid_t &vid, vid_t index)  {
    // 前面根据cores把顶点划分成了几个簇
    // 根据index和cores来判断当前index属于indices的[min, max]
    //随机选择一个节点
    vid_t min = index * cores;
    vid_t idx = sub_dis(gen); // 生成 1000 / 10 = 100；0-99
    // TODO 是否需要 - 1
    vid_t max = (index + 1) * cores - 1;
    // 选择数据范围在[min, max]

    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    vid = indices[min + idx];

    //TODO 什么叫已经超出平衡范围
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices_each_cores &&
           (adj_out[vid].size() + adj_in[vid].size() == 0 ||
            adj_out[vid].size() + adj_in[vid].size() >
            2 * average_degree ||
            is_cores[index].get(vid))) { // 此时候选集跟核心集是一致的，只需要判断一个即可
        // 应该是一个链式寻找的过程
        idx = (idx + ++count) % num_vertices_each_cores;
        vid = indices[min + idx];
    }
    if (count == num_vertices_each_cores)
        return false;
    return true;
}

bool Model4Partitioner::get_free_vertex(vid_t &vid) {
    //随机选择一个节点
    vid = dis(gen);
    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    //TODO 什么叫已经超出平衡范围
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices &&
           (adj_out[vid].size() + adj_in[vid].size() == 0 ||
            adj_out[vid].size() + adj_in[vid].size() >
            2 * average_degree ||
            is_cores[current_partition].get(vid))) {
        vid = (vid + ++count) % num_vertices;
    }
    if (count == num_vertices)
        return false;
    return true;
}

void Model4Partitioner::sub_occupy_vertex(vid_t vid, vid_t d, size_t index)  {
    CHECK(!is_cores[index].get(vid)) << "add " << vid << " to core again";
    // 核心集是vector<dense_bitset>，dense_bitset是一个稠密位图
    // 对位图的vid位置置1，表示vid被分配到num_partitions分区

    // 1. 把顶点加入核心集和边界集，不需要考虑重复设置的场景
    is_cores[index].set_bit_unsync(vid);
    // is_boundaries[partition].set_bit_unsync(vid);
    // 顶点的来源有两种情况，一是从V，二是从S\C
    // 如果顶点的度为0，不需要处理
    if (d == 0) return;
    sub_add_boundary(vid, index);

    // 2. 遍历vid的邻居顶点，把不属于边界集的顶点加入到边界集
    for (auto &i: adj_out[vid]) {
        if (edges[i.v].valid()) { // 因为是边分区算法，所以我们不能引入重复的边，只需要考虑未分配的边
            sub_add_boundary(edges[i.v].second, index);
        }
    }

    adj_out[vid].clear();

    for (auto &i: adj_in[vid]) {
        if (edges[i.v].valid())
            sub_add_boundary(edges[i.v].first, index);
    }
    adj_in[vid].clear();
}

void Model4Partitioner::occupy_vertex(vid_t vid, vid_t d)  {
    CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
    // 核心集是vector<dense_bitset>，dense_bitset是一个稠密位图
    // 对位图的vid位置置1，表示vid被分配到num_partitions分区
    is_cores[current_partition].set_bit_unsync(vid);
    // 如果顶点的度为0，不需要处理
    if (d == 0) return;

    // 这一步的意义是什么？
    add_boundary(vid);

    // 将顶点的邻居加入到边界集，这里需要过滤掉核心集中的顶点
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

void Model4Partitioner::sub_add_boundary(vid_t vid, size_t index) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[index];
    auto &is_boundary = is_boundaries[index];
    // 当从S\C引入顶点时，它的邻居顶点可能已经加入到边界集
    if (is_boundary.get(vid)) return;
    // 2. 如果顶点不在边界集，把顶点加入到边界集
    is_boundary.set_bit_unsync(vid);
    // 当我们引入新的顶点到边界集的时候，需要把它的度信息加入到最小堆
    if (!is_core.get(vid)) {
        min_heaps[index].insert(adj_out[vid].size() + adj_in[vid].size(), vid);

    }

    // 3. 把顶点的在边界集的邻居顶点的边加入到当前边集
    rep (direction, 2) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        // 遍历顶点的邻边
        for (size_t i = 0; i < neighbors.size();) {
            if (edges[neighbors[i].v].valid()) {
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                if (is_core.get(u)) { // 如果顶点在核心集中
                    assign_edge(index, direction ? vid : u,
                                direction ? u : vid);
                    min_heaps[index].decrease_key(vid); // 默认移除一条边
                    edges[neighbors[i].v].remove();
                    //TODO 交换到最后位置，然后长度减1
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else if (is_boundary.get(u) &&
                           occupied[index] < capacity * CAPACITY_RATIO) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                    // 将边加入边集
                    assign_edge(index, direction ? vid : u, direction ? u : vid);
                    min_heaps[index].decrease_key(vid);
                    min_heaps[index].decrease_key(u);
                    edges[neighbors[i].v].remove();
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else
                    i++;
            } else {
                //swap是pop的前提，先交换到最后位置然后把长度减1
                std::swap(neighbors[i], neighbors.back());
                neighbors.pop_back();
            }
        }
    }
}

void Model4Partitioner::add_boundary(vid_t vid) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[current_partition];
    auto &is_boundary = is_boundaries[current_partition];

    // 如果已经被加入到边界集，直接返回
    if (is_boundary.get(vid))
        return;
    // 把顶点加入到边界集
    is_boundary.set_bit_unsync(vid);

    // 如果顶点没有在核心集中，直接把顶点的度数据加入到最小堆
    // 因为在核心集中的顶点已经不需要处理了
    if (!is_core.get(vid)) {
        min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
    }

    //仅支持无向图，在计算neighbor的时候有向和无向会导致邻居的差别从而影响分割

    rep (direction, 2) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        // 遍历顶点的邻边
        for (size_t i = 0; i < neighbors.size();) {
            if (edges[neighbors[i].v].valid()) {
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                if (is_core.get(u)) { // 如果顶点在核心集中
                    assign_edge(current_partition, direction ? vid : u,
                                direction ? u : vid);
                    // TODO 这个要修改
                    min_heap.decrease_key(vid); // 默认移除一条边
                    edges[neighbors[i].v].remove();
                    //TODO 交换到最后位置，然后长度减1
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else if (is_boundary.get(u) &&
                           occupied[current_partition] < capacity) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                    // 将边加入边集
                    assign_edge(current_partition, direction ? vid : u, direction ? u : vid);
                    min_heap.decrease_key(vid);
                    min_heap.decrease_key(u);
                    edges[neighbors[i].v].remove();
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else
                    i++;
            } else {
                //swap是pop的前提，先交换到最后位置然后把长度减1
                std::swap(neighbors[i], neighbors.back());
                neighbors.pop_back();
            }
        }
    }
}

void Model4Partitioner::assign_edge(size_t index, vid_t from, vid_t to)  {
    // save_edge(from, to, num_partitions);
    true_vids.set_bit_unsync(from);
    true_vids.set_bit_unsync(to);
    is_mirrors[from].set_bit_unsync(index);
    is_mirrors[to].set_bit_unsync(index);
    assigned_edges++;
    occupied[index]++;
    degrees[from]--;
    degrees[to]--;
}

size_t Model4Partitioner::check_edge(const edge_t *e) {
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