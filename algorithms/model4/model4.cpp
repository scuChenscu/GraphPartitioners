#include "model4.hpp"
#include <thread>
//固定随机数
// 构造函数
Model4Partitioner::Model4Partitioner(const std::string& input, const std::string& algorithm, int num_partition)
        : input(input), gen(985) {
    p = num_partition;
    config_output_files(input, algorithm, num_partition);

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

    average_degree = (double) num_edges * 2 / num_vertices;
    assigned_edges = 0;
    capacity = (double) num_edges * BALANCE_RATIO / p + 1;
    occupied.assign(p, 0);
    num_vertices_in_partition.assign(p, 0);
    num_p_v = num_vertices / p;
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    visited = dense_bitset(num_vertices);
    indices.resize(num_vertices);
    is_cores.assign(p, dense_bitset(num_vertices));
    is_boundarys.assign(p, dense_bitset(num_vertices));
    true_vids.resize(num_vertices);
    is_mirrors.assign(num_vertices, dense_bitset(p));
    master.assign(num_vertices, -1);
    dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
    sub_dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_p_v - 1));

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

    part_degrees.assign(num_vertices, vector<vid_t>(p));
    balance_vertex_distribute.resize(num_vertices);
    fin.close();
}

//最后一个子图就是剩下边组合而成
void Model4Partitioner::assign_remaining() {
    auto &is_boundary = is_boundarys[p - 1], &is_core = is_cores[p - 1];
    repv(u, num_vertices) for (auto &i: adj_out[u])
            if (edges[i.v].valid()) {
                assign_edge(p - 1, u, edges[i.v].second);
                is_boundary.set_bit_unsync(u);
                is_boundary.set_bit_unsync(edges[i.v].second);
            }
    repv(i, num_vertices) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            //在其他子图中是核心集的不予理睬，不能设置为本子图的核心集
            rep(j, p - 1) if (is_cores[j].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
        }
    }
}

void Model4Partitioner::assign_master() {
    std::vector<vid_t> count_master(p, 0);
    std::vector<vid_t> quota(p, num_vertices);
    long long sum = p * num_vertices;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::vector<dense_bitset::iterator> pos(p);
    rep(b, p) pos[b] = is_boundarys[b].begin();
    vid_t count = 0;
    while (count < num_vertices) {
        long long r = distribution(gen) * sum;
        //随机选择哪个子图先赋值
        int k;
        for (k = 0; k < p; k++) {
            if (r < quota[k])
                break;
            r -= quota[k];
        }
        //选出当前位置还未被赋值的结点
        while (pos[k] != is_boundarys[k].end() && master[*pos[k]] != -1)
            pos[k]++;
        if (pos[k] != is_boundarys[k].end()) {
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
    rep(i, p) result += is_boundarys[i].popcount();
    return result;
}

void Model4Partitioner::split() {
    // 构造顶点的邻接表
    LOG(INFO) << "construct_adj_list" << endl;
    // TODO，感觉这个不一定需要，因为adj_out和adj_in已经存储了所有边的信息
    // 这里的adj_list是直接根据邻接表拿到邻居顶点
    construct_adj_list(edges);
    // 重新索引
    LOG(INFO) << "re_index" << endl;
    re_index();
    // 每个分区的理论顶点数
    // 0, num_vertices/p, 2*num_vertices/p, 3*num_vertices/p,...
    // p个分区，每个分区的顶点范围是[p * i, p * (i + 1) - 1]

    // 初始化最小堆，用于存储S\C的顶点信息
    min_heap.reserve(num_vertices);

    min_heaps.resize(p, MinHeap<vid_t, vid_t>(num_vertices));
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
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    const int numThreads = p;
    std::thread threads[numThreads];


    // 启动线程
//    for (int i = 0; i < numThreads; ++i) {
//        threads[i] = thread(&Model4Partitioner::sub_split,this, i);
//    }
//    for (int i = 0; i < numThreads; ++i) {
//        threads[i].join();
//    }


    // TODO 测试只处理一个分区
    sub_split(0);
    LOG(INFO) << "Model4 multi-partitioning finished!";
    min_heaps.clear();
    vertex_ofstream.close();
    edge_ofstream.close();
    adj_list.clear();
    visited.clear();

    part_degrees.clear();
    is_cores.clear();
    is_boundarys.clear();
    is_mirrors.clear();
    true_vids.clear();
    min_heap.clear();
    edges.clear();
    degrees.clear();
    delete &adj_out;
    // delete &adj_in;
    // v_set.clear();

    return;
    // min_heap = min_heaps[0];
    // 前p-1个分区
    for (bucket = 0; bucket < p - 1; bucket++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[bucket] < capacity) {
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
    bucket = p - 1;
    // 把剩余的边放入最后一个分区
    assign_remaining();

    CHECK_EQ(assigned_edges, num_edges);

    //根据结点平衡性、随机分配的重叠度以及结点的度大小来判断
    vector<vid_t> buckets(p);
    capacity = (double) true_vids.popcount() * 1.05 / p + 1;

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
        repv(j, p) {
            if (is_mirrors[i].get(j)) {
                num_vertices_in_partition[j]++; // 每个分区的顶点数
                double score = (part_degrees[i][j] / (degrees[i] + 1)) + (buckets[j] < capacity ? 1 : 0);
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
        ++buckets[which_p];
        save_vertex(i, which_p);
        balance_vertex_distribute[i] = which_p;
    }
    vertex_ofstream.close();
    repv(j, p) {
        LOG(INFO) << "each partition node count: " << buckets[j];
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
    edge_ofstream.close();

    total_time.stop();

    LOG(INFO) << "total partition time: " << total_time.get_time() << endl;
    stringstream result;

    calculate_replication_factor();
    calculate_alpha();
    calculate_rho();
    result << "Cost Time: " << total_time.get_time()
           << " | Replication Factor: " << replication_factor
           << " | Alpha: " << alpha
           << " | Replicas: " << replicas
           << " | Rho: " << rho
           << " | Rho Exclude Last Partition: " << rho_1
           << " | Max Edge: " << max_edge
           << " | Min Edge: " << min_edge
           << " | Avg Edge: " << num_edges / p
           << " | Edges: " << num_edges
           << " | Vertices: " << num_vertices
           //        << " | Max Vertex: " << max_vertex
           //        << " | Min Vertex: " << min_vertex
           << " | Avg Vertex: " << avg_vertex
           << " | Avg Vertex Exclude Last Partition: " << avg_vertex_1
           // << " | Beta: " << beta
           << endl;
    appendToFile(result.str());


}

void Model4Partitioner::sub_split(const int p_i) {
    // 根据i的值从v_set中获取指定范围的顶点集合，选择不同的起始点
    // 然后在该顶点集合上执行NE算法
    // 当前分区的边数小于负载上限时，添加顶点到核心集C
    // MinHeap<vid_t, vid_t> cur_min_heap = min_heaps[p_i];
    // 设定一个上限
    // occupied存储的是每个分区边的数量
    while (occupied[p_i] <= capacity * CAPACITY_RATIO) {
        vid_t degree, vid;
        if (!min_heaps[p_i].get_min(degree, vid)) { // 当S\C为空时，从V\C中随机选择顶点
            if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                break;
            }
            // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
            degree = adj_out[vid].size() + adj_in[vid].size();
        } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
            min_heaps[p_i].remove(vid);
        }
        // 把顶点加入到C，即核心集
        sub_occupy_vertex(vid, degree, p_i);
    }

    LOG(INFO) << "Occupy vertex finished." << endl;
    //TODO 因为每个分区独立使用min_heap,不需要清空最小堆
    min_heaps[p_i].clear();
    is_cores[p_i].clear();
    is_boundarys[p_i].clear();
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

