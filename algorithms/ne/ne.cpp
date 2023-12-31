#include "ne.hpp"

//固定随机数
// 构造函数
NePartitioner::NePartitioner(std::string input, std::string algorithm, int num_partition)
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
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_cores.assign(p, dense_bitset(num_vertices));
    is_boundarys.assign(p, dense_bitset(num_vertices));
    true_vids.resize(num_vertices);
    is_mirrors.assign(num_vertices, dense_bitset(p));
    master.assign(num_vertices, -1);
    dis.param(
            std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));


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
void NePartitioner::assign_remaining() {
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

void NePartitioner::assign_master() {
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

size_t NePartitioner::count_mirrors() {
    size_t result = 0;
    rep(i, p) result += is_boundarys[i].popcount();
    return result;
}

void NePartitioner::split() {
    // 初始化最小堆，用于存储S\C的顶点信息
    min_heap.reserve(num_vertices);

    LOG(INFO) << "Start NE partitioning...";
    // 前p-1个分区
    for (bucket = 0; bucket < p - 1; bucket++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[bucket] < capacity) {
            vid_t d, vid;
            if (!min_heap.get_min(d, vid)) { // 当S\C为空时，从V\C中随机选择顶点
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                    break;
                }
                // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
                d = adj_out[vid].size() + adj_in[vid].size();
            } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
                min_heap.remove(vid);
            }
            // 把顶点加入到C，即核心集
            occupy_vertex(vid, d);
        }
        //TODO 清空最小堆
        min_heap.clear();

        //TODO 为什么要全局扫描，移除已经分配的邻边
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
    LOG(INFO) << "total partition time: " << total_time.get_time();
}

