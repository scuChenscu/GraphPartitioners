#include "offstreamNH.hpp"

using namespace std;

//固定随机数
// 构造函数
OffstreamNHPartitioner::OffstreamNHPartitioner(BaseGraph& baseGraph, const string &input, const string &algorithm,
                             size_t num_partitions)
        : EdgePartitioner(baseGraph, algorithm, num_partitions), input(input), gen(985) {
    config_output_files();
    current_partition = 0;
    // average_degree = (double) num_edges * 2 / (double) num_vertices;
    assigned_edges = 0;

    partial_degree = baseGraph.partial_degree;

   // num_vertices_each_partition.assign(num_partitions, 0);

    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundaries.assign(num_partitions, dense_bitset(num_vertices));
    true_vids = baseGraph.true_vids;
    // master.assign(num_vertices, -1);
    dis.param(
            uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));
//
//    degrees.resize(num_vertices);
//    ifstream degree_file(degree_name(input), ios::binary);
//    degree_file.read((char *) &degrees[0], num_vertices * sizeof(vid_t));
//    degree_file.close();

    // part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    // balance_vertex_distribute.resize(num_vertices);

    stream_part = baseGraph.stream_part;
    off_part = baseGraph.off_part;

    vertex_partitions.assign(num_vertices, set<size_t>());

//    partial_degrees.resize(num_vertices);

    capacity =  BALANCE_RATIO * (off_part.size() / (double)num_partitions) + 1;

    max_partition_load = (uint64_t) BALANCE_RATIO * num_edges / num_partitions;
}

//最后一个子图就是剩下边组合而成
void OffstreamNHPartitioner::assign_remaining() {
    auto &is_boundary = is_boundaries[num_partitions - 1], &is_core = is_cores[num_partitions - 1];
    repv(u, num_vertices) for (auto &i: adj_out[u])
            if (edges[i.v].valid()) {
                assign_edge(num_partitions - 1, u, edges[i.v].second);
                vertex_partitions[u].insert(num_partitions - 1);
                vertex_partitions[edges[i.v].second].insert(num_partitions - 1);
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

size_t OffstreamNHPartitioner::count_mirrors() {
    size_t result = 0;
    rep(i, num_partitions) result += is_boundaries[i].popcount();
    return result;
}

void OffstreamNHPartitioner::split() {
    total_time.start();
    // 初始化最小堆，用于存储S\C的顶点信息
    min_heap.reserve(num_vertices);


    LOG(INFO) << "Start NE partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "OffstreamNH" << endl
       << "BALANCE RATIO:" << BALANCE_RATIO
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    // 前p-1个分区
    for (current_partition = 0; current_partition < num_partitions - 1; current_partition++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid;
            // LOG(INFO) << "Min_heap " << current_partition << "  size: " << min_heap.size() << endl;
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
            // 前面是获取顶点的两种方式，要么是从S\C，要么是从V\C
            // 把顶点加入到C，即核心集
            // 此时degree要么是从min_heap获取，要么是一个新顶点直接计算得到
            occupy_vertex(vid, degree);
        }
        //TODO 清空最小堆
        min_heap.clear();

        //TODO 为什么要全局扫描，移除已经分配的邻边
        // TODO 为什么要检查一遍
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
    repv(j, num_partitions) {
        LOG(INFO) << "Partition " << j << " Edge Count: " << occupied[j];
    }

    min_load = *min_element(occupied.begin(), occupied.end());
    max_load = *max_element(occupied.begin(), occupied.end());

    // 使用贪心来划分
    LOG(INFO) << "Start hdrf partitioning" << endl;
    LOG(INFO) << "max_partition_load: " << max_partition_load;


    for (auto &edge: stream_part) {
        int partition = find_max_score_partition(edge);
        is_mirrors[edge.first].set_bit_unsync(partition);
        is_mirrors[edge.second].set_bit_unsync(partition);
        occupied[partition]++;
        assigned_edges++;
        update_min_max_load(partition);
    }



    CHECK_EQ(assigned_edges, num_edges);
    total_time.stop();
}


int OffstreamNHPartitioner::find_max_score_partition(edge_t &e) {
    auto degree_u = ++partial_degree[e.first];
    auto degree_v = ++partial_degree[e.second];

    uint32_t sum;
    double max_score = 0;
    // TODO 这里默认分区为0
    uint32_t max_p = 0;
    double bal, gv, gu;

    for (int j = 0; j < num_partitions; j++) {
        // 如果分区已经超过了最大分区负载，直接跳过
        if (occupied[j] >= max_partition_load) {
            continue;
        }
        if (max_p == num_partitions) {
            max_p = j;
        }
        // 以下对应着hdrf算法的核心实现
//        gu = 0, gv = 0;
//        sum = degree_u + degree_v;
//        // 归一化
//        if (is_mirrors[e.first].get(j)) {
//            gu = degree_u;
//            gu /= sum;
//            gu = 1 + (1 - gu);
//        }
//        if (is_mirrors[e.second].get(j)) {
//            gv = degree_v;
//            gv /= sum;
//            gv = 1 + (1 - gv);
//        }
//        double rep = gu + gv; // rep值
//        bal = max_load - occupied[j];
//        if (min_load != UINT64_MAX) {
//            bal /= (epsilon + max_load - min_load);
//        }
        // 计算结果应该有两部分组成，rep和bal
        // LOG(INFO) << "rep: " << rep << " bal: " << lambda * bal;
        double score_p = calculate_rf_score(e.first, e.second, j) + lambda * calculate_lb_score(j);
        // LOG(INFO) << "score_p: " << score_p;
        CHECK_GE(score_p, 0) << "score_p: " << score_p;
        if (score_p > max_score) {
            max_score = score_p;
            max_p = j;
        }
    }
    return max_p;
}

double OffstreamNHPartitioner::calculate_rf_score(vid_t u, vid_t v, size_t partition_id) {
    double gu = 0, gv = 0;
    size_t u_degree = partial_degree[u];
    size_t v_degree = partial_degree[v];
    size_t sum = v_degree + u_degree;
    // 归一化
    if (is_mirrors[u].get(partition_id)) {
        gu = (double)u_degree;
        gu /= (double)sum;
        gu = 1 + (1 - gu);
    }
    if (is_mirrors[v].get(partition_id)) {
        gv = (double)v_degree;
        gv /= (double)sum;
        gv = 1 + (1 - gv);
    }
    double rf_score = gu + gv;
    return rf_score;
}

double OffstreamNHPartitioner::calculate_lb_score(size_t partition_id) {

    double  lb_score = (double)max_load - occupied[partition_id];
    if (min_load != UINT64_MAX) {
        lb_score /= (epsilon + (double)max_load - (double)min_load);
    }

    return lb_score;
}


bool OffstreamNHPartitioner::get_free_vertex(vid_t &vid) {
    //随机选择一个节点
    size_t index = dis(gen);
    vid  = index;
    vid_t count = 0; // 像是一个随机数，用来帮助选择随机顶点
    //TODO 什么叫已经超出平衡范围
    //如果是孤立节点直接跳过，或者当前结点在当前分割图中已超出平衡范围继续寻找，或者已经是核心集的结点
    while (count < num_vertices &&
           (adj_out[vid].size() + adj_in[vid].size() == 0 ||
            adj_out[vid].size() + adj_in[vid].size() >
            2 * avg_degree ||
            is_cores[current_partition].get(vid))) {
        vid = (index + ++count) % num_vertices;
    }
    if (count == num_vertices)
        return false;
    return true;
}

void OffstreamNHPartitioner::occupy_vertex(vid_t vid, vid_t d) {
    CHECK(!is_cores[current_partition].get(vid)) << "add " << vid << " to core again";
    // 核心集是vector<dense_bitset>，dense_bitset是一个稠密位图
    // 对位图的vid位置置1，表示vid被分配到current_partition分区
    is_cores[current_partition].set_bit_unsync(vid);
    // 如果顶点的度为0，不需要处理
    // d==0有两种可能，一是随机选择顶点，该顶点的度就是0，但是在前面选择的时候过滤掉这种顶点了；
    // 所以顶点为0只能是在S\C中选的顶点，但是该顶点的邻居已经被加入到S\C中，把它的度给更新了
    if (d == 0) return;

    // 把顶点x加入到边界集中，因为不知道当前顶点来自S还是C，走一遍逻辑
    // add_boundary的核心逻辑就是，把顶点加入到边界集，然后判断它的邻居是否在边界集，如果在就加边
    // add_boundary操作只会把当前定带你加入到边界集，不会引入新的顶点到边界集
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

size_t OffstreamNHPartitioner::check_edge(const edge_t *e) {
    rep (i, current_partition) {
        auto &is_boundary = is_boundaries[i];
        if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
            occupied[i] < capacity) {
            return i;
        }
    }

    rep (i, current_partition) {
        auto &is_core = is_cores[i], &is_boundary = is_boundaries[i];
        if ((is_core.get(e->first) || is_core.get(e->second)) &&
            occupied[i] < capacity) {
            if (is_core.get(e->first) && degrees[e->second] > avg_degree)
                continue;
            if (is_core.get(e->second) && degrees[e->first] > avg_degree)
                continue;
            is_boundary.set_bit(e->first);
            is_boundary.set_bit(e->second);
            return i;
        }
    }
    return num_partitions;
}

void OffstreamNHPartitioner::assign_edge(size_t partition, vid_t from, vid_t to) {
    // TODO 记录哪条边被分配了
    // save_edge(from, to, current_partition);
    true_vids.set_bit_unsync(from);
    true_vids.set_bit_unsync(to);

    is_mirrors[from].set_bit_unsync(partition);
    is_mirrors[to].set_bit_unsync(partition);

    assigned_edges++;
    occupied[partition]++;
    degrees[from]--;
    degrees[to]--;
}

void OffstreamNHPartitioner::add_boundary(vid_t vid) {
    // 获取到当前分区的核心集和边界集
    auto &is_core = is_cores[current_partition];
    auto &is_boundary = is_boundaries[current_partition];

    // TODO 什么情况顶点会在边界集中：它的邻居被加入到核心集，它被加入到边界集
    // 如果已经被加入到边界集，直接返回
    // TODO 为什么在边界集中就直接返回：因为如果它之前不在边界集中，在执行add_boundary操作时，已经把下面的步骤都执行完了
    if (is_boundary.get(vid)) {
        return;
    }
    // TODO 下面的操作是针对第一次把顶点加入到边界集
    // 这里不需要关心顶点是已经加入到核心集还是未加入到核心集，只要将顶点加入到边界集，都符合这个逻辑
    is_boundary.set_bit_unsync(vid);

    // 如果顶点没有在核心集中，直接把顶点的度数据加入到最小堆
    // 能够走到这一步，说明是将顶点加入到核心集中，将它的邻居加入到边界集，此时需要将它的度数计算出来加入到min heap
    if (!is_core.get(vid)) {
        min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
    }
    //仅支持无向图，在计算neighbor的时候有向和无向会导致邻居的差别从而影响分割
    // TODO 下面的步骤就是：把顶点x的邻居y = N(x)\S -> S
    // TODO 把z = N(y) x S -> e(y,z)加入到Ei

    // TODO 这里是加边的操作，不会把邻居加入到边界集
    rep (direction, 2) {
        adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
        // 遍历顶点的邻边
        for (size_t i = 0; i < neighbors.size();) {
            // 判断邻居edges[neighbors[i].v]是否已经分配
            if (edges[neighbors[i].v].valid()) {
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                if (is_core.get(u)) { // 如果顶点在核心集中
                    assign_edge(current_partition, direction ? vid : u,
                                direction ? u : vid);
                    vertex_partitions[u].insert(current_partition);
                    vertex_partitions[vid].insert(current_partition);
                    min_heap.decrease_key(vid); // 默认移除一条边
                    edges[neighbors[i].v].remove();
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else if (is_boundary.get(u) &&
                           occupied[current_partition] < capacity) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                    // 将边加入边集
                    assign_edge(current_partition, direction ? vid : u, direction ? u : vid);
                    vertex_partitions[u].insert(current_partition);
                    vertex_partitions[vid].insert(current_partition);
                    min_heap.decrease_key(vid);
                    min_heap.decrease_key(u);
                    edges[neighbors[i].v].remove();
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else
                    i++;
            } else { // TODO 对端顶点会走到这里
                //swap是pop的前提，先交换到最后位置然后把长度减1
                std::swap(neighbors[i], neighbors.back());
                neighbors.pop_back();
            }
        }
    }
}

size_t OffstreamNHPartitioner::leastLoad(set<size_t> partition_set) {
    // 遍历集合元素，找出最小occupied负载
    int min = INT_MAX;
    size_t partition_id;
    for (auto &partition: partition_set) {
        if (occupied[partition] < min) {
            min = occupied[partition];
            partition_id = partition;
        }
    }
    return partition_id;
}