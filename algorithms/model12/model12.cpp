#include "model12.hpp"
#include <unordered_map>
using namespace std;

//固定随机数
// 构造函数
Model12Partitioner::Model12Partitioner(BaseGraph& baseGraph, const string &input, const string &algorithm,
                             size_t num_partitions)
        : EdgePartitioner(baseGraph, algorithm, num_partitions), input(input), gen(985) {
    config_output_files();
    current_partition = 0;
    average_degree = (double) num_edges * 2 / (double) num_vertices;
    assigned_edges = 0;
    capacity = num_edges * BALANCE_RATIO / num_partitions + 1;
    num_vertices_each_partition.assign(num_partitions, 0);

    is_cores.assign(num_partitions, dense_bitset(num_vertices));
    is_boundaries.assign(num_partitions, dense_bitset(num_vertices));
    true_vids.resize(num_vertices);
    master.assign(num_vertices, -1);
    dis.param(
            uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    degrees.resize(num_vertices);
    ifstream degree_file(degree_name(input), ios::binary);
    degree_file.read((char *) &degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();

    part_degrees.assign(num_vertices, vector<vid_t>(num_partitions));
    balance_vertex_distribute.resize(num_vertices);

    average_factor  = 1.5;
    front_factor = 1;
    front_partition = num_partitions * front_factor;

    N = 5;

}

//最后一个子图就是剩下边组合而成
void Model12Partitioner::assign_remaining() {
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

void Model12Partitioner::assign_master() {
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

size_t Model12Partitioner::count_mirrors() {
    size_t result = 0;
    rep(i, num_partitions) result += is_boundaries[i].popcount();
    return result;
}

void Model12Partitioner::split() {
    total_time.start();
    // 初始化最小堆，用于存储S\C的顶点信息
    min_heap.reserve(num_vertices);
    min_hd.reserve(num_vertices);
    d.reserve(num_vertices);
    repv(vid, num_vertices) {
        d.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
    }
    LOG(INFO) << "Start Model12 partitioning...";
    // 把参数写入文件，NE是边分割算法，计算复制因子
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "Model12" << endl
       << "BALANCE RATIO: " << BALANCE_RATIO
       << "| Average Factor: " << average_factor
       << "| Front Factor: " << front_factor
       << "，引入新的数据结构，实现在选择新的核心顶点的时候，考虑相同剩余度数顶点的邻居顶点带来的影响"
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    // 前p-1个分区
    LOG(INFO) << "front_partition: " << front_partition << endl;

    for (current_partition = 0; current_partition < num_partitions - 1; current_partition++) {
        // 当前分区的边数小于负载上限时，添加顶点到核心集C
        while (occupied[current_partition] < capacity) {
            vid_t degree, vid;
            // LOG(INFO) << "Min_heap " << current_partition << "  size: " << min_heap.size() << endl;
            // !min_heap.get_min(degree, vid) && !min_hd.get_min(degree, vid)
            if (degree_map.empty()) { // 当S\C为空时，从V\C中随机选择顶点
                // TODO 这里的逻辑不需要变
                if (!get_free_vertex(vid)) { // 当V\C已经没有顶点，结束算法
                    break;
                }
                // 计算顶点的出度和入度，该顶点之前没有被加入S\C，所以它的邻边必然没有被加入过Ei
                degree = adj_out[vid].size() + adj_in[vid].size();
            } else { // 当S\C不为空时，从S\C，即最小堆的堆顶移出顶点
                // TODO 这里是移除顶点的逻辑，需要更新
                auto degree_vids = degree_min_heap.begin();
                degree = degree_vids->first;
                set<vid_t>& vids = degree_vids->second;
                int size = vids.size();
                // LOG(INFO) << "Size: " << size << " in degree: " << degree << endl;
                // TODO 使用一些方法，避免重复计算
                if (size > 1) { // 执行选取顶点策略
                    // LOG(INFO) << "Execute Select Vid" << endl;
                    // 遍历set, 从n个顶点中选取最优
                    int n = N;
                    if (N > 0 && size < N) {
                      n = size;
                    }
                    // n = size; // 全量选取
                    // 1. 遍历顶点的一阶邻居，计算一阶邻居有多少个二阶邻居在候选集，选择候选集邻居最多的顶点
                    // vids是度数相同的顶点，对每个vid去找它不在边界集中的顶点
                    auto iterator = vids.begin();
                    int max = 0;
                    vid = *iterator;
                    for(int i = 0; i < n; i++){  // 这一步是遍历n个最小相同剩余度数的顶点
                        vid_t cur = *iterator;
                        int count = 0;
                        set<vid_t> same_cluster;
                        int same_cluster_size = 0;
                        // 对每个候选顶点，遍历他们的邻居
                        rep (direction, 2) {
                            adjlist_t &neighbors = direction ? adj_out[cur] : adj_in[cur];
                            // 遍历顶点的邻边
                            // LOG(INFO) << neighbors.size();
                            // TODO 判断彼此之间的边，这种是一个定值，可以一次遍历完成计算，不用除以2
                            for (auto & neighbor : neighbors) { // 对每个候选顶点的一阶邻居
                                // 判断邻居edges[neighbors[i].v]是否已经分配
                                if (edges[neighbor.v].valid()) {
                                    // 对这个顶点去找邻居，邻居在边界集合，则计数+1
                                    // TODO u是一阶邻居
                                    vid_t &u = direction ? edges[neighbor.v].second : edges[neighbor.v].first;
                                    same_cluster.insert(u);
                                    // 对这个V中的顶点去找在候选顶点中的顶点数量
                                    if (edge_pre_allocation.contains(u)) { // 存在，不需要计算
                                        // TODO 引入一个新的候选顶点到核心顶点，新增的边来自 候选顶点和一阶邻居 以及 一阶邻居和 候选集合/一阶邻居之间的边
                                        // 正常情况下，是遍历候选顶点的在V中的一阶顶点，判断一阶顶点的邻居有多少是在S中的，要注意彼此之间可能会存在边
                                        count = edge_pre_allocation[u].size();  // pair中的key-value是什么，key是顶点，value是边还是对端vid
                                    } else { // 不存在，计算，更新到edge_pre_allocation
                                        set<vid_t> neighbor_set = edge_pre_allocation[u];
                                        rep(u_direction, 2) {
                                            adjlist_t &u_neighbors = u_direction ? adj_out[u] : adj_in[u];
                                            for (auto & u_neighbor : u_neighbors) {
                                                if (edges[u_neighbor.v].valid()) {
                                                    vid_t &u_n = u_direction ? edges[u_neighbor.v].second : edges[u_neighbor.v].first;
                                                    if (same_cluster.contains(u_n)) {
                                                        same_cluster_size++;
                                                    }
                                                    neighbor_set.insert(u_n);
                                                    // 这些边都是可以预分配的，不然每次都要重复计算
                                                    // 可是要怎么去实现呢
                                                    // 通过一个map来维护？把已经预计算过的顶点都维护起来，感觉整个过程挺复杂
                                                    // 首先对所有相同最小剩余度数的候选顶点计算二阶邻居数；然后当确定把某个候选顶点加入到核心顶点之后，
                                                    // 更新一遍在计算集合中的候选顶点的二阶邻居数；那么这样的话，在找到最大边数的顶点时，只需要去更新map，
                                                    // 但是如果某个顶点不在map中，则不需要去更新，因为后续它加入是会被更新
                                                    // TODO 是要引入最多边的候选顶点
                                                    if (is_boundaries[current_partition].get(u_n)) {
                                                        count++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (count + same_cluster_size > max) {
                            vid = cur;
                            max = count + same_cluster_size;
                        }
                        iterator++;
                    }
                } else {
                    vid = *vids.begin();
                    // 删除
                    // LOG(INFO) << "Vids first to Remove: " << *degree_vids->second.begin() << endl;
                    // LOG(INFO) << "Remove vid from S: " << vid << ", degree: " << degree << endl;
                }
                // TODO vid为被选择加入核心集合，1. 将vid的邻居都从预分配中移除，加入到边集合，2. 判断vid的邻居未被引入的邻居，更新预分配集合；
                rep (direction, 2) {
                    adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                }
                degree_map.erase(vid);
                vids.erase(vid);
                // LOG(INFO) << "Vids first after Remove: " << *degree_vids->second.begin() << endl;
                if (vids.empty()) {
                    degree_min_heap.erase(degree);
                }

            }
            // 前面是获取顶点的两种方式，要么是从S\C，要么是从V\C
            // 把顶点加入到C，即核心集
            // 此时degree要么是从min_heap获取，要么是一个新顶点直接计算得到
            // LOG(INFO) << "vid: " << vid << ", degree: " << degree << endl;
            occupy_vertex(vid, degree);
        }
        degree_min_heap.clear();
        degree_map.clear();

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
                num_vertices_each_partition[j]++; // 每个分区的顶点数
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
//    repv(j, num_partitions) {
//        LOG(INFO) << "Partition " << j << " Vertex Count: " << current_partitions[j];
//    }

    ifstream fin(binary_edgelist_name(input), ios::binary | ios::ate);
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), ios::beg);
    edges.resize(num_edges);
    fin.read((char *) &edges[0], sizeof(edge_t) * num_edges);
    for (auto &e: edges) {
        vid_t sp = balance_vertex_distribute[e.first], tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
    }
    total_time.stop();
}

void Model12Partitioner::reindex() {
    // 随机选择顶点，进行广度遍历，重新索引
    vid_t index = 0;
    vid_t vid = dis(gen);
    // 基于该顶点进行深度遍历，对每个顶点重新索引
    v_queue.push(vid);
    while (!v_queue.empty()) {
        vid_t v = v_queue.front();
        v_queue.pop();
        // 将v加入到indices,重新索引
        indices.insert(std::pair<vid_t, vid_t>(index++, v));

        // 获取v的邻居顶点
        set < vid_t > neighbor_set = adjacency_list.find(v)->second;
        // 将顶点v的邻居加入到队列中，注意去重
        std::set<int> differenceSet;

        // 使用 std::set_difference 求差值
        std::set_difference(neighbor_set.begin(), neighbor_set.end(),
                            v_set.begin(), v_set.end(),
                            std::inserter(differenceSet, differenceSet.begin()));
        // 将neighbor_set加入v_queue和v_set中
        for (auto &i: differenceSet) {
            v_queue.push(i);
            v_set.insert(i);
        }
    }
}

bool Model12Partitioner::get_target_vertex(vid_t &vid) {
    // TODO 将随机选择顶点改成选择度最小的顶点，或者是距离当前分区所有节点距离最近的顶点
    // TODO 以上这个计算不太现实
    // 选择度最小的顶点，因为这样跨分区的边从一定概率来说是最小的
    if (d.size() == 0) return false;
    vid_t degree;
    d.get_min(degree, vid);
    d.remove(vid);
    return true;
}

bool Model12Partitioner::get_free_vertex(vid_t &vid) {
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

void Model12Partitioner::occupy_vertex(vid_t vid, vid_t d) {
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

size_t Model12Partitioner::check_edge(const edge_t *e) {
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

void Model12Partitioner::assign_edge(size_t partition, vid_t from, vid_t to) {
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

void Model12Partitioner::add_boundary(vid_t vid) {
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
        vid_t degree = adj_out[vid].size() + adj_in[vid].size();
        // LOG(INFO) << "Insert vid: " << vid << ", degree: " << degree << endl;
        degree_map.insert({vid, degree});
        if (!degree_min_heap.contains(degree)) {
            degree_min_heap.insert({degree, {}});
        }
        degree_min_heap[degree].insert(vid);
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
                    // 更新顶点的度数
                    if (degree_map.contains(vid)) {
                        vid_t old_degree = degree_map[vid];
                        auto& old_keys = degree_min_heap[old_degree];
                        vid_t new_degree = degree_map[vid] - 1;
                        if (new_degree == 0) {
                             degree_map.erase(vid);
                        } else {
                            degree_map[vid] = new_degree;
                            if (!degree_min_heap.contains(new_degree)) {
                                degree_min_heap.insert({new_degree, {}});
                            }
                            degree_min_heap[new_degree].insert(vid);
                        }
                        old_keys.erase(vid);
                        if (old_keys.empty()) {
                            degree_min_heap.erase(old_degree);
                        }
                    }
                    edges[neighbors[i].v].remove();
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                } else if (is_boundary.get(u) &&
                           occupied[current_partition] < capacity) { // 如果顶点在边界集中，并且当前分区负载没有达到上限
                    // 将边加入边集
                    assign_edge(current_partition, direction ? vid : u, direction ? u : vid);
                    // 更新顶点的度数
                    if (degree_map.contains(vid)) {
                        vid_t old_degree = degree_map[vid];
                        auto& old_keys = degree_min_heap[old_degree];
                        vid_t new_degree = degree_map[vid] - 1;
                        if (new_degree == 0) {
                            degree_map.erase(vid);
                        } else {
                            degree_map[vid] = new_degree;
                            if (!degree_min_heap.contains(new_degree)) {
                                degree_min_heap.insert({new_degree, {}});
                            }
                            degree_min_heap[new_degree].insert(vid);
                        }
                        old_keys.erase(vid);
                        if (old_keys.empty()) {
                            degree_min_heap.erase(old_degree);
                        }
                    }
                    // 更新顶点的度数
                    if (degree_map.contains(u)) {
                        vid_t old_degree = degree_map[u];
                        auto& old_keys = degree_min_heap[old_degree];
                        vid_t new_degree = degree_map[u] - 1;
                        if (new_degree == 0) {
                            degree_map.erase(u);
                        } else {
                            degree_map[u] = new_degree;
                            if (!degree_min_heap.contains(new_degree)) {
                                degree_min_heap.insert({new_degree, {}});
                            }
                            degree_min_heap[new_degree].insert(u);
                        }
                        old_keys.erase(u);
                        if (old_keys.empty()) {
                            degree_min_heap.erase(old_degree);
                        }
                    }

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