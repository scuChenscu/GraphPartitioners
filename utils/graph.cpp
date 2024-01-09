#include "graph.hpp"


// build方法，传入边集合
void graph_t::build(const std::vector<edge_t> &edges) {
    LOG(INFO) << "building graph" << endl;
    // 构造过程
    //可能存在孤立节点，在有边的时候才创建
    // LOG(INFO) << nedges << endl;
    if (edges.size() > nedges) {
        // 重新分配空间，neighbors存储所有边的数据，通过vdata[vid]来获取到对应vid的邻边
        neighbors = (uint40_t *) realloc(neighbors, sizeof(uint40_t) * edges.size());
    }
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();
    // 初始化操作
    vector<size_t> count(num_vertices, 0);
    // count计算每个顶点有多少个邻居
    for (size_t i = 0; i < nedges; i++) {
        count[edges[i].first]++;

    }
    // uint40_t *neighbors;
    vdata[0] = adjlist_t(neighbors);
    // 通过前缀和创建出紧凑的邻接表
    for (vid_t v = 1; v < num_vertices; v++) {
        // count是一个数组，记录每个顶点有多少个邻居，这里转换成前缀和数组
        count[v] += count[v - 1];
        // 这里neighbors + count[v-1]是地址偏移，&neighbors是取值，没有&neighbors是取地址
        // vdata存储的是vector<adjlist_t>，adjlist_t包含一个指向uint40_t的指针和vid_t类型的len变量
        // 所以这里是对每个顶点的邻居进行初始化，初始化的时候是存储邻居的地址偏移量，也就是uint40_t的指针
        vdata[v] = adjlist_t(neighbors + count[v - 1]);
    }
    // 遍历边
    for (size_t i = 0; i < edges.size(); i++) {
        // void push_back(size_t data) { adj[len++].v = data; }
        // 这里是将边的索引值存储到邻接表中
        // 传入边的id？为啥是存边，不是存点？
        // 所以vdata存储的是每个顶点的邻边数据，用于在ne算法允许时，能够快速获取到邻边，者就是build方法的作用？
        vdata[edges[i].first].push_back(i);
    }
    // build方法先分配空间，在对每个顶点v，记录它的邻边id，邻边数目len，初始化完成vdata
}

void graph_t::build_reverse(const std::vector<edge_t> &edges) {
    LOG(INFO) << "building reverse graph " << endl;
    if (edges.size() > nedges)
        neighbors = (uint40_t *) realloc(neighbors, sizeof(uint40_t) * edges.size());
    CHECK(neighbors) << "allocation failed";
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++)
        count[edges[i].second]++;

    vdata[0] = adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v - 1];
        vdata[v] = adjlist_t(neighbors + count[v - 1]);
    }
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].second].push_back(i);
}
