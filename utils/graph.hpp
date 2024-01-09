//
// Created by muzongshen on 2021/9/23.
//
#pragma once

#include <vector>
// #include <parallel/algorithm>

#include "util.hpp"

//按照紧凑排列的方式，占用内存
struct uint40_t {
    //占用位段
    uint64_t v: 40; // 定义了64，但是只使用40，这个v存储的是边的id
} __attribute__((packed));

class adjlist_t {
private:
    uint40_t *adj; // 一个指向uint40_t类型数据的指针，存储顶点的邻边
    vid_t len;

public:
    adjlist_t() : adj(NULL), len(0) {}

    adjlist_t(uint40_t *adj, vid_t len = 0) : adj(adj), len(len) {}

    uint40_t *begin() { return adj; }

    uint40_t *end() { return adj + len; }

    void increment() { len++; }

    void push_back(size_t data) { adj[len++].v = data; }

    size_t size() const { return len; }

    uint40_t &operator[](size_t idx) { return adj[idx]; };

    const uint40_t &operator[](size_t idx) const { return adj[idx]; };

    uint40_t &back() { return adj[len - 1]; };

    const uint40_t &back() const { return adj[len - 1]; };

    void pop_back() { len--; }

    void clear() { len = 0; }
};

// 图结构
class graph_t {
private:
    vid_t num_vertices;
    size_t nedges;
    //neighbors是构造vdata的一个初始化数组
    uint40_t *neighbors; // 一个指向uint40_t类型数据的指针
    //TODO vdata使用adjlist_t来存储每个顶点的邻边
    //
    std::vector<adjlist_t> vdata;

public:
    graph_t() : num_vertices(0), nedges(0), neighbors(nullptr) {}

    // 析构函数
    ~graph_t() {
        if (neighbors) {
            // LOG(INFO) << neighbors << endl;
            free(neighbors);
        }
    }

    void resize(vid_t _num_vertices) {
        num_vertices = _num_vertices;
        vdata.resize(num_vertices);
    }

    [[nodiscard]] size_t num_edges() const { return nedges; }

    void build(const std::vector<edge_t> &edges);

    void build_reverse(const std::vector<edge_t> &edges);

    // 重载运算符
    adjlist_t &operator[](size_t idx) { return vdata[idx]; };

    const adjlist_t &operator[](size_t idx) const { return vdata[idx]; };
};
