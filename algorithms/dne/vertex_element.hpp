#include <iostream>
#include <queue>

struct Vertex_element {
    int key;
    int value;

    // 构造函数
    Vertex_element(int k, int v) : key(k), value(v) {}

    // 定义比较函数，用于 priority_queue 的排序
    bool operator<(const Vertex_element& other) const {
        return value < other.value;
    }
};