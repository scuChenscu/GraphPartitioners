#pragma once

#include <vector>
#include <unordered_map>
#include "util.hpp"

using namespace std;
// 模板函数，有点类似Java中的泛型
template<typename ValueType, typename KeyType, typename IdxType = vid_t>
class MinHeap {
public:
    IdxType n;
    std::vector<std::pair<ValueType, KeyType> > heap;
    // 记录每个顶点vid在heap中对应的下标
    std::vector<IdxType> key2idx;
// 构造方法，初始化成员变量n为0，heap为空，key2idx为空。
public:
    MinHeap() : n(0), heap(), key2idx() {}
    MinHeap(IdxType nelements) {
        reserve(nelements);
    }

    // 这是一个调整堆的过程
    IdxType shift_up(IdxType cur) {
        if (cur == 0) return 0;
        IdxType p = (cur - 1) / 2;
        // value是顶点的度，key是顶点vid
        if (heap[cur].first < heap[p].first) {
            std::swap(heap[cur], heap[p]);
            std::swap(key2idx[heap[cur].second], key2idx[heap[p].second]);
            return shift_up(p);
        }
        return cur;
    }

    void shift_down(IdxType cur) {
        IdxType l = cur * 2 + 1;
        IdxType r = cur * 2 + 2;

        if (l >= n)
            return;

        IdxType m = cur;
        if (heap[l].first < heap[cur].first)
            m = l;
        if (r < n && heap[r].first < heap[m].first)
            m = r;

        if (m != cur) {
            std::swap(heap[cur], heap[m]);
            std::swap(key2idx[heap[cur].second], key2idx[heap[m].second]);
            shift_down(m);
        }
    }

    void insert(ValueType value, KeyType key) {
        heap[n] = std::make_pair(value, key);
        key2idx[key] = n++;
        IdxType cur = shift_up(n - 1);
        shift_down(cur);
    }

    int size() {
        return n;
    }

    bool contains(KeyType key) {
        return key2idx[key] < n && heap[key2idx[key]].second == key;
    }

    void decrease_key(KeyType key, ValueType d = 1) {
        if (d == 0) return;
        IdxType cur = key2idx[key];
        CHECK(cur < n && heap[cur].second == key) << "key not found";
        CHECK_GE(heap[cur].first, d) << "value cannot be negative";
        heap[cur].first -= d;
        shift_up(cur);
    }

    bool remove(KeyType key) {
        IdxType cur = key2idx[key];
        if (cur >= n || heap[cur].second != key)
            return false;

        n--;
        if (n > 0) {
            heap[cur] = heap[n];
            key2idx[heap[cur].second] = cur;
            cur = shift_up(cur);
            shift_down(cur);
        }
        return true;
    }

    bool get_min(ValueType &value, KeyType &key) {
        // 如果堆里面有元素，选择堆顶元素，因为堆顶元素引入的新顶点最少
        if (n > 0) {
            value = heap[0].first; // 顶点的度数
            key = heap[0].second; // 顶点的vid
            return true;
        } else
            // 堆为空时，返回false，随机选择顶点
            return false;
    }
    // TODO 递归爆栈
    bool get_pure_min(ValueType &value, KeyType &key, dense_bitset* dirty_vertices, vector<unordered_map<size_t, edge_t*>>* vertex_adjacent_edges, size_t index) {
        // 如果堆里面有元素，选择堆顶元素，因为堆顶元素引入的新顶点最少
        // LOG(INFO) << index << std::endl;

        if (n > 0) {
            value = heap[0].first; // 顶点的度数
            key = heap[0].second; // 顶点的vid
            if (!dirty_vertices->empty() && index > 0) {
                if (dirty_vertices->get(key)) {
                    value = vertex_adjacent_edges[key].size(); // 顶点的度数
                    if (value > 0) return true;
                    get_pure_min(value, key, dirty_vertices, vertex_adjacent_edges, index);
                }
            }
            return true;
        } else
            // 堆为空时，返回false，随机选择顶点
            return false;
    }

    void reserve(IdxType nelements) {
        n = 0;
        heap.resize(nelements);
        key2idx.resize(nelements);
    }

    void clear() {
        n = 0;
    }
};