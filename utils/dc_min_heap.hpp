#include "min_heap.hpp"
#include <queue>

using namespace  std;

class DcMinHeap {
private:
    int low_size;
    int high_size;
    int cap;
    double threshold;
    // 哈希结构
    vector<int> indices; // 顶点在堆中的下标
    vector<pair<int, vid_t>> rest_degrees;
public:
    DcMinHeap(int size, double avg_degree) {
        cap = size;
        low_size = 0;
        high_size = 0;
        threshold = avg_degree;
        reserve(size);
    }
    void reserve(int size) {
        cap = 0;
        indices.resize(size);
        rest_degrees.resize(size);
    }
    void clear() {
        cap = 0;
        low_size = 0;
        high_size = 0;
    }
    /**
     *
     * @param value
     * @param vid
     * @return
     */
    bool get_min(int& value, vid_t& vid, int& position) {
        if (cap > 0) {
            if (high_size == 0 || (low_size > 0 && rest_degrees[0].first < rest_degrees[cap - 1].first)) { // 只有左边有元素，两边都有元素，但是左边小
                value = rest_degrees[0].first;
                vid = rest_degrees[0].second;
                position = 0;
                return true;
            }
            int index = cap - 1;
            value = rest_degrees[index].first;
            vid = rest_degrees[index].second;
            position = index;
            return true;
        }
        return false;
    }
    int size() {
        return cap;
    }

    bool contains(vid_t vid) {
        return indices[vid] != - 1;
    }

    int low_shift_up(int index) {
        while(index > 0 && rest_degrees[index].first < rest_degrees[(index - 1) / 2].first) {
            int parent = (index - 1) / 2;
            swap(rest_degrees[index], rest_degrees[parent]);
            swap((indices[rest_degrees[index].second]), indices[rest_degrees[parent].second]);
            index = parent;
        }
        return index;
    }

    void low_shift_down(int index) {
        while (true) {
            int minIndex = index; // 假定当前节点是最小节点

            // 获取左孩子的索引
            int leftChildIndex = 2 * index + 1;
            if (leftChildIndex < low_size && rest_degrees[leftChildIndex].first < rest_degrees[minIndex].first) {
                minIndex = leftChildIndex; // 如果左孩子小于当前节点，则更新最小节点索引
            }

            // 获取右孩子的索引
            int rightChildIndex = 2 * index + 2;
            if (rightChildIndex < low_size && rest_degrees[rightChildIndex].first < rest_degrees[minIndex].first) {
                minIndex = rightChildIndex; // 如果右孩子小于当前节点，则更新最小节点索引
            }

            // 如果当前节点已经是两个子节点中的最小值，则无需继续调整
            if (minIndex == index) {
                break;
            } else {
                swap(rest_degrees[index], rest_degrees[minIndex]); // 交换当前节点与找到的最小节点
                index = minIndex; // 继续对新的节点执行下滤操作
            }
        }
    }


    int high_shift_up(int index) {
        while(index > 0 && rest_degrees[index].first < rest_degrees[(index - 1) / 2].first) {
            int parent = (index - 1) / 2;
            swap(rest_degrees[index], rest_degrees[parent]);
            swap((indices[rest_degrees[index].second]), indices[rest_degrees[parent].second]);
            index = parent;
        }
        return index;
    }



    void high_shift_down(int index) {

    }

    void insert(int value, vid_t vid) {
        auto pair = make_pair(value, vid);
        if (value < threshold) {
            rest_degrees[low_size] = pair;
            indices[vid] = low_size;
            int index = low_shift_up(low_size++);
            indices[vid] = index;
        } else {
            rest_degrees[cap - high_size - 1] = pair;
            indices[vid] = cap - high_size - 1;
            int index = high_shift_up(cap - high_size++);
            indices[vid] = index;
        }
    }

    bool remove(vid_t vid) {
        int index = indices[vid];
        if (index == -1) return false;
        cap--;
        if (index < low_size) { // 向上调整
            // 把最后一个元素前移
            rest_degrees[index] = rest_degrees[low_size-1];

            low_size--;

        } else {
            rest_degrees[index] = rest_degrees[cap - high_size];

            high_size--;
        }
        return true;
    }


};