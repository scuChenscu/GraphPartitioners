#include <iostream>
#include <mutex>
#include <condition_variable>
#include <thread>

class Barrier {
public:
    //explicit Barrier(int count) : count_(count), current_count_(0), generation_(0) {}
    explicit Barrier() {}

    void init(int count) {
        count_ = count;
        current_count_ = 0;
        generation_ = 0;
    }


    void Wait() {
        std::unique_lock<std::mutex> lock(mutex_); // 这里有个锁

        int gen = generation_;

        ++current_count_;

        if (current_count_ < count_) {
            // 当前线程等待
            cv_.wait(lock, [this, gen] { return gen != generation_; });
        } else {
            // 最后一个线程到达，重置计数器，唤醒所有等待的线程
            current_count_ = 0;
            ++generation_; // 迭代轮次更新
            cv_.notify_all();
        }
    }

    void NotifyExit() {
        std::unique_lock<std::mutex> lock(mutex_);
        // 设置退出标志，唤醒所有等待的线程，避免死锁
        // 检查当前线程数，如果除了当前线程外的其他线程都已经到达，则唤醒所有线程
        count_ = count_ - 1;
        // 如果当前线程是最后一个到达
        if (current_count_ == count_) {
            current_count_ = 0;
            ++generation_;
            cv_.notify_all();
        }
    }

private:
    int count_;
    int current_count_;
    int generation_;
    std::mutex mutex_;
    std::condition_variable cv_;
};