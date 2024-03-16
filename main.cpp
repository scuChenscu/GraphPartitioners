//num_partition "number of partitions"
//memory_size "memory size in MB. memory_size is used as a chunk size if shuffling and batch size while partitioning"
//algorithm "partition algorithm"
//dataset "name of the dataset to construct graph."
//filename "name of the file to store edge list of a graph."
//shuffle "for streaming graph, influencing the partitioning quality"
#pragma once
#include <iostream>
#include <csignal>
#include "converter/conversions.hpp"
#include "converter/shuffler.hpp"
#include "utils/util.hpp"
#include "baseGraph/base_graph.hpp"
#include <filesystem>
namespace fs = filesystem;
using namespace std;
// #include <experimental/filesystem>
// namespace fs = experimental::filesystem;
// 原始图数据的路径

void signalHandler(int signum) {
    LOG(ERROR) << "Segmentation Fault (signal " << signum << ")" << endl;
    // 打印调用栈信息（需要编译时开启调试信息）
    // 注意: 打印调用栈信息可能需要开启编译器选项，如 "-g"。
    // 如果使用 g++，可以使用 "-g" 选项编译代码。
    LOG(ERROR) << "Call stack:" << endl;
    system("backtrace");
    // 退出程序
    exit(signum);
}

int main() {
    signal(SIGABRT, signalHandler);
    signal(SIGSEGV, signalHandler);
    google::InitGoogleLogging("main");  //参数为自己的可执行文件名
    FLAGS_logtostderr = true;

    // read the edges and change to binary format file
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "==============================================================================================================================================" << endl
       << "Time: " << current_time
       << " | Input: " << input
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());

    for (const auto &entry: fs::directory_iterator(input)) {
        // 判断entry是否为文件
        if (fs::is_regular_file(entry)) {
            string graph_name = entry.path().string();
            string name = graph_name.substr(0, graph_name.size() - 6);
            // 判断是否为.graph类型文件
            int len = graph_suffix.length();
            string substring = graph_name.substr(graph_name.size() - len, len);
            if (substring != graph_suffix) continue;
            LOG(INFO) << "Convert " << graph_name << " to binary edgelist" << endl;
            Converter *converter;
            if (need_to_shuffle) {
                LOG(INFO) << "Using shuffle dataset" << endl;
                converter = new Shuffler(graph_name);
            } else {
                converter = new Converter(graph_name);

            }
            convert(graph_name, converter, memory_size);
            delete converter;

            auto *baseGraph = new BaseGraph(graph_name);
//            baseGraph->partition();
//            LOG(INFO) << "Finish partition on " << graph_name << endl;
            delete baseGraph;
        }
    }
    google::ShutdownGoogleLogging();
    return 0;
}

