//num_partition "number of partitions"
//memory_size "memory size in MB. memory_size is used as a chunk size if shuffling and batch size while partitioning"
//algorithm "partition algorithm"
//dataset "name of the dataset to construct graph."
//filename "name of the file to store edge list of a graph."
//shuffle "for streaming graph, influencing the partitioning quality"
#pragma once

#include <filesystem>
#include <iostream>
#include "converter/conversions.hpp"
#include "utils/util.hpp"
#include "baseGraph/base_graph.hpp"

using namespace std;
namespace fs = filesystem;
// 原始图数据的路径
const static string input = "../graphs/input";

void signalHandler(int signum) {
    std::cerr << "Segmentation Fault (signal " << signum << ")" << std::endl;

    // 打印调用栈信息（需要编译时开启调试信息）
    // 注意: 打印调用栈信息可能需要开启编译器选项，如 "-g"。
    // 如果使用 g++，可以使用 "-g" 选项编译代码。
    std::cerr << "Call stack:\n";
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
    Converter *converter;
    // LOG(INFO) << "Using normal dataset, dont shuffle";
    // 遍历input下的.graph文件
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "==============================================================================================================================================" << endl
       << "Time: " << current_time
       << " | Input: " << input
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());
    ss.clear();

    for (const auto &entry: fs::directory_iterator(input)) {
        // 判断entry是否为文件
        if (fs::is_regular_file(entry)) {
            string graph_name = entry.path().string();
            if (!graph_name.ends_with(graph_suffix)) continue;
            LOG(INFO) << "Convert " << graph_name << " to binary edgelist" << endl;
            converter = new Converter(graph_name);
            convert(graph_name, converter, memory_size);
            delete converter;
            auto *baseGraph = new BaseGraph(graph_name);
            baseGraph->partition();
            LOG(INFO) << "Finish partition on " << graph_name << endl;
            delete baseGraph;
        }
    }
    google::ShutdownGoogleLogging();
    return 0;
}

