//num_partition "number of partitions"
//memory_size "memory size in MB. memory_size is used as a chunk size if shuffling and batch size while partitioning"
//algorithm "partition algorithm"
//dataset "name of the dataset to construct graph."
//filename "name of the file to store edge list of a graph."
//shuffle "for streaming graph, influencing the partitioning quality"

#include <iostream>
#include "converter/conversions.hpp"
#include "utils/util.hpp"
#include "algorithms/ne/ne.hpp"
#include "algorithms/dbh/dbh.hpp"
#include "algorithms/hdrf/hdrf.hpp"
#include "algorithms/ldg/ldg.hpp"
#include "algorithms/ours/ours.hpp"
#include "algorithms/fennel/fennel.hpp"
#include "algorithms/model1/model1.hpp"
#include "algorithms/model2/model2.hpp"
#include "algorithms/model3/model3.hpp"
#include "algorithms/model4/model4.hpp"
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

int main() {
    // int num_partition = 5;
    // int partitions[] = {2, 4, 8, 16, 32, 64};
    int partitions[] = { 8, 32};
    // int partitions[] = {64};
    int memory_size = 4096;
    double lambda = 1.1;
    double balance_ratio = 1.05;
    // string algorithms[] = {"ne", "dbh", "hdrf", "ldg", "fennel"};
    string algorithms[] = { "model4"};
    // 输入文件夹，存.graph文件，文件首行为顶点数 边数；其他行为邻接表
    string input = "../graphs/input";
    // TODO 需要一个文件，追加输出运行结果
    // TODO 计算负载均衡、边割率、复制因子，输出到index文件
    // 输出文件夹，每次输出以 timestamp/graphs/partitions/algorithms/ 命名
//    string output = "../graphs/output";
//    string graphname = "copter2.graph";
//    string filename = "../graphs/input/" + graphname;
//    string original ="../graphs/" + graphname;
    // TODO 随机打乱
    bool shuffle = false;

    google::InitGoogleLogging("main");  //参数为自己的可执行文件名
    FLAGS_logtostderr = true;

    // read the edges and change to binary format file
    Converter *converter;
    string binary_edgelist; // binary edgelist file
    LOG(INFO) << "Using normal dataset, dont shuffle";
    // 遍历input下的.graph文件
    string current_time = getCurrentTime();
    stringstream ss;
    ss << "==============================================================================================================================================" << endl
       << "Time:" << current_time
       << " lambda:" << lambda
       << " Balance ratio:" << balance_ratio
       << " Memory size:" << memory_size
       << " Shuffle:" << shuffle
       << endl;
    LOG(INFO) << ss.str();
    appendToFile(ss.str());
    ss.clear();
    for (const auto &entry: fs::directory_iterator(input)) {
        // 判断entry是否为文件
        if (fs::is_regular_file(entry)) {
            string filename = entry.path().string();
            if (!filename.ends_with("mdual.graph")) continue;
            LOG(INFO) << "Convert " << filename << " to binary edgelist" << endl;
            converter = new Converter(filename);
            convert(filename, converter, memory_size);
            LOG(INFO) << "Execute algorithms on: " << filename << endl;
            stringstream graph;
            graph << "Graph Name: " << filename << endl;
            appendToFile(graph.str());
            int size = sizeof(algorithms) / sizeof(algorithms[0]);
            int partition = sizeof(partitions) / sizeof(partitions[0]);
            for (int j = 0; j < partition; j++) {
                int num_partition = partitions[j];
                stringstream part;
                part << "Partition: " << num_partition << endl;
                LOG(INFO) << part.str();
                appendToFile(part.str());
                for (int i = 0; i < size; i++) {
                    string algorithm = algorithms[i];
                    LOG(INFO) << "execute " << algorithm << " on: " << filename << endl;
                    Partitioner *partitioner;
                    if (algorithm == "ne")
                        partitioner = new NePartitioner(filename, algorithm, num_partition);
                    else if (algorithm == "ours")
                        partitioner = new OursPartitioner(filename, algorithm, num_partition);
                    else if (algorithm == "model1")
                        partitioner = new Model1Partitioner(filename, algorithm, num_partition);
                    else if (algorithm == "model2")
                        partitioner = new Model2Partitioner(filename, algorithm, num_partition);
                    else if (algorithm == "model3")
                        partitioner = new Model3Partitioner(filename, algorithm, num_partition);
                    else if (algorithm == "model4")
                        partitioner = new Model4Partitioner(filename, algorithm, num_partition);
                    else if (algorithm == "dbh")
                        partitioner = new DbhPartitioner(filename, algorithm, num_partition, memory_size, shuffle);
                    else if (algorithm == "hdrf")
                        partitioner = new HdrfPartitioner(filename, algorithm, num_partition, memory_size,
                                                          balance_ratio,
                                                          lambda, shuffle);
                    else if (algorithm == "ldg")
                        partitioner = new LdgPartitioner(filename, algorithm, num_partition, memory_size, shuffle);
                    else if (algorithm == "fennel")
                        partitioner = new FennelPartitioner(filename, algorithm, num_partition, memory_size, shuffle);
                    else {
                        LOG(ERROR) << "Unknown algorithm: " << algorithm;
                        continue;
                    }
                    partitioner->split();
                    delete partitioner;
                }
            }
        }
    }
    google::ShutdownGoogleLogging();
    return 0;
}
