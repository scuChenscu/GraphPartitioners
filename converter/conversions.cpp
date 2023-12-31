#include <cstring>

#include "conversions.hpp"

using namespace std;

// Removes \n from the end of line
void FIXLINE(char *s) {
    int len = (int) strlen(s) - 1;
    if (s[len] == '\n')
        s[len] = 0;
}

void convert_edgelist(const string& filename, Converter *converter) {
    FILE *inf = fopen(filename.c_str(), "r");
    size_t bytesread = 0;
    size_t linenum = 0;
    if (inf == nullptr) {
        LOG(FATAL) << "Could not load:" << filename
                   << ", error: " << strerror(errno) << std::endl;
    }

    LOG(INFO) << "Reading in edge list format" << std::endl;
    char s[1024];
    while (fgets(s, 1024, inf) != nullptr) {
        linenum++;
        if (linenum % 10000000 == 0) {
            LOG(INFO) << "Read " << linenum << " lines, "
                      << bytesread / 1024 / 1024. << " MB" << std::endl;
        }
        FIXLINE(s);
        bytesread += strlen(s);
        if (s[0] == '#')
            continue; // Comment
        if (s[0] == '%')
            continue; // Comment

        char delims[] = "\t, ";
        char *t;
        t = strtok(s, delims);
        if (t == nullptr) {
            LOG(FATAL) << "Input file is not in right format. "
                       << "Expecting \"<from>\t<to>\". "
                       << "Current line: \"" << s << "\"\n";
        }
        vid_t from = atoi(t);
        t = strtok(nullptr, delims);
        if (t == nullptr) {
            LOG(FATAL) << "Input file is not in right format. "
                       << "Expecting \"<from>\t<to>\". "
                       << "Current line: \"" << s << "\"\n";
        }
        vid_t to = atoi(t);

        if (from != to) {
//            LOG(INFO) <<from<<" "<<to<<std::endl;
            converter->add_edge(from, to);
        }
    }
    fclose(inf);
}

void convert(const std::string& filename, Converter *converter, int memory_size = 4096) {
    if (filename.empty()) {
        LOG(FATAL) << "Empty filename";
    }
    if (converter->done()) {
        LOG(INFO) << "Skip convert, .binedgelist file exists";
        return;
    }
    converter->init(memory_size);
    convert_edgelist(filename, converter);
    converter->finalize();
}

// 将邻接表转为from to形式
void convert_adjacency_list(string graph_filename, string output) {
    LOG(INFO) << "converting adjacency list to [from to]" << graph_filename;
    ifstream inFile(graph_filename);
    if (!inFile) {
        LOG(FATAL) << "empty file name";
    }
    vector<pair<int, int>> edges;
    int v_size, e_size;
    string line;
    getline(inFile, line); // 跳过第一行
    istringstream _iss(line);
    _iss >> v_size >> e_size;
    LOG(INFO) << v_size << " " << e_size << endl;
    int num = 0;
    while (getline(inFile, line)) { // 读取一行
        int from;
        istringstream iss(line); // 获取当前行
        iss >> from; // 获取第一个数字
        LOG(INFO) << from << endl;
        int to;
        while (iss >> to) { // 获取剩下的数字
            num++;
            // 打印rest
            LOG(INFO) << from << " " << to << endl;
            std::pair<int, int> pair(from, to); // 将第一个数字与剩下的数字组成pair
            edges.push_back(pair); // 将pair添加到pairs容器中
        }
    }
    LOG(INFO) << "num: " << num << endl;
    sort(edges.begin(), edges.end());

    auto last = std::unique(edges.begin(), edges.end());
    edges.erase(last, edges.end());
    inFile.close();


    ofstream outFile(output);
    if (!outFile) {
        LOG(FATAL) << "Error: Cannot open file for writing!" << endl;
        return;
    }

    for (const auto &edge: edges) {
        outFile << edge.first << " " << edge.second << endl;
    }

    outFile.close();

    LOG(INFO) << "Conversion completed successfully." << endl;
}

