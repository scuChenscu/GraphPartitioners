#pragma once

#include <utility>
#include <chrono>
#include <cstdint>
#include <sys/stat.h>
#include <glog/logging.h>
#include <sstream>
#include <fstream>

using namespace std;

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define repv(i, n) for (vid_t i = 0; i < n; ++i)

typedef uint32_t vid_t;
const vid_t INVALID_VID = -1;

struct edge_t {
    vid_t first, second;

    edge_t() : first(0), second(0) {}

    edge_t(vid_t first, vid_t second) : first(first), second(second) {}

    [[nodiscard]] bool valid() const { return first != INVALID_VID; }

    void remove() { first = INVALID_VID; }
};

inline string change2tmpdir(const string &str) {
    string dir_path, file_name;
    size_t found = str.find_last_of('/');
    dir_path = str.substr(0, found);
    file_name = str.substr(found + 1);
    // dir_path += "/tmp_partitioning_dir";
    dir_path = "../graphs/output/";
    string cmd = "mkdir -p " + dir_path;
    int flag = system(cmd.c_str());
    if (flag == 1) {
        LOG(FATAL) << "make tmp dir error!";
    }
    return dir_path + "/" + file_name;
}

inline void appendToFile(const std::string& content) {
    string filename = "../graphs/index";
    std::ofstream file(filename, std::ios::app); // 打开文件以追加模式写入

    if (file.is_open()) {
        file << content; // 追加文件内容
        file.close(); // 关闭文件
    } else {
        LOG(INFO) << "Failed to append content to file." << endl;
    }
}

inline string getFormattedTime() {
    std::time_t result = std::time(0);
    std::tm* current_time = std::localtime(&result);
    char buffer[11];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d%H%M%S", current_time);
    return buffer;
}

inline std::string getCurrentTime() {
    std::time_t result = std::time(nullptr);
    std::tm* localTime = std::localtime(&result);

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << localTime->tm_year + 1900;
    ss << "/";
    ss << std::setw(2) << std::setfill('0') << localTime->tm_mon + 1;
    ss << "/";
    ss << std::setw(2) << std::setfill('0') << localTime->tm_mday;

    ss << " ";
    ss << std::setw(2) << std::setfill('0') << localTime->tm_hour;
    ss << ":";
    ss << std::setw(2) << std::setfill('0') << localTime->tm_min;
    ss << ":";
    ss << std::setw(2) << std::setfill('0') << localTime->tm_sec;

    return ss.str();
}



void writea(int f, char *buf, size_t nbytes);

inline string binary_edgelist_name(const string &filename) {
    // string new_filename = change2tmpdir(filename);
    string path = "../graphs/binedgelist/";
    size_t found = filename.find_last_of('/');
    string new_filename  = path + filename.substr(found + 1);
    stringstream ss;
    ss << new_filename << ".binedgelist";
    return ss.str();
}

inline string shuffled_binary_edgelist_name(const string &input) {
    string filename = change2tmpdir(input);
    stringstream ss;
    ss << filename << ".shuffled.binedgelist";
    return ss.str();
}

inline string degree_name(const string &input) {
    string filename = change2tmpdir(input);
    stringstream ss;
    ss << filename << ".degree";
    return ss.str();
}

inline string edge_partition_filename(const string &graph_name, const string &algorithm, int num_partition) {
    string filename = change2tmpdir(graph_name);
    stringstream ss;
    ss << filename << "." << algorithm << ".edge." << num_partition;
    return ss.str();
}

inline string vertex_partition_filename(const string &graph_name, const string &algorithm, int num_partition) {
    string filename = change2tmpdir(graph_name);
    stringstream ss;
    ss << filename << "." << algorithm << ".vertex." << num_partition;
    return ss.str();
}

inline bool is_exists(const string &name) {
    struct stat buffer{};
    return (stat(name.c_str(), &buffer) == 0);
}

class Timer {
private:
    chrono::system_clock::time_point t1, t2;
    double total;

public:
    Timer() : total(0) {}

    void reset() { total = 0; }

    void start() { t1 = chrono::system_clock::now(); }

    void stop() {
        t2 = chrono::system_clock::now();
        chrono::duration<double> diff = t2 - t1;
        total += diff.count();
    }

    [[nodiscard]] double get_time() const { return total; }
};
