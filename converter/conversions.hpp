#pragma once

#include <string>
#include <fstream>

#include <boost/unordered_map.hpp>
#include <utility>

#include "../utils/util.hpp"

using namespace std;

DECLARE_string(graphtype);

class Converter {
protected:
    string filename;
    vid_t num_vertices;
    size_t num_edges;
    std::vector<vid_t> degrees;
    std::ofstream fout;
    boost::unordered_map<vid_t, vid_t> name2vid;
    vid_t max_vid = 0;

    vid_t get_vid(vid_t v) {
        auto it = name2vid.find(v);
        if (it == name2vid.end()) {
            name2vid[v] = num_vertices;
            degrees.resize(num_vertices + 1);
            return num_vertices++;
        }
        return name2vid[v];
    }

public:
    Converter(string filename) : filename(std::move(filename)) {}

    virtual ~Converter() = default;

    virtual bool done() { return is_exists(binary_edgelist_name(filename)); }

    virtual void init(int memory_size) {
        num_vertices = 0;
        num_edges = 0;
        degrees.reserve(1 << 30);
        fout.open(binary_edgelist_name(filename), ios::binary);
        fout.write((char *) &num_vertices, sizeof(num_vertices));
        fout.write((char *) &num_edges, sizeof(num_edges));
    }

    virtual void add_edge(vid_t from, vid_t to) {
        if (to == from) {
            LOG(WARNING) << "Tried to add self-edge " << from << "->" << to
                         << std::endl;
            return;
        }

        if (from > to) {
            if (from > max_vid) {
                max_vid = from;
            }
        } else {
            if (to > max_vid) {
                max_vid = to;
            }
        }

        num_edges++;
//        from = get_vid(from);
//        to = get_vid(to);
//        get_vid(from);
//        get_vid(to);
        degrees[from]++;
        degrees[to]++;

        fout.write((char *) &from, sizeof(vid_t));
        fout.write((char *) &to, sizeof(vid_t));
    }

    virtual void finalize() {
        // TODO
        num_vertices = max_vid + 1;
        fout.seekp(0);
        LOG(INFO) << "num_vertices: " << num_vertices << ", num_edges: " << num_edges;
        fout.write((char *) &num_vertices, sizeof(num_vertices));
        fout.write((char *) &num_edges, sizeof(num_edges));
        fout.close();

        fout.open(degree_name(filename), std::ios::binary);
        fout.write((char *) &degrees[0], num_vertices * sizeof(vid_t));
        fout.close();
    }
};

void convert(const string &filename, Converter *converter, int memory_size);

void convert_adjacency_list(std::string graph_filename, std::string output);