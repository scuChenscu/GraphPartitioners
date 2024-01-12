#pragma once

#include "../utils/util.hpp"


class Vertex_t {
private:
    unordered_map<size_t , edge_t*>& adjacent_edges;
    vid_t vid;
public:
    Vertex_t(vid_t& vid, unordered_map<size_t , edge_t*>& adjacent_edges) : vid(vid), adjacent_edges(adjacent_edges) {
        this->adjacent_edges = adjacent_edges;
        this->vid = vid;
    }
    void remove_edge(size_t e_id) {
        adjacent_edges.erase(e_id);
    }

    size_t get_edge_count() {
        return adjacent_edges.size();
    }
};