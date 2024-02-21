

#pragma once

#include <string>
#include <vector>
#include <bitset>
#include <unordered_set>
#include "../../utils/dense_bitset.hpp"
#include "../../utils/util.hpp"
#include "../../partitioner/edgePartitioner.hpp"

using namespace std;

class AdwisePartitioner : public EdgePartitioner {
private:
    ifstream fin;
    vid_t num_vertices{};
    size_t num_edges{};
    size_t filesize;

    // batch processing globals
    uint32_t num_batches;
    uint32_t num_edges_per_batch;
    double lambda;

    int max_degree;

    vector<size_t> degrees;

    vector<size_t> partial_degrees;

    uint64_t max_partition_load;

    //balance vertex partition distribute
    vector<vector<vid_t> > part_degrees; //each partition node degree
    vector<int> balance_vertex_distribute; //each node belongs to which unique partition

    vector<uint64_t> edge_load;
    // vector<dense_bitset> is_mirrors;
    dense_bitset true_vids;
    uint64_t min_load = UINT64_MAX;
    uint64_t max_load = 0;
    const double epsilon = 1;

    int index = 0;
    unordered_set<int> C;
    unordered_set<int> Q;
    unordered_set<int> W;
    int w = 1;

    double threshold = 0.0;




protected:
    void batch_adwise(vector<edge_t> &edges);

    int find_max_score_partition(edge_t &e);

    void update_is_mirrors(edge_t &e, int max_p);

    void update_min_max_load(int max_p);

    void batch_node_assign_neighbors(vector<edge_t> &edges);

    void read_and_do(const string& opt_name);

    void calculate_set();
    void get_best_assignment(int &e_id, int &p_id);
    double calculate_score(int e_id);

public:
    AdwisePartitioner(BaseGraph& baseGraph, const string& input, const string& algorithm,
                    size_t num_partitions,
                    int memory_size,
                    double balance_ratio,
                    double balance_lambda,
                    bool shuffle);
    void split();

};

