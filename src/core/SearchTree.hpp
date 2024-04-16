#include <cstddef>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <array>
#include <vector>
#include <atomic>
#include <random>

#define DIM 128

struct Hyperplane
{

    float offset;
    float *normal_vec = new float[DIM];
};

class SearchTree
{
public:
    SearchTree(int n_points, size_t dimension, const float **full_data, int n_threads, int K);

    ~SearchTree();

    void construct_tree();

    void write_results();
    std::vector<float *> search_results(float *data_point, int k, unsigned int region_id);

    // private:
    // data
    int n_points;        // Number of all points in the dataset
    size_t D;            // dimension of points
    const int n_threads; // Number of threads available
    const size_t K;         // Maximum number of points in a terminal region
    size_t curr_length;  // Current length of array holding non-terminal points
    
    // random number generation
    std::random_device rd;
    std::mt19937 gen;

    std::unordered_map<unsigned int, Hyperplane *> region_hyperplane_map; // Map from region ID to hyperplane
    std::vector<unsigned int> active_regions;                             // TODO: add description
    std::unordered_set<unsigned int> terminal_regions;                    // region IDs of all terminal regions
    std::unordered_set<unsigned int> new_terminal_regions;                // region IDs of terminal regions found in a single iteration
    std::pair<float *, float *> *sampled_points;                          // pairs of points sample from active regions

    bool all_points_terminal = false;     // Flag indicating whether any non-terminal points remain
    unsigned int n_active_points_;        // Number of active points
    const float **full_data;              // Grid containing the data to train on
    float **data;                         // Array of pointers to each data point
    unsigned int *regions_curr;           // Array of current region assignments
    unsigned int *regions_next;           // Array of next region assignments
    unsigned int **regions_counts;        // Grid of number of points in each region (per-thread)
    unsigned int *regions_reduced_counts; // Array of number of points in each region (total)
    std::unordered_map<unsigned int, std::vector<float *>> terminal_region_data_store;

    // Hyperplane *hyperplanes; // Array of hyperplanes
    std::vector<Hyperplane> hyperplanes; // Array of hyperplanes


    // methods for building tree
    void classify(int start_idx, int end_idx_exclusive, int tid);
    void handle_load_imbalance(int *reduced_region_counts);
    void calculate_load();
    void sample_points(const std::vector<unsigned int> &active_regions, const unsigned int *regions_reduced_counts, std::pair<float *, float *> *sampled_points);
    inline float signed_distance(Hyperplane &hp, const float *point);
    inline void inc_region_count_(const int region, const int tid);
    inline int get_region_count(const int tid, const int region);
    void reduce_counts(unsigned int *reduced_region_counts);
    void make_hyperplane_(std::pair<float *, float *> &sampled_point, Hyperplane &hp);
    void update_active_regions(std::vector<unsigned int> &active_regions, unsigned int *reduced_region_counts);
    void move_terminal_points_(unsigned int start_idx, unsigned int end_idx_exclusive,
                               std::unordered_set<unsigned int> &terminal_regions);
    float l2_distance(float *v1, float *v2);

    // Methods for inference
    std::vector<float *> search_results(float *data_point, size_t k, unsigned int region_id);
};
