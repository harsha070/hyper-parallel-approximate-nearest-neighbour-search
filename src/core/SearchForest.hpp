#include <vector>
#include "SearchTree.hpp"

class SearchForest
{
public:
    SearchForest(int n_trees, int n_points, size_t dimension, const float **full_data, int n_threads, int K,  int argc, char* argv[]);
    ~SearchForest();
    std::vector<std::vector<float>> search_results(float *query_point, size_t k);

private:
    int n_trees_;
    std::vector<SearchTree *> trees_;
    int dimension;
    bool is_mpi_finalised = false;
};
