#include "core/SearchForest.hpp"
#include "utils/generate_data.cpp"
#include <cstddef>

int main(int argc, char* argv[])
{
    // Define number of trees in the forest
    int n_trees = 2;
    // Define number of test data points
    int n_points = 100;
    // Dimensionality of vectors
    size_t dimension = 2;
    // Synthesize test data
    const float **full_data = generateRandomData(n_points, dimension);
    // Number of OpenMP threads
    int n_threads = 1;
    // Number of nearest neighbours to fetch
    int K = 10;

    // Build the search forest
    SearchForest forest(n_trees, n_points, dimension, full_data, n_threads, K, argc, argv);

    // Free memory
    deallocateData(n_points, (float **)full_data);

    return 0;
}
