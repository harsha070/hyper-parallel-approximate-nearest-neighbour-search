#include "core/SearchForest.hpp"
#include "utils/generate_data.cpp"
#include <cstddef>

int main(int argc, char* argv[])
{
    int n_trees = 2;
    int n_points = 100;
    size_t dimension = 2;
    const float **full_data = generateRandomData(n_points, dimension);
    int n_threads = 1;
    int K = 10;

    SearchForest forest(n_trees, n_points, dimension, full_data, n_threads, K, argc, argv);

    deallocateData(n_points, (float **)full_data);

    return 0;
}
