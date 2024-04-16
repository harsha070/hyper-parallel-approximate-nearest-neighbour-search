#include "core/SearchForest.hpp"
#include "utils/generate_data.cpp"
#include <cstddef>

int main(int argc, char* argv[])
{
    int n_trees = 2;
    int n_points = 100;
    size_t dimension = 2;
    const float **full_data = generateRandomData(n_points, dimension);
    int n_threads = 2;
    int K = 10;

    SearchForest forest(n_trees, n_points, dimension, full_data, n_threads, K, argc, argv);

    std::cout << "Building forest complete!" << std::endl;

    size_t k_neighbours = 10;
    float *query_point = new float[dimension];
    for (int j = 0; j < (int) dimension; j++)
        query_point[j] = random_float(-5, 5);
    std::vector<std::vector<float>> result = forest.search_results(query_point, k_neighbours);
    std::cout << "Search complete! results found: " << result.size() << std::endl;

    deallocateData(n_points, (float **)full_data);

    return 0;
}
