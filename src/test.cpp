#include <cassert>
#include <iostream>
#include "core/SearchTree.hpp"
#include <cstdlib>
#include <random>
#include <fstream>
#include <omp.h>

float random_float(float a, float b)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);

    return dis(gen);
}

int main()
{

    // Define the input data
    int n_threads = 8;
    int n_points = 1000;
    int dimension = 128;
    int K = 100;

    std::cout << "Making data" << std::endl;

    // Create some example data
    float **full_data = new float *[n_points];

    for (int i = 0; i < n_points; i++)
    {
        full_data[i] = new float[dimension];
    }

    // Fill the array with random values
    std::srand(time(nullptr));
    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            full_data[i][j] = random_float(-5, 5);
        }
    }

    std::cout << "Making tree" << std::endl;

    // Start the timer
    double start = omp_get_wtime();
    SearchTree search_tree = SearchTree(n_points, dimension, (const float **)full_data, n_threads, K);
    // Stop the timer
    double end = omp_get_wtime();

    // Print the time taken
    std::cout << "Time taken: " << end - start << " seconds" << std::endl;

    std::cout << "Test passed!" << std::endl;
    int total_points = 0;

    for (auto region_id : search_tree.terminal_regions)
    {
        // std::cout << "terminal region id: " << region_id << " count: " << search_tree.regions_reduced_counts[region_id] << std::endl;
        total_points += search_tree.regions_reduced_counts[region_id];
    }

    std::unordered_map<unsigned int, std::vector<float *>> results;
    for (std::unordered_set<unsigned int>::const_iterator i = search_tree.terminal_regions.begin(); i != search_tree.terminal_regions.end(); ++i)
    {
        results.insert(std::make_pair(*i, std::vector<float *>()));
    }

    std::cout << "total points: " << total_points << std::endl;

    float *query_point = new float[dimension];
    for (int j = 0; j < dimension; j++)
    {
        query_point[j] = random_float(-5, 5);
    }
    std::vector<float *> result = search_tree.search_results(query_point, (size_t) 2, 0);

    for (int i = 0; i < (int)result.size(); i++)
    {
        std::cout << "point: " << i << " dist: " << search_tree.l2_distance(query_point, result[i]) << std::endl;
    }

    for (int i = 0; i < n_points; i++)
    {
        delete[] full_data[i];
    }
    delete[] full_data;
    delete[] query_point;

    return 0;
}
