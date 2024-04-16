#include <mpi.h>
#include <vector>
#include "SearchUtil.cpp"
#include "SearchForest.hpp"


SearchForest::SearchForest(int n_trees, int n_points, size_t dimension, const float **full_data, int max_threads, int K,  int argc, char* argv[])
    : n_trees_(n_trees), dimension(dimension)
{
    // resize trees_ vector to number of trees
    trees_.resize(n_trees_);
    // set mpi finalised to false
    is_mpi_finalised = false;
    // Initialize MPI with thread support
    int provided_thread_level;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided_thread_level);

    // Check if the required level of support is provided
    if (provided_thread_level < MPI_THREAD_FUNNELED) {
        std::cerr << "Error: MPI implementation does not provide the required level of thread support (MPI_THREAD_FUNNELED)." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int rank_, size_;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    // Create and distribute trees across MPI ranks
    for (int i = 0; i < n_trees; ++i)
    {
        if (i % size_ == rank_)
        {
            SearchTree* new_tree = new SearchTree(n_points, dimension, (const float **)full_data, max_threads, K);
            trees_[i] = new_tree;
        }
    }
}

std::vector<std::vector<float>> SearchForest::search_results(float *query_point, size_t k)
{
    int rank_, size_;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    std::cout << "world size: " << size_ << " rank: " << rank_ << " num trees: " << trees_.size() << std::endl;

    std::vector<std::vector<float>> aggregated_search_results(trees_.size() * k, std::vector<float>(dimension));

    // search trees in a distributed way across MPI ranks
    for (int i = 0; i < trees_.size(); ++i)
    {
        if (i % size_ == rank_)
        {
            std::vector<float *> result = trees_[i]->search_results(query_point, (size_t) k, 0);
            std::cout << "searched tree " << i << " found results: " << result.size() << std::endl;
            for(int j = 0; j < k; j++){
                for(int ii = 0; ii < dimension; ii++)
                    aggregated_search_results[i * k + j][ii] = result[j][ii];
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    is_mpi_finalised = true;

    // vector for final results
    std::vector<std::vector<float>> results;

    // drop duplicate results
    results = drop_duplicates(aggregated_search_results);
    // get top k candidates using a priority queue
    results = get_top_k(results, query_point, k, dimension);

    return results;
}

SearchForest::~SearchForest()
{
    // Deallocate the trees when we are done
    for (SearchTree* tree : trees_) {
        delete tree;
    }

    // Finalize MPI if not called before
    if (!is_mpi_finalised){
        is_mpi_finalised = true;
        MPI_Finalize();
    }
}

