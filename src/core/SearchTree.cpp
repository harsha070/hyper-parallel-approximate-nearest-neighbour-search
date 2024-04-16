#include <unordered_set>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <unordered_set>
#include <random>
#include <queue>
#include <cassert>
#include <cstring>
#include <cassert>

#include "SearchTree.hpp"

size_t hp_sz = 0;

SearchTree::SearchTree(const int n_points, const size_t dimension, const float **full_data, const int n_threads, int K)
    : n_points(n_points), D(dimension), full_data(full_data), n_threads(n_threads), K(K)
{
    size_t regions_reduced_counts_sz = 6 * n_points;
    size_t hyperplanes_sz = 6 * n_points;

    hp_sz = hyperplanes_sz;

    data = new float *[n_points];
    regions_curr = new unsigned int[n_points];
    regions_next = new unsigned int[n_points];
    regions_counts = new unsigned int *[n_threads];
    regions_reduced_counts = new unsigned int[regions_reduced_counts_sz];
    
    for (int i = 0; i < n_threads; ++i)
    {
        regions_counts[i] = new unsigned int[6 * n_points]();
    }

    std::mt19937 gen(rd());

    // Allocate memory for the sampled points vector
    sampled_points = new std::pair<float *, float *>[n_points];

    // Preallocate the memory for the hyperplanes
    hyperplanes.resize(hyperplanes_sz); // Allocate memory for the hyperplanes vector

    // initialize buffers
    std::fill(regions_curr, regions_curr + n_points, 0);
    std::fill(regions_reduced_counts, regions_reduced_counts + regions_reduced_counts_sz, 0);

    // TODO: this should be done with omp parallel for
    for (int i = 0; i < n_threads; ++i)
    {
        std::fill(regions_counts[i], regions_counts[i] + n_points, 0);
    }

    // initialize active regions. R_0 is the only active region at the beginning
    active_regions.push_back(0);
    regions_reduced_counts[0] = n_points;

    // initialze number of active points
    n_active_points_ = n_points;

    // TODO: initialize data
    std::memcpy(data, full_data, n_points * sizeof(float *));

    // Set up variables for keeping track of number of calls to the signed distance function
    // and the wall time spent in the signed distance function

    construct_tree();
}

SearchTree::~SearchTree()
{
    delete[] data;
    delete[] regions_curr;
    delete[] regions_next;
    delete[] regions_reduced_counts;

    for (int i = 0; i < n_threads; ++i)
    {
        delete[] regions_counts[i];
    }
    delete[] regions_counts;

    for (auto hp : hyperplanes)
    {
        delete[] hp.normal_vec;
    }

    delete[] sampled_points;
}

void SearchTree::construct_tree()
{
    // iterate until all points are terminal
    while (true)
    {
        sample_points(active_regions, regions_reduced_counts, sampled_points);

        // create hyperplanes in parallel for each of the active regions
#pragma omp parallel for
        for (auto region_id : active_regions)
        {

            assert(region_id < hp_sz);
            make_hyperplane_(sampled_points[region_id], hyperplanes[region_id]);
            // associate region with hyperplane
#pragma omp critical
            {
                region_hyperplane_map.insert(std::make_pair(region_id, &hyperplanes[region_id]));
            }
        }

        // classify points in parallel
#pragma omp parallel num_threads(n_threads)
        {
            int m = omp_get_num_threads();            // number of threads
            int tid = omp_get_thread_num();           // my thread number
            int partial_terms = n_active_points_ / m; // number of terms this rank has to sum
            int my_first = tid * partial_terms;       // the first term of my partial sum
            int my_last = my_first + partial_terms;   // the last term of my partial sum

            // distribute the work imbalance evenly among ranks.  Hint: you need to
            // modify my_first and my_last for specific ranks.
            const int rem = n_active_points_ % m; // remainder

            if (tid < rem)
            {
                my_first += tid;
                my_last += tid + 1;
            }
            else
            {
                my_first += rem;
                my_last += rem;
            }

            classify(my_first, my_last, tid);
        }

        // reduce the per-thread d_ir +ego * n counts and mark any terminal regions
        reduce_counts(regions_reduced_counts);

        // Check which regions are active
        update_active_regions(active_regions, regions_reduced_counts);

        // prepare for next iteration
        for (const auto &tr : new_terminal_regions)
        {
            terminal_regions.insert(tr);
        }
        new_terminal_regions.clear();

        // if no active regions remain, then algorithm terminates
        if (active_regions.size() == 0)
        {
            break;
        }
        std::swap(regions_curr, regions_next);
    }

    for (auto region : terminal_regions)
    {
        terminal_region_data_store[region] = {};
    }

    for (int i = 0; i < n_points; i++)
    {
        int region_id = regions_next[i];
        terminal_region_data_store[region_id].push_back(data[i]);
    }
}

/**
 * @brief
 *
 * Classify the points in the dataset and assign them to regions by
 * setting the corresponding region ID in the regions_next array
 *
 * @param start_idx
 * @param end_idx
 * @return void
 */

void SearchTree::classify(int start_idx, int end_idx, int tid)
{

    for (int i = start_idx; i < end_idx; ++i)
    {
        const float *point = full_data[i];

        // get the region ID for this point
        int region_id = regions_curr[i];
        if (terminal_regions.find(region_id) != terminal_regions.end())
        {
            regions_next[i] = region_id;
            continue;
        }

        const float region_above = 2 * region_id + 1;
        const float region_below = 2 * region_id + 2;
        Hyperplane *hp = region_hyperplane_map[region_id];

        // get the sign of the distance from the hyperplane for this point
        float distance = signed_distance(*hp, point);

        // classify the point and update count
        if (distance > 0)
        {
            regions_next[i] = region_above;
            inc_region_count_(tid, region_above);
        }
        else
        {
            regions_next[i] = region_below;
            inc_region_count_(tid, region_below);
        }
    }
}

/**
 * @brief
 *
 * @param reduced_region_counts
 * @param terminal_regs
 * @return void
 */

void SearchTree::reduce_counts(unsigned int *reduced_region_counts)
{
    for (int r = 0; r < active_regions.size(); r++)
    {
        // reduce the per-thread counts
        for (int t = 0; t < n_threads; t++)
        {
            reduced_region_counts[2 * active_regions[r] + 1] += get_region_count(t, 2 * active_regions[r] + 1);
            reduced_region_counts[2 * active_regions[r] + 2] += get_region_count(t, 2 * active_regions[r] + 2);
        }

        // mark a region as terminal if there are <= K points in the region
        if (reduced_region_counts[2 * active_regions[r] + 1] <= K)
        {
            new_terminal_regions.insert(2 * active_regions[r] + 1);
        }
        if (reduced_region_counts[2 * active_regions[r] + 2] <= K)
        {
            new_terminal_regions.insert(2 * active_regions[r] + 2);
        }
    }
}

/**
 * @brief thread-safe load balancing
 *
 *
 * @return void
 */
void SearchTree::move_terminal_points_(unsigned int start_idx, unsigned int end_idx_exclusive,
                                       std::unordered_set<unsigned int> &terminal_regions)
{
    for (unsigned int i = start_idx; i < end_idx_exclusive; i++)
    {
        if (terminal_regions.find(regions_next[i]) != terminal_regions.end())
        {
#pragma omp critical
            {
                // swap terminal and non-terminal points and associated regions
                std::swap(data[i], data[n_active_points_ - 1]);
                std::swap(regions_next[i], regions_next[n_active_points_ - 1]);

                // decrement active points counter. Everything beyond this index is terminal
                n_active_points_--;
            }
        }
    }
}

/**
 * @brief Update vector of active regions
 *
 * Active regions are ones that we may be classifying points into.
 * At iteration i, if region r was active and is non-terimal, then it should
 * become inactive and regions 2i + 1 and 2i + 2 should become active
 *
 */
void SearchTree::update_active_regions(std::vector<unsigned int> &active_regions, unsigned int *reduced_region_counts)
{

    // clear the active regions vector
    // int _n_active_regions = active_regions.size();
    std::vector<unsigned int> new_active_regions;

    for (auto region_id : active_regions)
    {
        if (reduced_region_counts[2 * region_id + 1] > K)
        {
            new_active_regions.push_back(2 * region_id + 1);
        }
        if (reduced_region_counts[2 * region_id + 2] > K)
        {
            new_active_regions.push_back(2 * region_id + 2);
        }
    }
    active_regions.clear();
    for (auto region_id : new_active_regions)
    {
        active_regions.push_back(region_id);
    }
}

/**
 * @brief Sample a pair of points from each active region
 *
 * The sampled points will be used to construct the hyperplanes for each of the active regions.
 *
 * Suppose region r has 10 points in it. Sampling is done by selescting two non-identical values i and j from
 * integers (0, 10). Then we do a pass through all points and select the ith and jth points that we've seen.
 *
 */
void SearchTree::sample_points(const std::vector<unsigned int> &active_regions, const unsigned int *regions_reduced_counts, std::pair<float *, float *> *sampled_points)
{
    int n_active_regions = active_regions.size();

    // Allocate int array of size 2 * n_active_regions
    auto sampled_sequence_nums = new std::pair<int, int>[6 * n_points];

    std::unordered_map<unsigned int, unsigned int> seen;

    // For each region, sample two points from the region in the range(0, reduced_region_counts[region])
#pragma omp parallel for
    for (int i = 0; i < n_active_regions; i++)
    {
        int n_points_in_region = regions_reduced_counts[active_regions[i]];
        std::uniform_int_distribution<> dis(0, n_points_in_region - 1);

        int index1, index2;
        do
        {
            index1 = dis(gen);
            index2 = dis(gen);
        } while (index1 == index2); // Repeat until two unique indices are generated

        sampled_sequence_nums[active_regions[i]] = std::make_pair(index1, index2);
        seen[active_regions[i]] = 0;
    }

    // Iterate over all points and check if they are one of the two sampled points
    // for their region. If so, add them to the sampled points array.
    for (int i = 0; i < n_points; i++)
    {
        int region = regions_curr[i];

        if (terminal_regions.find(region) != terminal_regions.end())
        {
            continue;
        }

        int n_seen = seen[region];

        if (n_seen == sampled_sequence_nums[region].first)
        {
            sampled_points[region].first = data[i];
        }
        else if (n_seen == sampled_sequence_nums[region].second)
        {
            sampled_points[region].second = data[i];
        }
        ++seen[region];
    }

    delete[] sampled_sequence_nums;
}

/**
 * @brief Calculates the sign of the distance between the a hyperplane and a point
 *
 * This is used to classify which side of the hyperplane the point is on
 *
 */
inline float SearchTree::signed_distance(Hyperplane &hp, const float *point)
{
    float distance = hp.offset;
#pragma omp simd reduction(+ : distance)
    for (size_t d = 0; d < D; d++)
    {
        distance += hp.normal_vec[d] * point[d];
    }

    return distance;
}

/**
 * @brief Construct a hyperplane from 2 points
 *
 */
void SearchTree::make_hyperplane_(std::pair<float *, float *> &sampled_point, Hyperplane &hp)
{
    // TODO: reimplement with intel intrinsics or at minimum use loop unrolling
    const float *p1 = sampled_point.first;
    const float *p2 = sampled_point.second;

    std::vector<float> midpoint(D, 0);
    hp.offset = 0.0;
    for (unsigned int i = 0; i < D; i++)
    {
        // Calculate the midpoint between the two points
        midpoint[i] = (p1[i] + p2[i]) / 2.0;

        // Calculate the normal vector of the hyperplane
        hp.normal_vec[i] = (p2[i] - p1[i]);

        // Calculate the offset of the hyperplane from the origin
        hp.offset -= hp.normal_vec[i] * midpoint[i];
    }
}

/**
 * @brief Increments the per-thread region count for some region
 *
 * This is an interface for updating whatever data structure will be used
 * to store the per-thread region counts
 *
 * @param tid the thread ID
 * @param region the region ID
 */
inline void SearchTree::inc_region_count_(const int tid, const int region)
{
    regions_counts[tid][region]++;
}

/**
 * @brief Gets the number of points a thread has seen in some region
 *
 * This is an interface for getting counts from whatever data structure will be used
 * to store the per-thread region counts
 *
 * @param tid
 * @param region the region ID
 * @return
 */
inline int SearchTree::get_region_count(const int tid, const int region)
{
    return regions_counts[tid][region];
}

float SearchTree::l2_distance(float *v1, float *v2)
{
    float dist = 0.0;
    for (size_t i = 0; i < D; i++)
    {
        float diff = v1[i] - v2[i];
        dist += diff * diff;
    }
    return sqrt(dist);
}

/**
 * @brief ANN Search tree inference
 *
 * Given a point, find the k approximate nearest neighbors
 *
 */
std::vector<float *> SearchTree::search_results(float *data_point, size_t k, unsigned int region_id)
{

    if (k > regions_reduced_counts[region_id])
    {
        k = regions_reduced_counts[region_id];
        std::cout << "Requested too many neighbours! limiting to " << k << " neighbours." << std::endl;
    }

    std::vector<float *> result(0);
    if (k == 0)
    {
        return result;
    }

    std::cout << "Searching region id: " << region_id << std::endl;

    // check if we are in a terminal region
    if (terminal_region_data_store.find(region_id) != terminal_region_data_store.end())
    {

        // we are in a terminal region
        std::cout << "Entered a terminal region!" << std::endl;

        // get all points in the terminal region
        std::vector<float *> points = terminal_region_data_store[region_id];

        // check assertion
        assert(k <= points.size());

        // intialise a priority queue
        std::priority_queue<std::pair<float, float *>> pq;

        // iterate over points and add them to priority queue
        for (auto vec : points)
        {

            // compute l2 distance
            float dist = l2_distance(vec, data_point);

            // Add the current vector to the priority queue if it is closer than the kth closest vector seen so far
            if (pq.size() < k)
            {
                pq.push(std::make_pair(dist, vec));
            }
            else if (dist < pq.top().first)
            {
                pq.pop();
                pq.push(std::make_pair(dist, vec));
            }
        }

        // Extract the k closest vectors from the priority queue and return them as a vector
        std::vector<float *> closest_vectors(k);
        for (size_t i = 0; i < k; i++)
        {
            closest_vectors[k - i - 1] = pq.top().second;
            pq.pop();
        }

        return closest_vectors;
    }
    Hyperplane hp = *(region_hyperplane_map[region_id]);
    float dist = signed_distance(hp, data_point);

    size_t k_left_subtree = 0;
    size_t k_right_subtree = 0;

    if (dist < 0)
    {
        size_t left_subtree_count = regions_reduced_counts[2 * region_id + 1];
        k_left_subtree = std::min(left_subtree_count, k);
        k_right_subtree = k - k_left_subtree;
    }
    else
    {
        size_t right_subtree_count = regions_reduced_counts[2 * region_id + 2];
        k_right_subtree = std::min(right_subtree_count, k);
        k_left_subtree = k - k_right_subtree;
    }

    std::vector<float *> result_left = search_results(data_point, k_left_subtree, 2 * region_id + 1);
    std::vector<float *> result_right = search_results(data_point, k_right_subtree, 2 * region_id + 2);

    for (auto vec : result_left)
        result.push_back(vec);
    for (auto vec : result_right)
        result.push_back(vec);
    return result;
}
