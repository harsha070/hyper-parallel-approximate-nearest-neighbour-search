#include <iostream>
#include <queue>
#include <math.h>
#include <vector>
#include <algorithm>

float l2_distance(std::vector<float> v1, std::vector<float> v2, int dimension)
{
    float dist = 0.0;
    for (size_t i = 0; i < dimension; i++)
    {
        float diff = v1[i] - v2[i];
        dist += diff * diff;
    }
    return sqrt(dist);
}

std::vector<std::vector<float>> drop_duplicates(std::vector<std::vector<float>> aggregated_search_results){
    std::vector<std::vector<float>> results;
    for(int i = 0; i < aggregated_search_results.size(); i++){
        bool is_duplicate = false;
        for(int j = 0; j < i; j++){
            if(aggregated_search_results[i] == aggregated_search_results[j]){
                is_duplicate = true;
                break;
            }
        }
        if(!is_duplicate)
            results.push_back(aggregated_search_results[i]);
    }
    return results;
}


std::vector<std::vector<float>> get_top_k(std::vector<std::vector<float>> points, float* query_point, int k, int dimension)
{
    // intialise a priority queue
    std::priority_queue<std::pair<float, std::vector<float>>> pq;
    std::vector<float> query_vec(query_point, query_point + dimension);

    // iterate over points and add them to priority queue
    for (auto vec: points)
    {
        // compute l2 distance
        float dist = l2_distance(vec, query_vec, dimension);

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
    std::vector<std::vector<float>> closest_vectors(k);
    for (size_t i = 0; i < k; i++)
    {
        closest_vectors[k - i - 1] = pq.top().second;
        pq.pop();
    }

    return closest_vectors;
}
