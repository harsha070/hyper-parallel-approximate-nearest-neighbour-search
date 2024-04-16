#include <iostream>
#include <random>
#include <vector>

const float **generateRandomData(int n, int d)
{
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    // Allocate memory for full_data array
    float **full_data = new float *[n];

    // Generate random points in d dimensions
    for (int i = 0; i < n; i++)
    {
        full_data[i] = new float[d];
        for (int j = 0; j < d; j++)
        {
            full_data[i][j] = dis(gen);
        }
    }

    return const_cast<const float **>(full_data);
}

float random_float(float a, float b)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);

    return dis(gen);
}

void deallocateData(int n, float **data)
{
    for (int i = 0; i < n; i++)
    {
        delete[] data[i];
    }
    delete[] data;
}
