#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <fstream> // Add this line to include the header file for std::ofstream

#include "../src/core/SearchTree.cpp"
#include "../src/utils/generate_data.cpp"

int main()
{
  // Print a starting message
  std::cout << "Running WeakScaling.cpp" << std::endl;

  double start, end, time;
  int n_init = 2000;                    // Number of points
  int d = 4096;                         // Dimensionality
  int K = 100;                          // Number of points in terminal region
  int work_init = n_init * log(n_init); // Amount of work per thread
  int max_threads = 32;                  // Max number of threads
  int n_runs = 10;                      // Number of runs

  // Set up output csv file
  int now = omp_get_wtime();

  std::ofstream outfile;                                                    // Define an object of std::ofstream
  outfile.open("weak_output/weak_scaling_" + std::to_string(now) + ".csv"); // Create a file named WeakScaling.csv

  // Write the header of the csv file
  outfile << "n_threads,time,n_points,dimensionality" << std::endl;

  // Warm up the cache for fairness between the runs
  int threads = 1;
  const float **full_data = generateRandomData(n_init, d);
  SearchTree tree(n_init, d, full_data, threads, K);

  for (int threads = 1; threads <= max_threads; threads *= 2)
  {
    // Print the number of threads being used to the console immediately without waiting for std::endl
    std::cout << "Running with " << threads << " threads";

    for (int run = 0; run < n_runs; run++)
    {

      // Set the number of threads
      omp_set_num_threads(threads);

      // Generate random data
      // Have it n log n points for each thread since we are doing a weak scaling,
      // i.e., the amount of work per thread should remain constant as the number of threads increases
      int work = work_init * threads;
      int n = n_init;

      while (n * log2(n) <= work)
      {
        n *= 1.2;
      }
      // int n = n_init * threads;
      const float **full_data = generateRandomData(n, d);

      start = omp_get_wtime();

      // ==================== Run the SearchTree building ====================

      SearchTree tree(n, d, full_data, threads, K);

      // ==================== Run the SearchTree building ====================

      end = omp_get_wtime();
      time = end - start;

      // Write the data to the csv file
      outfile << threads << "," << time << "," << n << "," << d << std::endl;
    }

    std::cout << " - Time: " << time << " seconds" << std::endl;
  }

  return 0;
}
