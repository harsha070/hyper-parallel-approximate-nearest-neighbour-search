#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <fstream>

#include "../src/core/SearchTree.cpp"
#include "../src/utils/generate_data.cpp"

int main(int argc, char *argv[])
{

  double start, end, time;
  int n = 30000; // Number of points
  int d = 4096;  // Dimensionality
  int K = 100;   // Number of points in terminal region
  int n_runs = 50;
  int max_threads = 8;

  // Parse the command-line argument
  if (argc > 1)
  {
    max_threads = std::stoi(argv[1]);
    std::cout << "Using max_threads: " << max_threads << std::endl;
  }
  else
  {
    std::cout << "No command-line argument provided for max_threads. Using default value: " << max_threads << std::endl;
  }

  // Print a starting message
  std::cout << "Running Strong Scaling Test with " << n << " points in " << d << " dimensions" << std::endl;
  std::cout << "Number of runs: " << n_runs << " - Max number of threads: " << max_threads << std::endl;

  // Set up output csv file
  int now = omp_get_wtime();
  std::ofstream outfile;                                                        // Define an object of std::ofstream
  outfile.open("strong_output/strong_scaling_" + std::to_string(now) + ".csv"); // Create a file named WeakScaling.csv

  // Write the header of the csv file
  outfile << "n_threads,time,n_points,dimensionality" << std::endl;

  // Generate random data
  const float **full_data = generateRandomData(n, d);

  // Warm up the cache for fairness between the runs
  int threads = 1;
  SearchTree tree(n, d, full_data, threads, K);

  for (threads = 1; threads <= max_threads; threads *= 2)
  {
    // Print the number of threads being used to the console immediately without waiting for std::endl
    std::cout << "Running with " << threads << " threads";

    // Set the number of threads
    omp_set_num_threads(threads);

    for (int run = 0; run < n_runs; run++)
    {

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
