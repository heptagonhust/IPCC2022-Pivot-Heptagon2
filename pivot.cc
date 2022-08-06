#include "defs.hpp"
#include <cmath>
#include <cstdlib>
#include <map>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <thread>
#include <vector>

// Calculate sum of distance while combining different pivots. Complexity : O( n^2 )
double SumDistance(const int k, const int n, const int dim, double *coord, int *pivots) {
  double *rebuiltCoord = (double *)malloc(sizeof(double) * n * k);
  int i;
  for (i = 0; i < n * k; i++) {
    rebuiltCoord[i] = 0;
  }

  // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
  for (i = 0; i < n; i++) {
    int ki;
    for (ki = 0; ki < k; ki++) {
      double distance = 0;
      int pivoti = pivots[ki];
      int j;
      for (j = 0; j < dim; j++) {
        distance += pow(coord[pivoti * dim + j] - coord[i * dim + j], 2);
      }
      rebuiltCoord[i * k + ki] = sqrt(distance);
    }
  }

  // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
  double chebyshevSum = 0;
  for (i = 0; i < n; i++) {
    int j;
    for (j = 0; j < n; j++) {
      double chebyshev = 0;
      int ki;
      for (ki = 0; ki < k; ki++) {
        double dis = fabs(rebuiltCoord[i * k + ki] - rebuiltCoord[j * k + ki]);
        chebyshev = dis > chebyshev ? dis : chebyshev;
      }
      chebyshevSum += chebyshev;
    }
  }

  free(rebuiltCoord);

  return chebyshevSum;
}

// Recursive function recursive_combinations() : combine pivots and calculate the sum of distance while combining different pivots.
// ki  : current depth of the recursion
// k   : number of pivots
// n   : number of points
// dim : dimension of metric space
// M   : number of combinations to store
// coord  : coordinates of points
// pivots : indexes of pivots
// maxDistanceSum  : the largest M distance sum
// maxDisSumPivots : the top M pivots combinations
// minDistanceSum  : the smallest M distance sum
// minDisSumPivots : the bottom M pivots combinations
void recursive_combinations(int ki, const int k, const int n, const int dim, const int M, double *coord, int *pivots, double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots) {
  if (ki == k - 1) {
    int i;
    for (i = pivots[ki - 1] + 1; i < n; i++) {
      pivots[ki] = i;

      // Calculate sum of distance while combining different pivots.
      double distanceSum = SumDistance(k, n, dim, coord, pivots);

      // put data at the end of array
      maxDistanceSum[M] = distanceSum;
      minDistanceSum[M] = distanceSum;
      int kj;
      for (kj = 0; kj < k; kj++) {
        maxDisSumPivots[M * k + kj] = pivots[kj];
      }
      for (kj = 0; kj < k; kj++) {
        minDisSumPivots[M * k + kj] = pivots[kj];
      }
      // sort
      int a;
      for (a = M; a > 0; a--) {
        if (maxDistanceSum[a] > maxDistanceSum[a - 1]) {
          double temp = maxDistanceSum[a];
          maxDistanceSum[a] = maxDistanceSum[a - 1];
          maxDistanceSum[a - 1] = temp;
          int kj;
          for (kj = 0; kj < k; kj++) {
            int temp = maxDisSumPivots[a * k + kj];
            maxDisSumPivots[a * k + kj] = maxDisSumPivots[(a - 1) * k + kj];
            maxDisSumPivots[(a - 1) * k + kj] = temp;
          }
        }
      }
      for (a = M; a > 0; a--) {
        if (minDistanceSum[a] < minDistanceSum[a - 1]) {
          double temp = minDistanceSum[a];
          minDistanceSum[a] = minDistanceSum[a - 1];
          minDistanceSum[a - 1] = temp;
          int kj;
          for (kj = 0; kj < k; kj++) {
            int temp = minDisSumPivots[a * k + kj];
            minDisSumPivots[a * k + kj] = minDisSumPivots[(a - 1) * k + kj];
            minDisSumPivots[(a - 1) * k + kj] = temp;
          }
        }
      }
    }
    return;
  }

  // Recursively call Combination() to combine pivots
  int i;
  for (i = pivots[ki - 1] + 1; i < n; i++) {
    pivots[ki] = i;
    recursive_combinations(ki + 1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
    *** You can delete the logging code. **/
    // if(ki==k-2){
    //     int kj;
    //     for(kj=0; kj<k; kj++){
    //         printf("%d ", pivots[kj]);
    //     }
    //     putchar('\t');
    //     for(kj=0; kj<k; kj++){
    //         printf("%d ", maxDisSumPivots[kj]);
    //     }
    //     printf("%lf\t", maxDistanceSum[0]);
    //     for(kj=0; kj<k; kj++){
    //         printf("%d ", minDisSumPivots[kj]);
    //     }
    //     printf("%lf\n", minDistanceSum[0]);
    // }
  }
}

int next_comb(int *arr, int n, int k) {
  int pos = k - 1;
  while (pos >= 0 && n - arr[pos] == k - pos) {
    pos = pos - 1;
  }
  if (pos == -1) {
    return -1;
  }
  int ret = pos;
  arr[pos]++;
  for (int i = pos + 1; i < k; i++) {
    arr[i] = arr[i - 1] + 1;
  }
  return ret;
}

i64 choose(i32 n, i32 k) {
  if (n < k) {
    return 0;
  }
  i64 ans = 1;
  for (int i = n; i >= (n - k + 1); --i)
    ans *= i;
  while (k)
    ans /= k--;
  return ans;
}

void mth_comb(i32 *arr, i32 n, i32 k, i32 m) {
  i32 res_ptr = 0;
  i32 a = n;
  i32 b = k;
  i32 x = (choose(n, k) - 1) - m;
  for (i32 i = 0; i < k; ++i) {
    a -= 1;
    while (choose(a, b) > x) {
      a -= 1;
    }
    arr[res_ptr++] = n - 1 - a;
    x -= choose(a, b);
    b -= 1;
  }
}

float distance(const double *coord, int ndims, int x, int y) {
  double dist = .0;
  for (int i = 0; i < ndims; i++) {
    dist += (coord[ndims * x + i] - coord[ndims * y + i]) * (coord[ndims * x + i] - coord[ndims * y + i]);
  }
  return sqrt(dist);
}

double calc_value(const int prev, const int npoints, const int npivots, const int ndims, int *pivots, const double *coord, float *rebuilt_coord, float *mx) {
  // Part 1. Rebuild Coordintate System
  for (int k = prev; k < npivots; k++) {
    int p = pivots[k];
    for (int i = 0; i < npoints; i++) {
      rebuilt_coord[k * npoints + i] = distance(coord, ndims, p, i);
    }
  }

  // Part 2. Evaluate System Value

  // Part 2.1. Calculate Chebyshev Distance

  int points_pairs = npoints * (npoints + 1) / 2;

  for (int k = prev; k < npivots - 1; k++) {
    int idx_cnt = 0;
    for (int i = 0; i < npoints; i++) {
      for (int j = 0; j < i; j++) {
        float chebyshev_dim_dist = fabs(rebuilt_coord[k * npoints + i] - rebuilt_coord[k * npoints + j]);
        if (k > 0 && chebyshev_dim_dist < mx[(k - 1) * points_pairs + idx_cnt + j]) {
          chebyshev_dim_dist = mx[(k - 1) * points_pairs + idx_cnt + j];
        }
        mx[k * points_pairs + idx_cnt + j] = chebyshev_dim_dist;
      }
      idx_cnt += i + 1;
    }
  }

  // Part 2.2. Last loop and Get Sum
  // k == npivots - 1
  double chebyshev_dist_sum = .0;
  int last = npivots - 1;
  int idx_cnt = 0;
  for (int i = 0; i < npoints; i++) {
    for (int j = 0; j < i; j++) {
      float chebyshev_dim_dist = fabs(rebuilt_coord[last * npoints + i] - rebuilt_coord[last * npoints + j]);
      if (chebyshev_dim_dist < mx[(last - 1) * points_pairs + idx_cnt + j]) {
        chebyshev_dim_dist = mx[(last - 1) * points_pairs + idx_cnt + j];
      }
      mx[last * points_pairs + idx_cnt + j] = chebyshev_dim_dist;
      chebyshev_dist_sum += chebyshev_dim_dist;
    }
    idx_cnt += i + 1;
  }

  // Calculate Half of All Pairs, Then Double
  return chebyshev_dist_sum * 2;
}

struct MinMaxPivotPtrs {
  int *minDisSumPivots;
  int *maxDisSumPivots;
  double *minDistanceSum;
  double *maxDistanceSum;
};

// maxDisSum, minDisSum, maxDisSumPivots, minDisSumPivots
// run as a thread
void combinations(const int num_total_threads, const int blocks, const int cnk, const int thread_id, const int npoints, const int npivots, const int ndims, const int M, const double *coord, MinMaxPivotPtrs *ptrs) {
  int *minDisSumPivots = (int *)malloc(sizeof(int) * M * npivots);
  int *maxDisSumPivots = (int *)malloc(sizeof(int) * M * npivots);
  double *minDistanceSum = (double *)malloc(sizeof(double) * M);
  double *maxDistanceSum = (double *)malloc(sizeof(double) * M);

  ptrs->minDisSumPivots = minDisSumPivots;
  ptrs->maxDisSumPivots = maxDisSumPivots;
  ptrs->minDistanceSum = minDistanceSum;
  ptrs->maxDistanceSum = maxDistanceSum;

  int points_pairs = npoints * (npoints + 1) / 2;
  int chips = (blocks * num_total_threads);
  int workload = cnk / chips;
  if (cnk % chips != 0) {
    workload++;
  }

  struct timeval start, end;
  gettimeofday(&start, NULL);

  float *rebuilt_coord = (float *)malloc(sizeof(float) * npivots * npoints);
  float *mx = (float *)malloc(sizeof(float) * npivots * points_pairs);
  int *maxTmpPivots = (int *)malloc(sizeof(int) * M * npivots);
  int *minTmpPivots = (int *)malloc(sizeof(int) * M * npivots);
  std::map<double, int> mx_mp{};
  std::map<double, int> mn_mp{};

  int pivots_cnt = 0;

  for (int b = 0; b < blocks; b++) {
    int chip_id = b * num_total_threads + thread_id;
    int start_point = chip_id * workload;
    int end_point = start_point + workload;
    if (chip_id + 1 == num_total_threads * blocks) {
      end_point = cnk;
    }
    // fprintf(stderr, "thread: %d, chip: %d, block: %d\n", thread_id, chip_id, b);
    int pivots[npivots];
    mth_comb(pivots, npoints, npivots, start_point);
    int prev = 0;
    for (int comb_cnt = start_point; comb_cnt < end_point && prev != -1; ++comb_cnt) {
      double value = calc_value(prev, npoints, npivots, ndims, pivots, coord, rebuilt_coord, mx);
      pivots_cnt += npivots - prev;
      // Part 3. Get Top M and Bottom M

      // Part 3.1. Top M cases
      int size = (int)mx_mp.size();
      if (size < M) {
        int idx = (int)mx_mp.size();
        mx_mp.insert(std::map<double, int>::value_type(value, idx));
        for (int i = 0; i < npivots; i++) {
          maxTmpPivots[idx * npivots + i] = pivots[i];
        }
      } else {
        auto iter = mx_mp.begin();
        if (iter->first < value) {
          int idx = iter->second;
          mx_mp.erase(iter);
          mx_mp.insert(std::map<double, int>::value_type(value, idx));
          for (int i = 0; i < npivots; i++) {
            maxTmpPivots[idx * npivots + i] = pivots[i];
          }
        }
      }

      // Part 3.2. Bottom M cases
      size = (int)mn_mp.size();
      if (size < M) {
        int idx = (int)mn_mp.size();
        mn_mp.insert(std::map<double, int>::value_type(value, idx));
        for (int i = 0; i < npivots; i++) {
          minTmpPivots[idx * npivots + i] = pivots[i];
        }
      } else {
        auto iter = mn_mp.end();
        iter--;
        if (iter->first > value) {
          int idx = iter->second;
          mn_mp.erase(iter);
          mn_mp.insert(std::map<double, int>::value_type(value, idx));
          for (int i = 0; i < npivots; i++) {
            minTmpPivots[idx * npivots + i] = pivots[i];
          }
        }
      }

      prev = next_comb(pivots, npoints, npivots);
    }
  }

  // Part 4. Sort Answer_pivots arrays

  int idx = 0;
  for (auto iter = mx_mp.rbegin(); iter != mx_mp.rend(); iter++) {
    int tmp_idx = iter->second;
    for (int i = 0; i < npivots; i++) {
      maxDisSumPivots[idx * npivots + i] = maxTmpPivots[tmp_idx * npivots + i];
    }
    maxDistanceSum[idx] = iter->first;
    idx++;
  }

  idx = 0;
  for (auto iter = mn_mp.begin(); iter != mn_mp.end(); iter++) {
    int tmp_idx = iter->second;
    for (int i = 0; i < npivots; i++) {
      minDisSumPivots[idx * npivots + i] = minTmpPivots[tmp_idx * npivots + i];
    }
    minDistanceSum[idx] = iter->first;
    idx++;
  }
  // fprintf(stderr, "thread %d: pivots count: %d\n", thread_id, pivots_cnt);

  free(maxTmpPivots);
  free(minTmpPivots);
  free(rebuilt_coord);
  free(mx);

  gettimeofday(&end, NULL);
  fprintf(stderr, "Thead %d : %d pivots, %f ms\n", thread_id, pivots_cnt, (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);
}

void Combination(int ki, const int k, const int n, const int dim, const int M, const double *coord, int *pivots, double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots) {
  // 通过调用以下方法来得到所有可以工作的进程数量
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // 得到当前进程的秩
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  u32 threads_per_rank = 32;
  u32 blocks = 2;
  u32 num_total_threads = threads_per_rank * world_size;
  i32 cnk = choose(n, k);
  std::vector<MinMaxPivotPtrs> thread_data(threads_per_rank);
  std::vector<std::thread> threads(threads_per_rank);
  for (u32 i = 0; i < threads_per_rank; ++i) {
    i32 thread_id = i * world_size + world_rank;
    // i32 start_point = workload * thread_id;
    // i32 end_point = start_point + workload;
    // if (thread_id + 1 == num_total_threads) {
    //   end_point = cnk;
    // }
    printf("%d/%d: thread: %d\n", world_rank, world_size, thread_id);
    threads[i] = std::thread([&, i, thread_id, n, k, dim, M, coord] { combinations(num_total_threads, blocks, cnk, thread_id, n, k, dim, M, coord, &thread_data[i]); });
    // bind thread to core
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(thread_id, &cpuset);
    int rc = pthread_setaffinity_np(threads[i].native_handle(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
      printf("Error calling pthread_setaffinity_np: %d", rc);
    }
  }
  for (auto &i : threads) {
    i.join();
  }
  int *MPIMinDisSumPivots;
  int *MPIMaxDisSumPivots;
  double *MPIMinDistanceSum;
  double *MPIMaxDistanceSum;
  if (world_rank == 0) {
    // for rank 0, alloc extra buffer to hold data from other threads
    MPIMinDisSumPivots = (int *)malloc(sizeof(int) * world_size * M * k);
    MPIMaxDisSumPivots = (int *)malloc(sizeof(int) * world_size * M * k);
    MPIMinDistanceSum = (double *)malloc(sizeof(double) * world_size * M);
    MPIMaxDistanceSum = (double *)malloc(sizeof(double) * world_size * M);
  } else {
    MPIMinDisSumPivots = nullptr;
    MPIMaxDisSumPivots = nullptr;
    MPIMinDistanceSum = nullptr;
    MPIMaxDistanceSum = nullptr;
  }

  // reduce thread min max
  std::vector<u32> maxPtr(threads_per_rank), minPtr(threads_per_rank);
  for (int i = 0; i < M; ++i) {
    f64 max_value = -1 / 0.0, min_value = 1 / 0.0;
    i32 max_idx = -1, min_idx = -1;
    for (u32 j = 0; j < threads_per_rank; ++j) {
      if (thread_data[j].maxDistanceSum[maxPtr[j]] > max_value) {
        max_value = thread_data[j].maxDistanceSum[maxPtr[j]];
        max_idx = j;
      }
      if (thread_data[j].minDistanceSum[minPtr[j]] < min_value) {
        min_value = thread_data[j].minDistanceSum[minPtr[j]];
        min_idx = j;
      }
    }
    maxDistanceSum[i] = max_value;
    minDistanceSum[i] = min_value;
    for (int j = 0; j < k; ++j) {
      maxDisSumPivots[i * k + j] = thread_data[max_idx].maxDisSumPivots[maxPtr[max_idx] * k + j];
      minDisSumPivots[i * k + j] = thread_data[min_idx].minDisSumPivots[minPtr[min_idx] * k + j];
    }
    ++maxPtr[max_idx];
    ++minPtr[min_idx];
  }
  // reduce mpi min max
  MPI_Gather(minDistanceSum, M, MPI_DOUBLE, MPIMinDistanceSum, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(maxDistanceSum, M, MPI_DOUBLE, MPIMaxDistanceSum, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(minDisSumPivots, M * k, MPI_INT, MPIMinDisSumPivots, M * k, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(maxDisSumPivots, M * k, MPI_INT, MPIMaxDisSumPivots, M * k, MPI_INT, 0, MPI_COMM_WORLD);
  if (world_rank == 0) {
    std::vector<u32> maxPtr(world_size), minPtr(world_size);
    for (int i = 0; i < M; ++i) {
      f64 max_value = -1 / 0.0, min_value = 1 / 0.0;
      i32 max_idx = -1, min_idx = -1;
      for (u32 j = 0; j < world_size; ++j) {
        if (MPIMaxDistanceSum[j * M + maxPtr[j]] > max_value) {
          max_value = MPIMaxDistanceSum[j * M + maxPtr[j]];
          max_idx = j;
        }
        if (MPIMinDistanceSum[j * M + minPtr[j]] < min_value) {
          min_value = MPIMinDistanceSum[j * M + minPtr[j]];
          min_idx = j;
        }
      }
      maxDistanceSum[i] = max_value;
      minDistanceSum[i] = min_value;
      for (int j = 0; j < k; ++j) {
        maxDisSumPivots[i * k + j] = MPIMaxDisSumPivots[max_idx * M * k + maxPtr[max_idx] * k + j]; // MPIMaxDisSumPivots[max_idx][maxPtr[max_idx]][j]
        minDisSumPivots[i * k + j] = MPIMinDisSumPivots[min_idx * M * k + minPtr[min_idx] * k + j];
      }
      ++maxPtr[max_idx];
      ++minPtr[min_idx];
    }
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(NULL, NULL);
  // 通过调用以下方法来得到所有可以工作的进程数量
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // 得到当前进程的秩
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // filename : input file namespace
  char *filename = (char *)"uniformvector-2dim-5h.txt";
  if (argc == 2) {
    filename = argv[1];
  } else if (argc != 1) {
    printf("Usage: ./pivot <filename>\n");
    return -1;
  }
  // M : number of combinations to store
  const int M = 1000;
  // dim : dimension of metric space
  int dim;
  // n : number of points
  int n;
  // k : number of pivots
  int k;

  // Read parameter
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    printf("%s file not found.\n", filename);
    return -1;
  }
  fscanf(file, "%d", &dim);
  fscanf(file, "%d", &n);
  fscanf(file, "%d", &k);
  printf("dim = %d, n = %d, k = %d\n", dim, n, k);

  // Start timing
  struct timeval start;

  // Read Data
  double *coord = (double *)malloc(sizeof(double) * dim * n);
  int i;
  for (i = 0; i < n; i++) {
    int j;
    for (j = 0; j < dim; j++) {
      fscanf(file, "%lf", &coord[i * dim + j]);
    }
  }
  fclose(file);
  gettimeofday(&start, NULL);

  // maxDistanceSum : the largest M distance sum
  double *maxDistanceSum = (double *)malloc(sizeof(double) * (M + 1));
  for (i = 0; i < M; i++) {
    maxDistanceSum[i] = 0;
  }
  // maxDisSumPivots : the top M pivots combinations
  int *maxDisSumPivots = (int *)malloc(sizeof(int) * k * (M + 1));
  for (i = 0; i < M; i++) {
    int ki;
    for (ki = 0; ki < k; ki++) {
      maxDisSumPivots[i * k + ki] = 0;
    }
  }
  // minDistanceSum : the smallest M distance sum
  double *minDistanceSum = (double *)malloc(sizeof(double) * (M + 1));
  for (i = 0; i < M; i++) {
    minDistanceSum[i] = __DBL_MAX__;
  }
  // minDisSumPivots : the bottom M pivots combinations
  int *minDisSumPivots = (int *)malloc(sizeof(int) * k * (M + 1));
  for (i = 0; i < M; i++) {
    int ki;
    for (ki = 0; ki < k; ki++) {
      minDisSumPivots[i * k + ki] = 0;
    }
  }

  // temp : indexes of pivots with dummy array head
  int *temp = (int *)malloc(sizeof(int) * (k + 1));
  temp[0] = -1;

  // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
  Combination(0, k, n, dim, M, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);
  if (world_rank == 0) {
    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    // Store the result
    FILE *out = fopen("result.txt", "w");
    for (i = 0; i < M; i++) {
      int ki;
      for (ki = 0; ki < k - 1; ki++) {
        fprintf(out, "%d ", maxDisSumPivots[i * k + ki]);
      }
      fprintf(out, "%d\n", maxDisSumPivots[i * k + k - 1]);
    }
    for (i = 0; i < M; i++) {
      int ki;
      for (ki = 0; ki < k - 1; ki++) {
        fprintf(out, "%d ", minDisSumPivots[i * k + ki]);
      }
      fprintf(out, "%d\n", minDisSumPivots[i * k + k - 1]);
    }
    fclose(out);

    // Log
    int ki;
    printf("max : ");
    for (ki = 0; ki < k; ki++) {
      printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for (ki = 0; ki < k; ki++) {
      printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);
    // for(i=0; i<M; i++){
    // int ki;
    // for(ki=0; ki<k; ki++){
    // printf("%d\t", maxDisSumPivots[i*k + ki]);
    // }
    // printf("%lf\n", maxDistanceSum[i]);
    // }
    // for(i=0; i<M; i++){
    // int ki;
    // for(ki=0; ki<k; ki++){
    // printf("%d\t", minDisSumPivots[i*k + ki]);
    // }
    // printf("%lf\n", minDistanceSum[i]);
    // }
  }
  MPI_Finalize();
  return 0;
}
