#include "defs.hpp"
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <math.h>
#include <pthread.h>
#include <iostream>
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

double distance(const double *coord, int ndims, int x, int y) {
  double dist = .0;
  for (int i = 0; i < ndims; i++) {
    dist += (coord[ndims * x + i] - coord[ndims * y + i]) * (coord[ndims * x + i] - coord[ndims * y + i]);
  }
  return sqrt(dist);
}

double calc_value(const int prev, const int npoints, const int npivots, const int ndims, int *pivots, const double *coord, double *rebuilt_coord, double *mx) {
  // Part 1. Rebuild Coordintate System
  for (int k = prev; k < npivots; k++) {
    int p = pivots[k];
    for (int i = 0; i < npoints; i++) {
      rebuilt_coord[k * npoints + i] = distance(coord, ndims, p, i);
    }
  }

  // Part 2. Evaluate System Value

  // Part 2.1. Calculate Chebyshev Distance
  for (int k = prev; k < npivots; k++) {
    for (int i = 0; i < npoints; i++) {
      for (int j = i + 1; j < npoints; j++) {
        double chebyshev_dim_dist = fabs(rebuilt_coord[k * npoints + i] - rebuilt_coord[k * npoints + j]);
        if (k > 0 && chebyshev_dim_dist < mx[(k - 1) * npoints * npoints + i * npoints + j]) {
          chebyshev_dim_dist = mx[(k - 1) * npoints * npoints + i * npoints + j];
        }
        mx[k * npoints * npoints + i * npoints + j] = chebyshev_dim_dist;
      }
    }
  }

  // Part 2.2. Calculate Sum of Chebyshev Distance between Each Pair

  double *chebyshev_dist = &mx[(npivots - 1) * npoints * npoints];
  double chebyshev_dist_sum = .0;
  for (int i = 0; i < npoints; i++) {
    for (int j = i + 1; j < npoints; j++) {
      chebyshev_dist_sum += chebyshev_dist[i * npoints + j];
    }
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
void combinations(const int start_point, const int end_point, const int npoints, const int npivots, const int ndims, const int M, const double *coord, MinMaxPivotPtrs *ptrs) {
  int *minDisSumPivots = (int *)malloc(sizeof(int) * M * npivots);
  int *maxDisSumPivots = (int *)malloc(sizeof(int) * M * npivots);
  double *minDistanceSum = (double *)malloc(sizeof(double) * M);
  double *maxDistanceSum = (double *)malloc(sizeof(double) * M);
  ;
  ptrs->minDisSumPivots = minDisSumPivots;
  ptrs->maxDisSumPivots = maxDisSumPivots;
  ptrs->minDistanceSum = minDistanceSum;
  ptrs->maxDistanceSum = maxDistanceSum;

  double *rebuilt_coord = (double *)malloc(sizeof(double) * npivots * npoints);
  double *mx = (double *)malloc(sizeof(double) * npivots * npoints * npoints);
  int *maxTmpPivots = (int *)malloc(sizeof(int) * M * npivots);
  int *minTmpPivots = (int *)malloc(sizeof(int) * M * npivots);
  std::map<double, int> mx_mp{};
  std::map<double, int> mn_mp{};
  int pivots[npivots];
  mth_comb(pivots, npoints, npivots, start_point);
  int prev = 0;
  for (int comb_cnt = start_point; comb_cnt < end_point && prev != -1; ++comb_cnt) {
    double value = calc_value(prev, npoints, npivots, ndims, pivots, coord, rebuilt_coord, mx);

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

  free(maxTmpPivots);
  free(minTmpPivots);
  free(rebuilt_coord);
  free(mx);
}

std::vector<int> affinity = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};

struct PThreadArg {
  int start_point;
  int end_point;
  int npoints;
  int npivots;
  int ndims;
  int M;
  const double *coord;
  MinMaxPivotPtrs *ptrs;
};

void *pthread_combination(void *vargs) {
  PThreadArg *args = (PThreadArg *)vargs;
  // auto t1 = std::chrono::high_resolution_clock::now();
  combinations(args->start_point, args->end_point, args->npoints, args->npivots, args->ndims, args->M, args->coord, args->ptrs);
  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << "f() took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us\n";
  return nullptr;
}

void Combination(int ki, const int k, const int n, const int dim, const int M, const double *coord, int *pivots, double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots) {
  u32 num_threads = 64;
  i32 cnk = choose(n, k);
  std::vector<MinMaxPivotPtrs> thread_data(num_threads);
  std::vector<pthread_t> threads(num_threads);
  std::vector<PThreadArg> args(num_threads);
  for (u32 i = 0; i < num_threads; ++i) {
    i32 start_point = (cnk / num_threads) * i;
    i32 end_point = start_point + cnk / num_threads;
    if (i == num_threads - 1) {
      end_point = cnk;
    }
    args[i].start_point = start_point;
    args[i].end_point = end_point;
    args[i].npoints = n;
    args[i].npivots = k;
    args[i].ndims = dim;
    args[i].M = M;
    args[i].coord = coord;
    args[i].ptrs = &thread_data[i];
    // threads[i] = std::thread([&, i, start_point, end_point, n, k, dim, M, coord] { combinations(start_point, end_point, n, k, dim, M, coord, &thread_data[i]); });
    // bind thread to core
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(affinity[i], &cpuset);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
    pthread_create(&threads[i], &attr, pthread_combination, &args[i]);
  }
  for (auto &i : threads) {
    pthread_join(i, NULL);
  }
  std::vector<u32> maxPtr(M), minPtr(M);
  for (int i = 0; i < M; ++i) {
    f64 max_value = -1 / 0.0, min_value = 1 / 0.0;
    i32 max_idx = -1, min_idx = -1;
    for (u32 j = 0; j < num_threads; ++j) {
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
}

int main(int argc, char *argv[]) {
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

  return 0;
}
