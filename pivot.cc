#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

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

void kth_comb(int *arr, int n, int k) {}

double distance(double *coord, int ndims, int x, int y) {
  double dist = .0;
  for (int i = 0; i < ndims; i++) {
    dist += (coord[ndims * x + i] - coord[ndims * y + i]) * (coord[ndims * x + i] - coord[ndims * y + i]);
  }
  return sqrt(dist);
}

double calc_value(const int prev, const int npoints, const int npivots, const int ndims, int *pivots, double *coord, double *rebuilt_coord, double *mx) {
  // Part 1. Rebuild Coordintate System
  for (int k = prev; k < npivots; k++) {
    int p = pivots[k];
    for (int i = 0; i < npoints; i++) {
      rebuilt_coord[k * npoints + i] = distance(coord, ndims, p, i);
    }
  }

  // Part 2. Evaluate System Value
  const int point_pairs = npoints * (npoints - 1) / 2;

  // Part 2.1. Calculate Chebyshev Distance
  for (int k = prev; k < npivots; k++) {
    int idx_cnt = 0;
    for (int i = 0; i < npoints; i++) {
      for (int j = 0; j < i; j++) {
        double chebyshev_dim_dist = fabs(rebuilt_coord[k * npoints + i] - rebuilt_coord[k * npoints + j]);
        if (k > 0 && chebyshev_dim_dist < mx[(k - 1) * point_pairs + idx_cnt + j]) {
          chebyshev_dim_dist = mx[(k - 1) * point_pairs + idx_cnt + j];
        }
        mx[k * point_pairs + idx_cnt + j] = chebyshev_dim_dist;
      }
      idx_cnt += i;
    }
  }

  // Part 2.2. Calculate Sum of Chebyshev Distance between Each Pair

  double *chebyshev_dist = &mx[(npivots - 1) * point_pairs];
  double chebyshev_dist_sum = .0;
  int idx_cnt = 0;
  for (int i = 0; i < npoints; i++) {
    for (int j = 0; j < i; j++) {
      chebyshev_dist_sum += chebyshev_dist[idx_cnt + j];
    }
    idx_cnt += i;
  }

  // Calculate Half of All Pairs, Then Double
  return chebyshev_dist_sum * 2;
}

// maxDisSum, minDisSum, maxDisSumPivots, minDisSumPivots
void combinations(const int npoints, const int npivots, const int ndims, const int M, double *coord, int *pivots, double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots) {
  const int point_pairs = npoints * (npoints - 1) / 2;
  double *rebuilt_coord = (double *)malloc(sizeof(double) * npivots * npoints);
  double *mx = (double *)malloc(sizeof(double) * npivots * point_pairs);
  int *maxTmpPivots = (int *)malloc(sizeof(int) * M * npivots);
  int *minTmpPivots = (int *)malloc(sizeof(int) * M * npivots);
  std::map<double, int> mx_mp{};
  std::map<double, int> mn_mp{};
  for (int i = 0; i < npivots; i++) {
    pivots[i] = i;
  }
  int prev = 0;
  while (prev != -1) {
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

void Combination(int ki, const int k, const int n, const int dim, const int M, double *coord, int *pivots, double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots) { combinations(n, k, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots); }

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
