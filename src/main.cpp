#include "ac_float.h"
#include "ac_int.h"
#include "read_txt.h"
#include "top.h"
#include <iostream>

#define A_SIZE 5
#define H_COLS 500
#define W_COLS 16

template <typename T, int N> void create_random_graph(T A[N][N]) {
  srand(time(0));
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      int edge = rand() % 2;
      A[i][j] = edge;
      A[j][i] = edge;
    }
    A[i][i] = 1;
  }
}

template <typename T, typename I, int N>
void compute_sqrt_degree_matrix(T A[N][N], I D[N][N]) {
  for (int i = 0; i < N; ++i) {
    T degree = 0.0;
    for (int j = 0; j < N; ++j) {
      degree += A[i][j];
    }
    for (int j = 0; j < N; ++j) {
      D[i][j] = (i == j) ? 1 / sqrt(degree) : 0.0;
    }
  }
}

template <typename T, int N>
void compute_normalized_adjacency_matrix(T A[N][N], T D[N][N],
                                         T A_tilde[N][N]) {
  T tmp[N][N];
  for (int k = 0; k < A_SIZE; k++) {
    for (int m = 0; m < A_SIZE; m++) {
      for (int n = 0; n < A_SIZE; n++) {
        tmp[m][n] += D[m][k] * A[k][n];
      }
    }
  }

  for (int k = 0; k < A_SIZE; k++) {
    for (int m = 0; m < A_SIZE; m++) {
      for (int n = 0; n < A_SIZE; n++) {
        A_tilde[m][n] += tmp[m][k] * D[k][n];
      }
    }
  }
}

template <typename T, int A, int B> void create_feature_matrix(T H_in[A][B]) {
  srand(time(0));
  for (int i = 0; i < A; ++i) {
    for (int j = 0; j < B; ++j) {
      H_in[i][j] = abs(rand() % 256);
    }
  }
}

int main() {
  std::cout << "Init...." << std::endl;
  double A[A_SIZE][A_SIZE] = {{1, 1, 1, 0, 0},
                              {1, 1, 1, 1, 0},
                              {1, 1, 1, 1, 1},
                              {0, 1, 1, 1, 1},
                              {0, 0, 1, 1, 1}};
  double D[A_SIZE][A_SIZE];
  double A_tilde[A_SIZE][A_SIZE];
  double H_in[A_SIZE][H_COLS];

  std::cout << "Starting init data..." << std::endl;
  // create_random_graph(A);
  std::cout << "Graph (Adjacency Matrix):\n";
  for (int i = 0; i < A_SIZE; ++i) {
    for (int j = 0; j < A_SIZE; ++j) {
      std::cout << A[i][j] << " ";
    }
    std::cout << "\n";
  }
  compute_sqrt_degree_matrix(A, D);
  std::cout << "D (Adjacency Matrix):\n";
  for (int i = 0; i < A_SIZE; ++i) {
    for (int j = 0; j < A_SIZE; ++j) {
      std::cout << D[i][j] << " ";
    }
    std::cout << "\n";
  }
  compute_normalized_adjacency_matrix(A, D, A_tilde);
  std::cout << "Á_tilde (Adjacency Matrix):\n";
  for (int i = 0; i < A_SIZE; ++i) {
    for (int j = 0; j < A_SIZE; ++j) {
      std::cout << A_tilde[i][j] << " ";
    }
    std::cout << "\n";
  }
  create_feature_matrix<double, A_SIZE, H_COLS>(H_in);
  std::cout << "Finished init data" << std::endl;

  std::cout << "Starting calculating H_out..." << std::endl;
  double W[H_COLS][W_COLS];
  read_data<double, H_COLS, W_COLS>(W, "weights.txt");

  double temp[A_SIZE][H_COLS]; // A_tilde * H_in
  double H_out[A_SIZE][W_COLS];

  for (int k = 0; k < A_SIZE; k++) {
    for (int m = 0; m < A_SIZE; m++) {
      for (int n = 0; n < H_COLS; n++) {
        temp[m][n] += A_tilde[m][k] * H_in[k][n];
      }
    }
  }
  for (int k = 0; k < H_COLS; k++) {
    for (int m = 0; m < A_SIZE; m++) {
      for (int n = 0; n < H_COLS; n++) {
        H_out[m][n] += temp[m][k] * W[k][n];
      }
    }
  }
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < W_COLS; j++) {
      H_out[i][j] = (H_out[i][j] > 0) ? H_out[i][j] : 0;
    }
  }
  std::cout << "Finished calculating H_out" << std::endl;

  std::cout << "Starting calculating ac_H_out..." << std::endl;
  ac_fixed<16, 8, true> ac_A[A_SIZE][A_SIZE];
  ac_fixed<16, 8, true> ac_H_out[A_SIZE][W_COLS];
  ac_fixed<16, 8, true> ac_H_in[A_SIZE][H_COLS];

  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < A_SIZE; j++) {
      ac_A[i][j] = ac_fixed<16, 8, true>(A_tilde[i][j]);
      std::cout << A_tilde[i][j] << " " << ac_A[i][j] << std::endl;
    }
  }
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < H_COLS; j++) {
      ac_H_in[i][j] = ac_fixed<16, 8, true>(H_in[i][j]);
      std::cout << H_in[i][j] << " " << ac_H_in[i][j] << std::endl;
    }
  }

  calc(ac_A, ac_H_in, ac_H_out);
  std::cout << "Finished calculating ac_H_out" << std::endl;

  std::cout << "Starting testing..." << std::endl;
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < W_COLS; j++) {
      if (ac_H_out[i][j] != H_out[i][j]) {
        std::cout << "At [" << i << "," << j << "]" << std::endl;
        std::cout << "Expected: " << H_out[i][j] << std::endl;
        std::cout << "Got: " << ac_H_out[i][j] << std::endl;
        return -1;
      }
    }
  }

  std::cout << "Finished Testing!!!" << std::endl;

  return 0;
}
