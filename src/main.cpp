#include "ac_float.h"
#include "ac_int.h"
#include "read_txt.h"
#include "top.h"
#include <iomanip>
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
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      tmp[i][j] = 0;
      A_tilde[i][j] = 0;
    }
  }
  for (int k = 0; k < N; k++) {
    for (int m = 0; m < N; m++) {
      for (int n = 0; n < N; n++) {
        tmp[m][n] += D[m][k] * A[k][n];
      }
    }
  }
  for (int k = 0; k < N; k++) {
    for (int m = 0; m < N; m++) {
      for (int n = 0; n < N; n++) {
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
  std::ofstream out_file("feature_matrix.txt");
  if (!out_file) {
    std::cerr << "Error opening file!" << std::endl;
  }
  for (int i = 0; i < A; ++i) {
    for (int j = 0; j < B; ++j) {
      out_file << H_in[i][j];
      if (j != B - 1) {
        out_file << ", ";
      }
    }
    out_file << std::endl;
  }
  out_file.close();
  std::cout << "H_in written to feature_matrix.txt" << std::endl;
}

template <typename T>
void matmul(T *A, T *B, T *C, int A_rows, int A_cols, int B_cols) {
  for (int i = 0; i < A_rows; ++i) {
    for (int j = 0; j < B_cols; ++j) {
      C[i * B_cols + j] = 0;
      for (int k = 0; k < A_cols; ++k) {
        C[i * B_cols + j] += A[i * A_cols + k] * B[k * B_cols + j];
      }
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
  matmul(*A_tilde, *H_in, *temp, A_SIZE, A_SIZE, H_COLS);

  double H_out[A_SIZE][W_COLS];
  matmul(*temp, *W, *H_out, A_SIZE, H_COLS, W_COLS);
  // relu
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < W_COLS; j++) {
      H_out[i][j] = (H_out[i][j] > 0) ? H_out[i][j] : 0;
    }
  }
  std::cout << "Finished calculating H_out" << std::endl;

  std::cout << "Starting calculating ac_H_out..." << std::endl;
  ac_fixed<45, 9, true> ac_A[A_SIZE][A_SIZE] = {0};
  ac_fixed<45, 9, true> ac_H_out[A_SIZE][W_COLS] = {0};
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < W_COLS; j++) {
      ac_H_out[i][j] = 0;
    }
  }

  ac_fixed<45, 9, true> ac_H_in[A_SIZE][H_COLS] = {0};

  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < A_SIZE; j++) {
      ac_A[i][j] = ac_fixed<45, 9, true>(A_tilde[i][j]);
      // std::cout << A_tilde[i][j] << " " << ac_A[i][j] << std::endl;
    }
  }
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < H_COLS; j++) {
      ac_H_in[i][j] = ac_fixed<45, 9, true>(H_in[i][j]);
      // std::cout << H_in[i][j] << " " << ac_H_in[i][j] << std::endl;
    }
  }
  calc<ac_fixed<45, 9, true>, A_SIZE, H_COLS, W_COLS>(ac_A, ac_H_in, ac_H_out);
  std::cout << "Finished calculating ac_H_out" << std::endl;

  std::cout << "Starting testing..." << std::endl;
  // for (int i = 0; i < A_SIZE; i++) {
  //   for (int j = 0; j < W_COLS; j++) {
  //     if (ac_H_out[i][j] != H_out[i][j]) {
  //       std::cout << "At [" << i << "," << j << "]" << std::endl;
  //       std::cout << "Expected: " << H_out[i][j] << std::endl;
  //       std::cout << "Got: " << ac_H_out[i][j] << std::endl;
  //       // return -1;
  //     }
  //   }
  // }

  std::cout << "H_out:" << std::endl;
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < W_COLS; ++j) {
      std::cout << std::setw(4) << std::scientific << std::setprecision(2)
                << H_out[i][j] << " ";
    }
    std::cout << "\n";
  }

  std::cout << "ac_H_out:" << std::endl;
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < W_COLS; ++j) {
      std::cout << std::setw(4) << std::scientific << std::setprecision(2)
                << ac_H_out[i][j].to_double() << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Finished Testing!!!" << std::endl;

  return 0;
}
