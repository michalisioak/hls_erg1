#include "mc_scverify.h"
#include "read_txt.h"
#include <ac_fixed.h>
#include <iostream>

#define C_A_SIZE 5
#define C_H_COLS 500
#define C_W_COLS 16

#define AC_INPUT ac_fixed<16, 9, true, AC_TRN, AC_SAT>
#define AC_OUTPUT ac_fixed<64, 32, true, AC_TRN, AC_SAT>

#define PERCISION 0.25e5

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

// #pragma hls_design block
template <typename T, int M, int K, int N, int ID>
void multiply_matrices(T A[M][K], T B[K][N], T C[M][N]) {
  for (int k = 0; k < K; k++) {
#pragma hls_unroll 8
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        C[m][n] += A[m][k] * B[k][n];
      }
    }
  }
}

// #pragma hls_design block
template <typename T, int M, int N> void relu(T in[M][N], T out[M][N]) {
  for (int m = 0; m < M; m++) {
#pragma hls_unroll 8
    for (int n = 0; n < N; n++) {
      out[m][n] = (in[m][n] > 0) ? in[m][n] : T(0);
    }
  }
}

#pragma hls_design top
template <typename INPUT, typename OUTPUT, int A_SIZE, int H_COLS, int W_COLS>
void CCS_BLOCK(calc)(INPUT A[A_SIZE][A_SIZE], INPUT H_in[A_SIZE][H_COLS],
                     INPUT weights[H_COLS][W_COLS],
                     OUTPUT H_out[A_SIZE][W_COLS]) {
  OUTPUT A_[A_SIZE][A_SIZE];
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < A_SIZE; j++) {
      A_[i][j] = OUTPUT(A[i][j]);
    }
  }
  OUTPUT H_in_[A_SIZE][H_COLS];
  for (int i = 0; i < A_SIZE; i++) {
    for (int j = 0; j < H_COLS; j++) {
      H_in_[i][j] = OUTPUT(H_in[i][j]);
    }
  }
  OUTPUT weights_[H_COLS][W_COLS];
  for (int i = 0; i < H_COLS; i++) {
    for (int j = 0; j < W_COLS; j++) {
      weights_[i][j] = OUTPUT(weights[i][j]);
    }
  }

  OUTPUT tmp1[A_SIZE][H_COLS];
  multiply_matrices<OUTPUT, A_SIZE, A_SIZE, H_COLS, 0>(A_, H_in_, tmp1);

  OUTPUT tmp2[A_SIZE][W_COLS];
  multiply_matrices<OUTPUT, A_SIZE, H_COLS, W_COLS, 1>(tmp1, weights_, tmp2);

  relu<OUTPUT, A_SIZE, W_COLS>(tmp2, H_out);
}
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

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

CCS_MAIN(int argc, char *argv) {
  std::cout << "Init...." << std::endl;
  double A[C_A_SIZE][C_A_SIZE] = {{1, 1, 1, 0, 0},
                                  {1, 1, 1, 1, 0},
                                  {1, 1, 1, 1, 1},
                                  {0, 1, 1, 1, 1},
                                  {0, 0, 1, 1, 1}};
  double D[C_A_SIZE][C_A_SIZE];
  double A_tilde[C_A_SIZE][C_A_SIZE];
  double H_in[C_A_SIZE][C_H_COLS];

  std::cout << "Starting init data..." << std::endl;
  // create_random_graph(A);
  std::cout << "Graph (Adjacency Matrix):\n";
  for (int i = 0; i < C_A_SIZE; ++i) {
    for (int j = 0; j < C_A_SIZE; ++j) {
      std::cout << A[i][j] << " ";
    }
    std::cout << "\n";
  }
  compute_sqrt_degree_matrix(A, D);
  std::cout << "D (Adjacency Matrix):\n";
  for (int i = 0; i < C_A_SIZE; ++i) {
    for (int j = 0; j < C_A_SIZE; ++j) {
      std::cout << D[i][j] << " ";
    }
    std::cout << "\n";
  }
  compute_normalized_adjacency_matrix(A, D, A_tilde);
  std::cout << "Á_tilde (Adjacency Matrix):\n";
  for (int i = 0; i < C_A_SIZE; ++i) {
    for (int j = 0; j < C_A_SIZE; ++j) {
      std::cout << A_tilde[i][j] << " ";
    }
    std::cout << "\n";
  }
  create_feature_matrix<double, C_A_SIZE, C_H_COLS>(H_in);

  double weights[C_H_COLS][C_W_COLS];
  read_data<double, C_H_COLS, C_W_COLS>(weights, "weights.txt");
  std::cout << "Finished init data" << std::endl;

  std::cout << "Starting calculating H_out..." << std::endl;

  double temp[C_A_SIZE][C_H_COLS]; // A_tilde * H_in
  matmul(*A_tilde, *H_in, *temp, C_A_SIZE, C_A_SIZE, C_H_COLS);

  double H_out[C_A_SIZE][C_W_COLS];
  matmul(*temp, *weights, *H_out, C_A_SIZE, C_H_COLS, C_W_COLS);
  // relu
  for (int i = 0; i < C_A_SIZE; i++) {
    for (int j = 0; j < C_W_COLS; j++) {
      H_out[i][j] = (H_out[i][j] > 0) ? H_out[i][j] : 0;
    }
  }
  std::cout << "Finished calculating H_out" << std::endl;

  std::cout << "Starting calculating ac_H_out..." << std::endl;
  AC_INPUT ac_A[C_A_SIZE][C_A_SIZE];
  AC_OUTPUT ac_H_out[C_A_SIZE][C_W_COLS] = {0};
  AC_INPUT ac_weights[C_H_COLS][C_W_COLS];

  AC_INPUT ac_H_in[C_A_SIZE][C_H_COLS] = {0};
  for (int i = 0; i < C_H_COLS; i++) {
    for (int j = 0; j < C_W_COLS; j++) {
      ac_weights[i][j] = AC_INPUT(weights[i][j]);
      // std::cout << H_in[i][j] << " " << ac_H_in[i][j] << std::endl;
    }
  }
  for (int i = 0; i < C_A_SIZE; i++) {
    for (int j = 0; j < C_A_SIZE; j++) {
      ac_A[i][j] = AC_INPUT(A_tilde[i][j]);
      // std::cout << A_tilde[i][j] << " " << ac_A[i][j] << std::endl;
    }
  }
  for (int i = 0; i < C_A_SIZE; i++) {
    for (int j = 0; j < C_H_COLS; j++) {
      ac_H_in[i][j] = AC_INPUT(H_in[i][j]);
      // std::cout << H_in[i][j] << " " << ac_H_in[i][j] << std::endl;
    }
  }
  calc<AC_INPUT, AC_OUTPUT, C_A_SIZE, C_H_COLS, C_W_COLS>(ac_A, ac_H_in,
                                                          ac_weights, ac_H_out);
  std::cout << "Finished calculating ac_H_out" << std::endl;

  std::cout << "Starting testing..." << std::endl;
  for (int i = 0; i < C_A_SIZE; i++) {
    for (int j = 0; j < C_W_COLS; j++) {
      if (abs(ac_H_out[i][j].to_double() - H_out[i][j]) > PERCISION) {
        std::cout << "At [" << i << "," << j << "]" << std::endl;
        std::cout << "Expected: " << H_out[i][j] << std::endl;
        std::cout << "Got: " << ac_H_out[i][j] << std::endl;
        CSS_RETURN(-1);
      }
    }
  }

  std::cout << "Finished Testing!!!" << std::endl;

  CSS_RETURN(0);
}
