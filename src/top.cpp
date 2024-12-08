#include "read_txt.h"

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
    for (int n = 0; n < N; n++) {
      out[m][n] = (in[m][n] > 0) ? in[m][n] : T(0);
    }
  }
}

#pragma hls_design top
template <typename INPUT, typename OUTPUT, int A_SIZE, int H_COLS, int W_COLS>
void calc(INPUT A[A_SIZE][A_SIZE], INPUT H_in[A_SIZE][H_COLS],
          INPUT weights[H_COLS][W_COLS], OUTPUT H_out[A_SIZE][W_COLS]) {
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
