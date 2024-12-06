#include "read_txt.h"

#pragma HLS target_frequency "500 MHz"

template <typename T, int M, int K, int N, int ID>
void multiply_matrices(T A[M][K], T B[K][N], T C[M][N]) {
  for (int k = 0; k < K; k++) {
#pragma hls_unroll
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        C[m][n] += A[m][k] * B[k][n];
      }
    }
  }
}

template <typename T, int M, int N> void relu(T in[M][N], T out[M][N]) {
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      out[m][n] = (in[m][n] > 0) ? in[m][n] : 0;
    }
  }
}

template <typename T, int A_SIZE, int H_COLS, int W_COLS>
void calc(T A[A_SIZE][A_SIZE], T H_in[A_SIZE][H_COLS],
          T H_out[A_SIZE][W_COLS]) {
  T weights[H_COLS][W_COLS];
  read_data<T, H_COLS, W_COLS>(weights, "weights.txt");

  T tmp1[A_SIZE][H_COLS];
  multiply_matrices<T, A_SIZE, A_SIZE, H_COLS, 0>(A, H_in, tmp1);

  T tmp2[A_SIZE][W_COLS];
  multiply_matrices<T, A_SIZE, H_COLS, W_COLS, 1>(tmp1, weights, tmp2);

  relu<T, A_SIZE, W_COLS>(tmp2, H_out);
}
