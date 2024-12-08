#include "read_txt.h"

#pragma HLS target_frequency "500 MHz"

// template <typename T>
// void matmul(T *A, T *B, T *C, int A_rows, int A_cols, int B_cols) {
//   for (int i = 0; i < A_rows; ++i) {
//     for (int j = 0; j < B_cols; ++j) {
//       C[i * B_cols + j] = 0;
//       for (int k = 0; k < A_cols; ++k) {
//         C[i * B_cols + j] += A[i * A_cols + k] * B[k * B_cols + j];
//       }
//     }
//   }
// }

template <typename T, int M, int K, int N, int ID>
void multiply_matrices(T A[M][K], T B[K][N], T C[M][N]) {
  for (int k = 0; k < K; k++) {
#pragma HLS block name = MUL_MATRIX
#pragma hls_unroll 8
    for (int m = 0; m < M; m++) {
      for (int n = 0; n < N; n++) {
        C[m][n] += A[m][k] * B[k][n];
      }
    }
  }
}

template <typename T, int M, int N> void relu(T in[M][N], T out[M][N]) {
#pragma HLS block name = RELU
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      out[m][n] = (in[m][n] > 0) ? in[m][n] : T(0);
      // ac_math::ac_relu(in[m][n], out[m][n])
    }
  }
}

template <typename T, int A_SIZE, int H_COLS, int W_COLS>
void calc(T A[A_SIZE][A_SIZE], T H_in[A_SIZE][H_COLS],
          T H_out[A_SIZE][W_COLS]) {
#pragma HLS top name = CALC
  T weights[H_COLS][W_COLS];
  read_data<T, H_COLS, W_COLS>(weights, "weights.txt");

  T tmp1[A_SIZE][H_COLS];
  // for (int i = 0; i < A_SIZE; i++) {
  //   for (int j = 0; j < H_COLS; j++) {
  //     tmp1[i][j] = T(0);
  //   }
  // }

  multiply_matrices<T, A_SIZE, A_SIZE, H_COLS, 0>(A, H_in, tmp1);

  T tmp2[A_SIZE][W_COLS];
  // for (int i = 0; i < A_SIZE; i++) {
  //   for (int j = 0; j < W_COLS; j++) {
  //     tmp2[i][j] = T(0);
  //   }
  // }

  multiply_matrices<T, A_SIZE, H_COLS, W_COLS, 1>(tmp1, weights, tmp2);

  relu<T, A_SIZE, W_COLS>(tmp2, H_out);
}
