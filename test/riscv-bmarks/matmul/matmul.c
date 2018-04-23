// See LICENSE for license details.

#include "dataset.h"
#include "util.h"
#include <stddef.h>

void matmul(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k;
  
  for (i = 0; i < lda; i++) {
    for (j = coreid; j < lda; j += ncores) {
      data_t sum = 0;
      for (k = 0; k < lda; k++)
        sum += A[j*lda + k] * B[k*lda + i];
      C[i + j*lda] = sum;
    }
  }
}

void matmul_opt(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  //block the code 
  size_t i, j, k;
  int start, end;

  start = coreid * lda / ncores;
  if(coreid == ncores - 1)
    end = lda;
  else
    //only 2 cores
    end = lda / ncores;

  // transpose b for spatial locality 
  size_t size = lda*lda;
  data_t tpd [size];
  for (size_t l = 0; l < size; l++) {
    int r = l / lda;
    int c = l % lda;
    tpd[c*lda + r] = B[r*lda + c];
  }
 
  for (j = start; j < end; j++) {
    for (i = 0; i < lda; i++) {
      data_t sum = 0;
      for (k = 0; k < lda; k += 8) {
        sum += A[j*lda + k] * tpd[i*lda + k];
        sum += A[j*lda + k + 1] * tpd[i*lda + k + 1];
        sum += A[j*lda + k + 2] * tpd[i*lda + k + 2];
        sum += A[j*lda + k + 3] * tpd[i*lda + k + 3];
        sum += A[j*lda + k + 4] * tpd[i*lda + k + 4];
        sum += A[j*lda + k + 5] * tpd[i*lda + k + 5];
        sum += A[j*lda + k + 6] * tpd[i*lda + k + 6];
        sum += A[j*lda + k + 7] * tpd[i*lda + k + 7];
      }
      C[i + j*lda] = sum;
    }
  }
}
