#include <benchmark/benchmark.h>
#include <blasfeo.h>
#include <stdio.h>
#include <stdlib.h>

#include "bdldl.h"

static void BM_chole(benchmark::State &state) {
  int ii;
  int n = 16;
  srand(100);

  double *A;
  d_zeros(&A, n, n);
  for (ii = 0; ii < n * n; ii++) A[ii] = rand();
  // d_print_mat(n, n, A, n);

  double *B;
  d_zeros(&B, n, n);
  double *D;
  d_zeros(&D, n, n);

  // work space enough for 6 matrix structs for size n times n
  int size_strmat = 6 * blasfeo_memsize_dmat(n, n);
  void *memory_strmat;
  v_zeros_align(&memory_strmat, size_strmat);
  char *ptr_memory_strmat = (char *)memory_strmat;

  struct blasfeo_dmat sA;
  blasfeo_create_dmat(n, n, &sA, ptr_memory_strmat);
  ptr_memory_strmat += sA.memsize;
  blasfeo_pack_dmat(n, n, A, n, &sA, 0, 0);
  // printf("\nA = \n");
  // blasfeo_print_dmat(n, n, &sA, 0, 0);

  struct blasfeo_dmat sB;
  blasfeo_create_dmat(n, n, &sB, ptr_memory_strmat);
  ptr_memory_strmat += sB.memsize;
  blasfeo_pack_dmat(n, n, B, n, &sB, 0, 0);
  // printf("\nB = \n");
  // blasfeo_print_dmat(n, n, &sB, 0, 0);

  struct blasfeo_dmat sL;
  blasfeo_create_dmat(n, n, &sL, ptr_memory_strmat);
  ptr_memory_strmat += sL.memsize;

  struct blasfeo_dmat sD;
  blasfeo_create_dmat(n, n, &sD, ptr_memory_strmat);
  ptr_memory_strmat += sD.memsize;

  blasfeo_dgemm_nt(n, n, n, 1.0, &sA, 0, 0, &sA, 0, 0, 0, &sB, 0, 0, &sD, 0, 0);
  // printf("\nA*A' = \n");
  // blasfeo_print_dmat(n, n, &sD, 0, 0);

  for (auto _ : state) {
    blasfeo_dpotrf_l(n, &sD, 0, 0, &sL, 0, 0);
    // printf("\nL = \n");
    // blasfeo_print_dmat(n, n, &sL, 0, 0);
  }

  d_free(A);
  d_free(B);
  d_free(D);
  v_free_align(memory_strmat);
}

static void BM_bdldl(benchmark::State &state) {
  initOCPData();
  for (auto _ : state) {
    factorizeOCPData();
  }
  freeOCPData();
}

// BENCHMARK(BM_chole);
BENCHMARK(BM_bdldl);
