#include <blasfeo.h>
#include <stdio.h>

#define M 2

struct blasfeo_dmat Im;

void inverse(struct blasfeo_dmat* A, struct blasfeo_dmat* B) {
  struct blasfeo_dmat L;
  blasfeo_allocate_dmat(M, M, &L);

  blasfeo_dpotrf_l(M, A, 0, 0, &L, 0, 0);
  blasfeo_dtrsm_rltn(M, M, 1.0, &L, 0, 0, &Im, 0, 0, B, 0, 0);
  blasfeo_print_dmat(M, M, B, 0, 0);

  blasfeo_dtrsm_rlnn(M, M, 1.0, &L, 0, 0, B, 0, 0, B, 0, 0);
  blasfeo_print_dmat(M, M, B, 0, 0);

  blasfeo_free_dmat(&L);
}

int main() {
  blasfeo_allocate_dmat(M, M, &Im);
  blasfeo_dgese(M, M, 0, &Im, 0, 0);
  blasfeo_ddiare(M, 1.0, &Im, 0, 0);

  double* A;
  d_zeros(&A, M, M);
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < M; i++) {
      A[i + j * M] = i + j * M + 1;
    }
  }
  struct blasfeo_dmat sA;
  blasfeo_allocate_dmat(M, M, &sA);
  blasfeo_pack_l_dmat(M, M, A, M, &sA, 0, 0);
  blasfeo_print_dmat(M, M, &sA, 0, 0);

  struct blasfeo_dmat sB;
  blasfeo_allocate_dmat(M, M, &sB);
  blasfeo_dgese(M, M, 0, &sB, 0, 0);
  blasfeo_dsyrk_ln(M, M, 1.0, &sA, 0, 0, &sA, 0, 0, 0, &sB, 0, 0, &sB, 0, 0);
  blasfeo_dtrtr_l(M, &sB, 0, 0, &sB, 0, 0);
  blasfeo_print_dmat(M, M, &sB, 0, 0);

  struct blasfeo_dmat sBInv;
  blasfeo_allocate_dmat(M, M, &sBInv);
  inverse(&sB, &sBInv);

  printf("inverse\n");
  blasfeo_print_dmat(M, M, &sBInv, 0, 0);

  struct blasfeo_dvec vb;
  blasfeo_allocate_dvec(M, &vb);
  blasfeo_dvecse(M, 1.0, &vb, 0);
  blasfeo_dgemv_n(M, M, 1.0, &sB, 0, 0, &vb, 0, 0.0, &vb, 0, &vb, 0);
  blasfeo_print_dvec(M, &vb, 0);

  struct blasfeo_dmat sL;
  blasfeo_allocate_dmat(M, M, &sL);
  blasfeo_dgese(M, M, 0, &sL, 0, 0);

  blasfeo_dpotrf_l(M, &sB, 0, 0, &sL, 0, 0);
  blasfeo_print_dmat(M, M, &sL, 0, 0);

  blasfeo_dtrsv_lnn(M, &sL, 0, 0, &vb, 0, &vb, 0);  // inv(L)b
  blasfeo_print_dvec(M, &vb, 0);

  blasfeo_dtrsv_ltn(M, &sL, 0, 0, &vb, 0, &vb, 0);  // inv(L')b
  blasfeo_print_dvec(M, &vb, 0);

  // blasfeo_free_dmat(&sA);
  // blasfeo_free_dmat(&sB);
  // blasfeo_free_dmat(&sL);

  return 0;
}