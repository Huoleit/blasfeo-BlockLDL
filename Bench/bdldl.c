#include <stdio.h>
#include <stdlib.h>

#include "subroutines.h"

#define NX 8
#define NU 4
#define NG 1
#define N 50

int nx[N + 1], nu[N + 1], ng[N + 1];
ldl_float A[NX * NX] = {1,   0, 0,   0, 0, 0, 0, 0, 0, 1,   0, 0,   0, 0, 0, 0,
                        0,   0, 1,   0, 0, 0, 0, 0, 0, 0,   0, 1,   0, 0, 0, 0,
                        0.1, 0, 0,   0, 1, 0, 0, 0, 0, 0.1, 0, 0,   0, 1, 0, 0,
                        0,   0, 0.1, 0, 0, 0, 1, 0, 0, 0,   0, 0.1, 0, 0, 0, 1};
ldl_float B[NX * NU] = {0.005, 0, 0,   0, 0.1, 0,     0, 0,     0, 0.005, 0,
                        0,     0, 0.1, 0, 0,   0,     0, 0.005, 0, 0,     0,
                        0.1,   0, 0,   0, 0,   0.005, 0, 0,     0, 0.1};
ldl_float Q[NX * NX] = {1, 0, 0, 0, 0,   0, 0,   0, 0, 1, 0, 0, 0, 0,   0, 0,
                        0, 0, 1, 0, 0,   0, 0,   0, 0, 0, 0, 1, 0, 0,   0, 0,
                        0, 0, 0, 0, 0.1, 0, 0,   0, 0, 0, 0, 0, 0, 0.1, 0, 0,
                        0, 0, 0, 0, 0,   0, 0.1, 0, 0, 0, 0, 0, 0, 0,   0, 0.1};
ldl_float R[NU * NU] = {0.3, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.3};

ldl_matrix AA[N], BB[N], QQ[N + 1], RR[N];

ldl_matrix LL[3 * N - 1];
ldl_matrix DD[3 * N];
ldl_matrix DDInv[3 * N];

ldl_matrix Ix, Iu, tmpLx, tmpLu;

void initOCPData() {
  int ii;

  for (ii = 0; ii < N + 1; ++ii) {
    nx[ii] = NX;
    nu[ii] = NU;
    ng[ii] = NG;
  }

  // for (ii = 0; ii < NX * NX; ++ii) A[ii] = ii % NX;
  // for (ii = 0; ii < NX * NU; ++ii) B[ii] = ii / NX;
  // for (ii = 0; ii < NX; ++ii) Q[(NX + 1) * ii] = 1.0;
  // for (ii = 0; ii < NU; ++ii) R[(NU + 1) * ii] = 2.0;

  ALLOCATE_MAT(NX, NX, &Ix);
  GESE(NX, NX, 0, &Ix, 0, 0);
  DIARE(NX, 1.0, &Ix, 0, 0);

  ALLOCATE_MAT(NU, NU, &Iu);
  GESE(NU, NU, 0, &Iu, 0, 0);
  DIARE(NU, 1.0, &Iu, 0, 0);

  ALLOCATE_MAT(NX, NX, &tmpLx);
  ALLOCATE_MAT(NU, NU, &tmpLu);

  for (ii = 0; ii < N; ++ii) {
    ALLOCATE_MAT(nx[ii + 1], nx[ii], AA + ii);
    PACK_MAT(nx[ii + 1], nx[ii], A, nx[ii + 1], AA + ii, 0, 0);  // A

    ALLOCATE_MAT(nx[ii + 1], nu[ii], BB + ii);
    PACK_MAT(nx[ii + 1], nu[ii], B, nx[ii + 1], BB + ii, 0, 0);  // B

    ALLOCATE_MAT(nx[ii], nx[ii], QQ + ii);
    PACK_MAT(nx[ii], nx[ii], Q, nx[ii], QQ + ii, 0, 0);  // Q

    ALLOCATE_MAT(nu[ii], nu[ii], RR + ii);
    PACK_MAT(nu[ii], nu[ii], R, nu[ii], RR + ii, 0, 0);  // R

    if (ii != 0) ALLOCATE_MAT(nx[ii + 1], nx[ii], LL + ii * 3 - 1);  // L_A
    ALLOCATE_MAT(nx[ii + 1], nu[ii], LL + ii * 3);                   // L_B
    ALLOCATE_MAT(nx[ii + 1], nx[ii + 1], LL + ii * 3 + 1);           // L_I

    ALLOCATE_MAT(nu[ii], nu[ii], DD + ii * 3);              // D_R
    ALLOCATE_MAT(nx[ii + 1], nx[ii + 1], DD + ii * 3 + 1);  // D_B
    ALLOCATE_MAT(nx[ii + 1], nx[ii + 1], DD + ii * 3 + 2);  // D_Q

    ALLOCATE_MAT(nu[ii], nu[ii], DDInv + ii * 3);              // DInv_R
    ALLOCATE_MAT(nx[ii + 1], nx[ii + 1], DDInv + ii * 3 + 1);  // DInv_B
    ALLOCATE_MAT(nx[ii + 1], nx[ii + 1], DDInv + ii * 3 + 2);  // DInv_Q
  }
  ALLOCATE_MAT(nx[ii], nx[ii], QQ + ii);
  PACK_MAT(nx[ii], nx[ii], Q, nx[ii], QQ + ii, 0, 0);  // Q_N
}

void freeOCPData() {
  int ii;

  for (ii = 0; ii < N; ++ii) {
    FREE_MAT(AA + ii);
    FREE_MAT(BB + ii);
    FREE_MAT(QQ + ii);
    FREE_MAT(RR + ii);
  }
  FREE_MAT(QQ + ii);

  for (ii = 0; ii < 3 * N - 1; ++ii) {
    FREE_MAT(LL + ii);
    FREE_MAT(DD + ii);
    FREE_MAT(DDInv + ii);
  }
  FREE_MAT(DD + ii);
  FREE_MAT(DDInv + ii);

  FREE_MAT(&Ix);
  FREE_MAT(&Iu);

  FREE_MAT(&tmpLx);
  FREE_MAT(&tmpLu);
}

void printWorkspace() {
  int ii;

  for (ii = 0; ii < N; ii++) {
    printf("%d  D_R = \n", ii);
    PRINT_MAT(nu[ii], nu[ii], DD + ii * 3, 0, 0);
    printf("%d  D_B = \n", ii);
    PRINT_MAT(nx[ii + 1], nx[ii + 1], DD + ii * 3 + 1, 0, 0);
    printf("%d  D_Q = \n", ii);
    PRINT_MAT(nx[ii + 1], nx[ii + 1], DD + ii * 3 + 2, 0, 0);
  }
  for (ii = 0; ii < N; ii++) {
    printf("%d  DInv_R = \n", ii);
    PRINT_MAT(nu[ii], nu[ii], DDInv + ii * 3, 0, 0);
    printf("%d  DInv_B = \n", ii);
    PRINT_MAT(nx[ii + 1], nx[ii + 1], DDInv + ii * 3 + 1, 0, 0);
    printf("%d  DInv_Q = \n", ii);
    PRINT_MAT(nx[ii + 1], nx[ii + 1], DDInv + ii * 3 + 2, 0, 0);
  }
  for (ii = 0; ii < N; ii++) {
    if (ii != 0) {
      printf("%d  L_A = \n", ii);
      PRINT_MAT(nx[ii + 1], nx[ii], LL + ii * 3 - 1, 0, 0);
    }
    printf("%d  L_B = \n", ii);
    PRINT_MAT(nx[ii + 1], nu[ii], LL + ii * 3, 0, 0);
    printf("%d  L_I = \n", ii);
    PRINT_MAT(nx[ii + 1], nx[ii + 1], LL + ii * 3 + 1, 0, 0);
  }
}

void printOCPData() {
  int ii;

  for (ii = 0; ii < N; ii++) {
    printf("\nk = %d nx = %d nu = %d ng = %d\n", ii, nx[ii], nu[ii], ng[ii]);
    printf("A = \n");
    PRINT_MAT(nx[ii + 1], nx[ii], AA + ii, 0, 0);
    printf("B = \n");
    PRINT_MAT(nx[ii + 1], nu[ii], BB + ii, 0, 0);

    printf("Q = \n");
    PRINT_MAT(nx[ii], nx[ii], QQ + ii, 0, 0);
    printf("R = \n");
    PRINT_MAT(nu[ii], nu[ii], RR + ii, 0, 0);
  }
  printf("\nk = %d nx = %d\n", ii, nx[ii]);
  printf("Q = \n");
  PRINT_MAT(nx[ii], nx[ii], QQ + ii, 0, 0);

  printf("Ix = \n");
  PRINT_MAT(NX, NX, &Ix, 0, 0);
  printf("Iu = \n");
  PRINT_MAT(NU, NU, &Iu, 0, 0);
}

// LLTInv <= inv(L)
void LLTInvNx(ldl_matrix* L, ldl_matrix* LLTInv) {
  TRSM_LLNN(NX, NX, 1.0, L, 0, 0, &Ix, 0, 0, LLTInv, 0, 0);
  TRSM_LLTN(NX, NX, 1.0, L, 0, 0, LLTInv, 0, 0, LLTInv, 0, 0);
}

void LLTInvNu(ldl_matrix* L, ldl_matrix* LLTInv) {
  TRSM_LLNN(NU, NU, 1.0, L, 0, 0, &Iu, 0, 0, LLTInv, 0, 0);
  TRSM_LLTN(NU, NU, 1.0, L, 0, 0, LLTInv, 0, 0, LLTInv, 0, 0);
}

void factorizeOCPData() {
  for (int k = 0; k < N; ++k) {
    // D = R
    GECP(nu[k], nu[k], RR + k, 0, 0, DD + k * 3, 0, 0);
    // DInv = Cholesky(D)
    POTRF_L(nu[k], DD + k * 3, 0, 0, &tmpLu, 0, 0);
    LLTInvNu(&tmpLu, DDInv + k * 3);

    if (k != 0) {
      // L = A*DInv
      GEMM_NT(nx[k + 1], nx[k], nx[k], 1.0, AA + k, 0, 0, DDInv + k * 3 - 1, 0,
              0, 0.0, LL + k * 3 - 1, 0, 0, LL + k * 3 - 1, 0, 0);
    }
    // L = B*DInv
    GEMM_NT(nx[k + 1], nu[k], nu[k], 1.0, BB + k, 0, 0, DDInv + k * 3, 0, 0,
            0.0, LL + k * 3, 0, 0, LL + k * 3, 0, 0);

    // D = B * L'
    GEMM_NT(nx[k + 1], nx[k + 1], nu[k], 1.0, BB + k, 0, 0, LL + k * 3, 0, 0,
            0.0, DD + k * 3 + 1, 0, 0, DD + k * 3 + 1, 0, 0);
    if (k != 0) {
      GEMM_NT(nx[k + 1], nx[k + 1], nx[k], 1.0, AA + k, 0, 0, LL + k * 3 - 1, 0,
              0, 1.0, DD + k * 3 + 1, 0, 0, DD + k * 3 + 1, 0, 0);
    }
    DIARE(nx[k + 1], 1e-5, DD + k * 3 + 1, 0, 0);  // regularization
    // DInv = Cholesky(D)
    POTRF_L(nx[k + 1], DD + k * 3 + 1, 0, 0, &tmpLx, 0, 0);
    LLTInvNx(&tmpLx, DDInv + k * 3 + 1);

    // L = -DInv
    GECP(nx[k + 1], nx[k + 1], DDInv + k * 3 + 1, 0, 0, LL + k * 3 + 1, 0, 0);

    GESC(nx[k + 1], nx[k + 1], -1.0, DDInv + k * 3 + 1, 0, 0);
    GESC(nx[k + 1], nx[k + 1], -1.0, DD + k * 3 + 1, 0, 0);

    // D = Q - inv(D)
    GECP(nx[k + 1], nx[k + 1], QQ + k + 1, 0, 0, DD + k * 3 + 2, 0, 0);
    GEAD(nx[k + 1], nx[k + 1], -1.0, DDInv + k * 3 + 1, 0, 0, DD + k * 3 + 2, 0,
         0);
    // DInv = Cholesky(D)
    POTRF_L(nx[k + 1], DD + k * 3 + 2, 0, 0, &tmpLx, 0, 0);
    LLTInvNx(&tmpLx, DDInv + k * 3 + 2);
  }
}

// int main() {
//   initOCPData();
//   printOCPData();

//   factorizeOCPData();

//   printWorkspace();

//   freeOCPData();
//   return 0;
// }
