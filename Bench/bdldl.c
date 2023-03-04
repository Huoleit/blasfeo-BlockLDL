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
void* mem;

void initOCPData() {
  int ii;

  for (ii = 0; ii < N + 1; ++ii) {
    nx[ii] = NX;
    nu[ii] = NU;
    ng[ii] = NG;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Calculate memory size for the matrices
  //////////////////////////////////////////////////////////////////////////////
  size_t size = 0;
  for (ii = 0; ii < N; ++ii) {
    size += SIZE_STRMAT(nx[ii + 1], nu[ii]) * 2;      // B L_B
    size += SIZE_STRMAT(nx[ii + 1], nx[ii]) * 2;      // A L_A
    size += SIZE_STRMAT(nu[ii], nu[ii]);              // R
    size += SIZE_STRMAT(nx[ii + 1], nx[ii + 1]) * 2;  // Q L_I

    size += SIZE_STRMAT(nu[ii], nu[ii]) * 2;          // D_R DInv_R
    size += SIZE_STRMAT(nx[ii + 1], nx[ii + 1]) * 2;  // D_B DInv_B
    size += SIZE_STRMAT(nx[ii + 1], nx[ii + 1]) * 2;  // D_Q DInv_Q
  }
  blasfeo_malloc_align(&mem, size);

  /////////////////////////////////////////////////////////////////////////////
  // Allocate auxiliary matrices internally
  //////////////////////////////////////////////////////////////////////////////
  ALLOCATE_MAT(NX, NX, &Ix);
  GESE(NX, NX, 0, &Ix, 0, 0);
  DIARE(NX, 1.0, &Ix, 0, 0);

  ALLOCATE_MAT(NU, NU, &Iu);
  GESE(NU, NU, 0, &Iu, 0, 0);
  DIARE(NU, 1.0, &Iu, 0, 0);

  ALLOCATE_MAT(NX, NX, &tmpLx);
  ALLOCATE_MAT(NU, NU, &tmpLu);

  /////////////////////////////////////////////////////////////////////////////
  // Map working data to contiguous memory
  //////////////////////////////////////////////////////////////////////////////
  char* ptr = (char*)mem;
  for (ii = 0; ii < N; ++ii) {
    CREATE_STRMAT(nx[ii + 1], nu[ii], BB + ii, ptr);  // B
    PACK_MAT(nx[ii + 1], nu[ii], B, nx[ii + 1], BB + ii, 0, 0);
    ptr += (BB + ii)->memsize;

    CREATE_STRMAT(nx[ii + 1], nx[ii], AA + ii, ptr);
    PACK_MAT(nx[ii + 1], nx[ii], A, nx[ii + 1], AA + ii, 0, 0);  // A
    ptr += (AA + ii)->memsize;

    CREATE_STRMAT(nu[ii], nu[ii], RR + ii, ptr);  // R
    PACK_MAT(nu[ii], nu[ii], R, nu[ii], RR + ii, 0, 0);
    ptr += (RR + ii)->memsize;

    CREATE_STRMAT(nx[ii], nx[ii], QQ + ii + 1, ptr);
    PACK_MAT(nx[ii], nx[ii], Q, nx[ii], QQ + ii + 1, 0, 0);  // Q
    ptr += (QQ + ii + 1)->memsize;
  }

  for (ii = 0; ii < N; ++ii) {
    if (ii != 0) {
      CREATE_STRMAT(nx[ii + 1], nx[ii], LL + ii * 3 - 1, ptr);  // L_A
      ptr += (LL + ii * 3 - 1)->memsize;
    }
    CREATE_STRMAT(nx[ii + 1], nu[ii], LL + ii * 3, ptr);  // L_B
    ptr += (LL + ii * 3)->memsize;

    CREATE_STRMAT(nx[ii + 1], nx[ii + 1], LL + ii * 3 + 1, ptr);  // L_I
    ptr += (LL + ii * 3 + 1)->memsize;

    CREATE_STRMAT(nu[ii], nu[ii], DD + ii * 3, ptr);  // D_R
    ptr += (DD + ii * 3)->memsize;
    CREATE_STRMAT(nx[ii + 1], nx[ii + 1], DD + ii * 3 + 1, ptr);  // D_B
    ptr += (DD + ii * 3 + 1)->memsize;
    CREATE_STRMAT(nx[ii + 1], nx[ii + 1], DD + ii * 3 + 2, ptr);  // D_Q
    ptr += (DD + ii * 3 + 2)->memsize;

    CREATE_STRMAT(nu[ii], nu[ii], DDInv + ii * 3, ptr);  // DInv_R
    ptr += (DDInv + ii * 3)->memsize;
    CREATE_STRMAT(nx[ii + 1], nx[ii + 1], DDInv + ii * 3 + 1, ptr);  // DInv_B
    ptr += (DDInv + ii * 3 + 1)->memsize;
    CREATE_STRMAT(nx[ii + 1], nx[ii + 1], DDInv + ii * 3 + 2, ptr);  // DInv_Q
    ptr += (DDInv + ii * 3 + 2)->memsize;
  }
}

void freeOCPData() {
  FREE_MAT(&Ix);
  FREE_MAT(&Iu);

  FREE_MAT(&tmpLx);
  FREE_MAT(&tmpLu);

  blasfeo_free_align(mem);
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
    if (ii != 0) {
      printf("Q = \n");
      PRINT_MAT(nx[ii], nx[ii], QQ + ii, 0, 0);
    }
    printf("R = \n");
    PRINT_MAT(nu[ii], nu[ii], RR + ii, 0, 0);
  }
  printf("\nk = %d(Final) nx = %d\n", ii, nx[ii]);
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
//   // printOCPData();
//   blasfeo_timer timer;
//   int nrep = 10000;
//   blasfeo_tic(&timer);
//   for (int i = 0; i < nrep; ++i) {
//     factorizeOCPData();
//   }
//   double time_elapsed = blasfeo_toc(&timer) / (double)nrep;
//   printf("Time %f [ms]\n", time_elapsed * 1e3);

//   // printWorkspace();

//   freeOCPData();
//   return 0;
// }
