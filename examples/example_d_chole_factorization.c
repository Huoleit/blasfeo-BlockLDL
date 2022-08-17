#include <blasfeo.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  printf("\nExample of LU factorization and backsolve\n\n");

#if defined(LA_HIGH_PERFORMANCE)

  printf("\nLA provided by BLASFEO\n\n");

#elif defined(LA_REFERENCE)

  printf("\nLA provided by REFERENCE\n\n");

#elif defined(LA_EXTERNAL_BLAS_WRAPPER)

  printf("\nLA provided by EXTERNAL_BLAS_WRAPPER\n\n");

#else

  printf("\nLA provided by ???\n\n");
  exit(2);

#endif

  printf("Testing processor\n");

  char supportString[50];
  blasfeo_processor_library_string(supportString);
  printf("Library requires processor features:%s\n", supportString);

  int features = 0;
  int procCheckSucceed = blasfeo_processor_cpu_features(&features);
  blasfeo_processor_feature_string(features, supportString);
  printf("Processor supports features:%s\n", supportString);

  if (!procCheckSucceed) {
    printf(
        "Current processor does not support the current compiled BLASFEO "
        "library.\n");
    printf("Please get a BLASFEO library compatible with this processor.\n");
    exit(3);
  }

  // //	blasfeo_dgetrf_nopivot(n, n, &sD, 0, 0, &sD, 0, 0);

  // 	blasfeo_dgetrf_rp(n, n, &sD, 0, 0, &sLU, 0, 0, ipiv);
  // 	printf("\nLU = \n");
  // 	blasfeo_print_dmat(n, n, &sLU, 0, 0);
  // 	printf("\nipiv = \n");
  // 	int_print_mat(1, n, ipiv, 1);

  // #if 0 // solve A X = P L U X = B  =>  L U X = P^T B
  // 	blasfeo_pack_dmat(n, n, I, n, &sI, 0, 0);
  // 	printf("\nI = \n");
  // 	blasfeo_print_dmat(n, n, &sI, 0, 0);

  // 	blasfeo_drowpe(n, ipiv, &sI);
  // 	printf("\nperm(I) = \n");
  // 	blasfeo_print_dmat(n, n, &sI, 0, 0);

  // 	blasfeo_dtrsm_llnu(n, n, 1.0, &sLU, 0, 0, &sI, 0, 0, &sD, 0, 0);
  // 	printf("\nperm(inv(L)) = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);
  // 	blasfeo_dtrsm_lunn(n, n, 1.0, &sLU, 0, 0, &sD, 0, 0, &sD, 0, 0);
  // 	printf("\ninv(A) = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);

  // 	// convert from strmat to column major matrix
  // 	blasfeo_unpack_dmat(n, n, &sD, 0, 0, D, n);
  // #elif 0 // solve X^T A^T = X^T (P L U)^T = B^T  =>  X^T U^T L^T = B^T P
  // 	blasfeo_pack_tran_dmat(n, n, I, n, &sI, 0, 0);
  // 	printf("\nI' = \n");
  // 	blasfeo_print_dmat(n, n, &sI, 0, 0);

  // 	blasfeo_dcolpe(n, ipiv, &sB);
  // 	printf("\nperm(I') = \n");
  // 	blasfeo_print_dmat(n, n, &sB, 0, 0);

  // 	blasfeo_dtrsm_rltu(n, n, 1.0, &sLU, 0, 0, &sB, 0, 0, &sD, 0, 0);
  // 	printf("\nperm(inv(L')) = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);
  // 	blasfeo_dtrsm_rutn(n, n, 1.0, &sLU, 0, 0, &sD, 0, 0, &sD, 0, 0);
  // 	printf("\ninv(A') = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);

  // 	// convert from strmat to column major matrix
  // 	blasfeo_unpack_tran_dmat(n, n, &sD, 0, 0, D, n);
  // #else // solve A^T X = (P L U)^T X = U^T L^T P^T X = B
  // 	blasfeo_dgetr(n, n, &sLU, 0, 0, &sLUt, 0, 0);

  // 	blasfeo_pack_dmat(n, n, I, n, &sI, 0, 0);
  // 	printf("\nI = \n");
  // 	blasfeo_print_dmat(n, n, &sI, 0, 0);

  // 	blasfeo_dtrsm_llnn(n, n, 1.0, &sLUt, 0, 0, &sI, 0, 0, &sD, 0, 0);
  // 	printf("\ninv(U^T) = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);

  // 	blasfeo_dtrsm_lunu(n, n, 1.0, &sLUt, 0, 0, &sD, 0, 0, &sD, 0, 0);
  // 	printf("\n(inv(L^T)*inv(U^T)) = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);

  // 	blasfeo_drowpei(n, ipiv, &sD);
  // 	printf("\nperm(inv(L^T)*inv(U^T)) = \n");
  // 	blasfeo_print_dmat(n, n, &sD, 0, 0);

  // 	// convert from strmat to column major matrix
  // 	blasfeo_unpack_dmat(n, n, &sD, 0, 0, D, n);
  // #endif

  // 	// print matrix in column-major format
  // 	printf("\ninv(A) = \n");
  // 	d_print_mat(n, n, D, n);

  //
  // free memory
  //

  d_free(A);
  d_free(B);
  d_free(D);
  d_free(I);
  int_free(ipiv);
  //	blasfeo_free_dmat(&sA);
  //	blasfeo_free_dmat(&sB);
  //	blasfeo_free_dmat(&sD);
  //	blasfeo_free_dmat(&sI);
  v_free_align(memory_strmat);

  return 0;
}
