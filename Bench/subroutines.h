#pragma once

#include <blasfeo.h>

// #define BDLDL_USE_SINGLE_PRECISION

#if defined(BDLDL_USE_SINGLE_PRECISION)
typedef float ldl_float;
typedef struct blasfeo_smat ldl_matrix;
typedef struct blasfeo_svec ldl_vector;

#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define CREATE_STRMAT blasfeo_create_smat
#define PACK_MAT blasfeo_pack_smat
#define PACK_TRAN_MAT blasfeo_pack_tran_smat
#define PRINT_MAT blasfeo_print_smat
#define PRINT_TRANS_MAT blasfeo_print_tran_smat
#define ALLOCATE_MAT blasfeo_allocate_smat
#define FREE_MAT blasfeo_free_smat
#define GECP blasfeo_sgecp
#define POTRF_L blasfeo_spotrf_l
#define TRCP_L blasfeo_strcp_l
#define POTRF_L blasfeo_spotrf_l
#define TRSM_RLTN blasfeo_strsm_rltn
#define TRSM_RLNN blasfeo_strsm_rlnn
#define GEMM_NT blasfeo_sgemm_nt
#define GESE blasfeo_sgese
#define DIAIN blasfeo_sdiain
#define DIARE blasfeo_sdiare
#define GEMM_NN blasfeo_sgemm_nn
#define GESC blasfeo_sgesc
#define GEAD blasfeo_sgead

#define TRSM_LLNN blasfeo_strsm_llnn
#define TRSM_LLTN blasfeo_strsm_lltn

#else
typedef double ldl_float;
typedef struct blasfeo_dmat ldl_matrix;
typedef struct blasfeo_dvec ldl_vector;

#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define CREATE_STRMAT blasfeo_create_dmat
#define PACK_MAT blasfeo_pack_dmat
#define PACK_TRAN_MAT blasfeo_pack_tran_dmat
#define PRINT_MAT blasfeo_print_dmat
#define PRINT_TRANS_MAT blasfeo_print_tran_dmat
#define ALLOCATE_MAT blasfeo_allocate_dmat
#define FREE_MAT blasfeo_free_dmat
#define GECP blasfeo_dgecp
#define TRCP_L blasfeo_dtrcp_l
#define POTRF_L blasfeo_dpotrf_l
#define TRSM_RLTN blasfeo_dtrsm_rltn
#define TRSM_RLNN blasfeo_dtrsm_rlnn
#define GEMM_NT blasfeo_dgemm_nt
#define GESE blasfeo_dgese
#define DIAIN blasfeo_ddiain
#define DIARE blasfeo_ddiare
#define GEMM_NN blasfeo_dgemm_nn
#define GESC blasfeo_dgesc
#define GEAD blasfeo_dgead

#define TRSM_LLNN blasfeo_dtrsm_llnn
#define TRSM_LLTN blasfeo_dtrsm_lltn

#endif