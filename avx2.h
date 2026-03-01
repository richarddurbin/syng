/*  File: avx2.h
 *  Description: declarations for AVX2 SIMD helper functions
 *               consolidated from seqio.c, kmerhash.c, syncmer_avx2.c
 */

#ifndef AVX2_DEFINED
#define AVX2_DEFINED

#include "utils.h"

void   seqPackAVX2 (const char *s, U8 *u, U64 len) ;
void   seqPackRevCompAVX2 (const char *s, U8 *u, U64 len) ;
void   seqPackIndex4AVX2 (const char *s, U8 *u, U64 len) ;
void   seqPackRevCompIndex4AVX2 (const char *s, U8 *u, U64 len) ;
size_t convertFilterAVX2 (char *seq, char *end, int *convert) ;
size_t convertFilterIndex4AVX2 (char *seq, char *end, int *convert) ;
bool   isMatchAVX2 (U64 *u, U64 *v, int n) ;

/* AVX2 branchless minimum of n U64 values.  Uses unsigned comparison
   (_mm256_cmpgt_epi64 with sign-bit flip) and _mm256_blendv_epi8.
   Eliminates branch mispredictions in the syncmer window minimum scan. */
U64 seqhash_min_u64 (const U64 *arr, int n) ;

/* AVX2 batch canonical hash: given n forward and RC k-mers, compute
   canonical hashes (min of fwd/RC) and direction flags.
   Uses _mm256_mul_epu32 with hash decomposition for k <= 16. */
void seqhash_hash_avx2 (const U64 *fwd, const U64 *rc, int n,
                         U64 factor1, int shift1, int k,
                         U64 *out_hash, bool *out_fwd) ;

#endif /* AVX2_DEFINED */
