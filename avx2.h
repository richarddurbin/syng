/*  File: avx2.h
 *  Description: declarations for AVX2 SIMD helper functions
 *               consolidated from seqio.c, kmerhash.c, syncmer_avx2.c
 */

#ifndef AVX2_DEFINED
#define AVX2_DEFINED

#include "utils.h"

void   seqPackAVX2 (const char *s, U8 *u, U64 len) ;
void   seqPackRevCompAVX2 (const char *s, U8 *u, U64 len) ;
size_t convertFilterAVX2 (char *seq, char *end, int *convert) ;
bool   isMatchAVX2 (U64 *u, U64 *v, int n) ;

#endif /* AVX2_DEFINED */
