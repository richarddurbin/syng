/*  File: avx2.h
 *  Description: declarations for the AVX2 SIMD helpers used by seqio.c and kmerhash.c
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
/* non-filtering variants for FASTQ, where the converted length must stay equal to the
   raw sequence line length (the quality line is checked and skipped against it) */
size_t convertIndex4AVX2 (char *seq, char *end, int *convert) ;
size_t convertTextN2AAVX2 (char *seq, char *end, int *convert) ;
bool   isMatchAVX2 (U64 *u, U64 *v, int n) ;

#endif /* AVX2_DEFINED */
