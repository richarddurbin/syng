/*  File: avx2.c
 *  Description: all AVX2 intrinsic code consolidated into one file.
 *               Only this file is compiled with -mavx2 -march=native.
 *               Contains: seqPackAVX2, seqPackRevCompAVX2, convertFilterAVX2
 *               (from seqio.c), isMatchAVX2 (from kmerhash.c), and
 *               syncmer functions (from syncmer_avx2.c).
 */

#include <immintrin.h>
#include "avx2.h"
#include <stdio.h>

/************ AVX2 helpers from seqio.c ************/

// AVX2 SIMD packing: process 32 ASCII bases at once -> 8 packed bytes
void seqPackAVX2 (const char *s, U8 *u, U64 len)
{
  // Lookup table indexed by (char & 0x0F): A=1->0, C=3->1, T=4->3, G=7->2
  const __m256i lut = _mm256_setr_epi8 (
    0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,  // low 128-bit lane
    0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0   // high 128-bit lane
  ) ;
  const __m256i mask0f = _mm256_set1_epi8 (0x0F) ;

  while (len >= 32)
    { __m256i chars = _mm256_loadu_si256 ((const __m256i*)s) ;
      __m256i nibbles = _mm256_and_si256 (chars, mask0f) ;
      __m256i bases = _mm256_shuffle_epi8 (lut, nibbles) ; // 32 x 2-bit values in low 2 bits of each byte

      // Now combine 4 bases per byte:
      // bases = [b0, b1, b2, b3, b4, b5, b6, b7, ...]
      // want: [b0 | b1<<2 | b2<<4 | b3<<6, b4 | b5<<2 | b6<<4 | b7<<6, ...]

      // Shift odd positions left by 2: multiply by [1, 4, 1, 4, ...]
      const __m256i mult1 = _mm256_set1_epi16 (0x0401) ; // low byte * 1, high byte * 4
      __m256i step1 = _mm256_maddubs_epi16 (bases, mult1) ; // pairs combined: 16 x 16-bit

      // Now step1 has pairs: [b0+b1*4, b2+b3*4, ...] in 16-bit values (only low 6 bits used)
      // Combine pairs: multiply by [1, 16] and add
      const __m256i mult2 = _mm256_set1_epi32 (0x00100001) ; // low 16-bit * 1, high 16-bit * 16
      __m256i step2 = _mm256_madd_epi16 (step1, mult2) ; // 8 x 32-bit values

      // Pack down to 8 bytes using shuffle
      // step2 has bytes: [packed, 0, 0, 0, packed, 0, 0, 0, ...] in each lane
      __m256i shuffled = _mm256_shuffle_epi8 (step2, _mm256_setr_epi8 (
        0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      )) ;

      // Extract low 4 bytes from each 128-bit lane
      U32 lo = _mm256_extract_epi32 (shuffled, 0) ;
      U32 hi = _mm256_extract_epi32 (shuffled, 4) ;
      *(U32*)u = lo ;
      *(U32*)(u+4) = hi ;

      len -= 32 ; s += 32 ; u += 8 ;
    }

  // Scalar fallback for remainder
  static const U8 p[] = { 0,0,0,1,3,0,0,2 } ;
  while (len >= 4)
    { *u++ = p[s[0]&7] | (p[s[1]&7] << 2) | (p[s[2]&7] << 4) | (p[s[3]&7] << 6) ;
      len -= 4 ; s += 4 ;
    }
  if (len >= 1) { *u = p[s[0]&7] ; }
  if (len >= 2) { *u |= p[s[1]&7] << 2 ; }
  if (len >= 3) { *u |= p[s[2]&7] << 4 ; }
}

// AVX2 reverse complement packing for ASCII input
void seqPackRevCompAVX2 (const char *s, U8 *u, U64 len)
{
  // Complement lookup: A=1->3, C=3->2, T=4->0, G=7->1
  const __m256i lut = _mm256_setr_epi8 (
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0
  ) ;
  const __m256i mask0f = _mm256_set1_epi8 (0x0F) ;
  const __m256i rev_idx = _mm256_setr_epi8 (
    15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
    15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
  ) ;

  const char *sOrig = s ;

  if (len >= 32)
    { s += len - 32 ; // start from end
      while (len >= 32)
        { __m256i chars = _mm256_loadu_si256 ((const __m256i*)s) ;
          // Reverse bytes within each 128-bit lane
          chars = _mm256_shuffle_epi8 (chars, rev_idx) ;
          // Swap the two 128-bit lanes
          chars = _mm256_permute2x128_si256 (chars, chars, 0x01) ;

          __m256i nibbles = _mm256_and_si256 (chars, mask0f) ;
          __m256i bases = _mm256_shuffle_epi8 (lut, nibbles) ;

          const __m256i mult1 = _mm256_set1_epi16 (0x0401) ;
          __m256i step1 = _mm256_maddubs_epi16 (bases, mult1) ;

          const __m256i mult2 = _mm256_set1_epi32 (0x00100001) ; // low 16-bit * 1, high 16-bit * 16
          __m256i step2 = _mm256_madd_epi16 (step1, mult2) ;

          __m256i shuffled = _mm256_shuffle_epi8 (step2, _mm256_setr_epi8 (
            0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
          )) ;

          U32 lo = _mm256_extract_epi32 (shuffled, 0) ;
          U32 hi = _mm256_extract_epi32 (shuffled, 4) ;
          *(U32*)u = lo ;
          *(U32*)(u+4) = hi ;

          len -= 32 ; s -= 32 ; u += 8 ;
        }
    }

  // Scalar fallback: process remaining len chars from sOrig[0..len-1] in reverse
  static const U8 p[] = { 0,3,0,2,0,0,0,1 } ;
  s = sOrig + len ;
  while (len >= 4)
    { s -= 4 ;
      *u++ = p[s[3]&7] | (p[s[2]&7] << 2) | (p[s[1]&7] << 4) | (p[s[0]&7] << 6) ;
      len -= 4 ;
    }
  if (len >= 1) { *u = p[sOrig[len-1]&7] ; }
  if (len >= 2) { *u |= p[sOrig[len-2]&7] << 2 ; }
  if (len >= 3) { *u |= p[sOrig[len-3]&7] << 4 ; }
}

// AVX2 packing for 0-3 encoded input (no LUT needed — bytes are already 2-bit values)
void seqPackIndex4AVX2 (const char *s, U8 *u, U64 len)
{
  while (len >= 32)
    { __m256i bases = _mm256_loadu_si256 ((const __m256i*)s) ;

      const __m256i mult1 = _mm256_set1_epi16 (0x0401) ;
      __m256i step1 = _mm256_maddubs_epi16 (bases, mult1) ;

      const __m256i mult2 = _mm256_set1_epi32 (0x00100001) ;
      __m256i step2 = _mm256_madd_epi16 (step1, mult2) ;

      __m256i shuffled = _mm256_shuffle_epi8 (step2, _mm256_setr_epi8 (
        0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      )) ;

      U32 lo = _mm256_extract_epi32 (shuffled, 0) ;
      U32 hi = _mm256_extract_epi32 (shuffled, 4) ;
      *(U32*)u = lo ;
      *(U32*)(u+4) = hi ;

      len -= 32 ; s += 32 ; u += 8 ;
    }

  while (len >= 4)
    { *u++ = s[0] | (s[1] << 2) | (s[2] << 4) | (s[3] << 6) ;
      len -= 4 ; s += 4 ;
    }
  if (len >= 1) { *u = s[0] ; }
  if (len >= 2) { *u |= s[1] << 2 ; }
  if (len >= 3) { *u |= s[2] << 4 ; }
}

// AVX2 reverse complement packing for 0-3 encoded input
void seqPackRevCompIndex4AVX2 (const char *s, U8 *u, U64 len)
{
  const __m256i vxor3 = _mm256_set1_epi8 (3) ;
  const __m256i rev_idx = _mm256_setr_epi8 (
    15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
    15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
  ) ;

  const char *sOrig = s ;

  if (len >= 32)
    { s += len - 32 ;
      while (len >= 32)
        { __m256i chars = _mm256_loadu_si256 ((const __m256i*)s) ;
          chars = _mm256_shuffle_epi8 (chars, rev_idx) ;
          chars = _mm256_permute2x128_si256 (chars, chars, 0x01) ;

          __m256i bases = _mm256_xor_si256 (chars, vxor3) ; // complement: 0<->3, 1<->2

          const __m256i mult1 = _mm256_set1_epi16 (0x0401) ;
          __m256i step1 = _mm256_maddubs_epi16 (bases, mult1) ;

          const __m256i mult2 = _mm256_set1_epi32 (0x00100001) ;
          __m256i step2 = _mm256_madd_epi16 (step1, mult2) ;

          __m256i shuffled = _mm256_shuffle_epi8 (step2, _mm256_setr_epi8 (
            0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
          )) ;

          U32 lo = _mm256_extract_epi32 (shuffled, 0) ;
          U32 hi = _mm256_extract_epi32 (shuffled, 4) ;
          *(U32*)u = lo ;
          *(U32*)(u+4) = hi ;

          len -= 32 ; s -= 32 ; u += 8 ;
        }
    }

  s = sOrig + len ;
  while (len >= 4)
    { s -= 4 ;
      *u++ = (3-s[3]) | ((3-s[2]) << 2) | ((3-s[1]) << 4) | ((3-s[0]) << 6) ;
      len -= 4 ;
    }
  if (len >= 1) { *u = (3-sOrig[len-1]) ; }
  if (len >= 2) { *u |= (3-sOrig[len-2]) << 2 ; }
  if (len >= 3) { *u |= (3-sOrig[len-3]) << 4 ; }
}

/* AVX2 convert + filter for FASTA: uppercase, N->A, strip newlines.
   Works in-place (t <= s always). Returns output length. */
size_t convertFilterAVX2 (char *seq, char *end, int *convert)
{
  char *s = seq, *t = seq ;
  const __m256i vmaskDF = _mm256_set1_epi8 ((char)0xDF) ;
  const __m256i vN      = _mm256_set1_epi8 ('N') ;
  const __m256i vA      = _mm256_set1_epi8 ('A') ;
  const __m256i vcA = _mm256_set1_epi8 ('A') , vcC = _mm256_set1_epi8 ('C') ,
                vcG = _mm256_set1_epi8 ('G') , vcT = _mm256_set1_epi8 ('T') ,
                vcN = _mm256_set1_epi8 ('N') ;

  while (s + 32 <= end)
    { __m256i d = _mm256_loadu_si256 ((const __m256i*)s) ;
      __m256i u = _mm256_and_si256 (d, vmaskDF) ;             /* uppercase */
      /* exactly ACGTN?  anything else (IUPAC codes, newlines) must go through
         convert[] so it is dropped, not silently folded onto a base by the LUT */
      __m256i ok = _mm256_or_si256
        (_mm256_or_si256 (_mm256_cmpeq_epi8 (u, vcA), _mm256_cmpeq_epi8 (u, vcC)),
         _mm256_or_si256 (_mm256_or_si256 (_mm256_cmpeq_epi8 (u, vcG), _mm256_cmpeq_epi8 (u, vcT)),
                          _mm256_cmpeq_epi8 (u, vcN))) ;
      if (_mm256_movemask_epi8 (ok) == (int)0xFFFFFFFF)
        { __m256i r = _mm256_blendv_epi8 (u, vA, _mm256_cmpeq_epi8 (u, vN)) ;
          _mm256_storeu_si256 ((__m256i*)t, r) ;
          s += 32 ; t += 32 ;
        }
      else                                                     /* newlines or non-alpha present */
        { char *e = s + 32 ;
          while (s < e) { int c = convert[(int)(unsigned char)*s++] ; if (c >= 0) *t++ = (char)c ; }
        }
    }
  while (s < end) { int c = convert[(int)(unsigned char)*s++] ; if (c >= 0) *t++ = (char)c ; }
  return (size_t)(t - seq) ;
}

/* AVX2 convert + filter for FASTA with dna2index4Conv: ACGT->0123, N->0, strip newlines.
   Works in-place (t <= s always). Returns output length. */
size_t convertFilterIndex4AVX2 (char *seq, char *end, int *convert)
{
  char *s = seq, *t = seq ;
  const __m256i vmaskDF = _mm256_set1_epi8 ((char)0xDF) ;
  const __m256i mask0f  = _mm256_set1_epi8 (0x0F) ;
  const __m256i vcA = _mm256_set1_epi8 ('A') , vcC = _mm256_set1_epi8 ('C') ,
                vcG = _mm256_set1_epi8 ('G') , vcT = _mm256_set1_epi8 ('T') ,
                vcN = _mm256_set1_epi8 ('N') ;
  // LUT indexed by (uppercased_char & 0x0F): A=1->0, C=3->1, T=4->3, G=7->2, N=14->0
  const __m256i lut = _mm256_setr_epi8 (
    0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0
  ) ;

  while (s + 32 <= end)
    { __m256i d = _mm256_loadu_si256 ((const __m256i*)s) ;
      __m256i u = _mm256_and_si256 (d, vmaskDF) ;             /* uppercase */
      /* exactly ACGTN?  anything else (IUPAC codes, newlines) must go through
         convert[] so it is dropped, not silently folded onto a base by the LUT */
      __m256i ok = _mm256_or_si256
        (_mm256_or_si256 (_mm256_cmpeq_epi8 (u, vcA), _mm256_cmpeq_epi8 (u, vcC)),
         _mm256_or_si256 (_mm256_or_si256 (_mm256_cmpeq_epi8 (u, vcG), _mm256_cmpeq_epi8 (u, vcT)),
                          _mm256_cmpeq_epi8 (u, vcN))) ;
      if (_mm256_movemask_epi8 (ok) == (int)0xFFFFFFFF)
        { __m256i nibbles = _mm256_and_si256 (u, mask0f) ;
          __m256i r = _mm256_shuffle_epi8 (lut, nibbles) ;
          _mm256_storeu_si256 ((__m256i*)t, r) ;
          s += 32 ; t += 32 ;
        }
      else                                                     /* newlines or non-alpha present */
        { char *e = s + 32 ;
          while (s < e) { int c = convert[(int)(unsigned char)*s++] ; if (c >= 0) *t++ = (char)c ; }
        }
    }
  while (s < end) { int c = convert[(int)(unsigned char)*s++] ; if (c >= 0) *t++ = (char)c ; }
  return (size_t)(t - seq) ;
}


/* Non-filtering in-place convert, for FASTQ.  Upstream converts a FASTQ sequence line
   with a plain `*s = convert[*s]` loop, so the converted length always equals the raw
   line length - the quality line is validated and skipped against it.  Dropping
   characters here (as the FASTA convertFilter* routines do) desynchronises the parser.
   Returns the length, which is unchanged, for symmetry with the filtering versions. */

static inline size_t convertNoFilterAVX2 (char *seq, char *end, int *convert, __m256i lut)
{
  char *s = seq ;
  const __m256i vmaskDF = _mm256_set1_epi8 ((char)0xDF) ;
  const __m256i mask0f  = _mm256_set1_epi8 (0x0F) ;
  const __m256i vcA = _mm256_set1_epi8 ('A') , vcC = _mm256_set1_epi8 ('C') ,
                vcG = _mm256_set1_epi8 ('G') , vcT = _mm256_set1_epi8 ('T') ,
                vcN = _mm256_set1_epi8 ('N') ;

  while (s + 32 <= end)
    { __m256i d = _mm256_loadu_si256 ((const __m256i*)s) ;
      __m256i u = _mm256_and_si256 (d, vmaskDF) ;             /* uppercase */
      __m256i ok = _mm256_or_si256
        (_mm256_or_si256 (_mm256_cmpeq_epi8 (u, vcA), _mm256_cmpeq_epi8 (u, vcC)),
         _mm256_or_si256 (_mm256_or_si256 (_mm256_cmpeq_epi8 (u, vcG), _mm256_cmpeq_epi8 (u, vcT)),
                          _mm256_cmpeq_epi8 (u, vcN))) ;
      if (_mm256_movemask_epi8 (ok) == (int)0xFFFFFFFF)
        _mm256_storeu_si256 ((__m256i*)s, _mm256_shuffle_epi8 (lut, _mm256_and_si256 (u, mask0f))) ;
      else /* something outside ACGTN - take the table verbatim, keeping the length */
        { int i ; for (i = 0 ; i < 32 ; ++i) s[i] = (char) convert[(int)(unsigned char)s[i]] ; }
      s += 32 ;
    }
  while (s < end) { *s = (char) convert[(int)(unsigned char)*s] ; ++s ; }
  return (size_t)(end - seq) ;
}

size_t convertIndex4AVX2 (char *seq, char *end, int *convert)
{ // A,N->0 C->1 G->2 T->3, indexed by (uppercased char & 0x0F)
  const __m256i lut = _mm256_setr_epi8 (
    0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0) ;
  return convertNoFilterAVX2 (seq, end, convert, lut) ;
}

size_t convertTextN2AAVX2 (char *seq, char *end, int *convert)
{ // A,N->'A' C->'C' G->'G' T->'T', indexed by (uppercased char & 0x0F)
  const __m256i lut = _mm256_setr_epi8 (
    0, 'A', 0, 'C', 'T', 0, 0, 'G', 0, 0, 0, 0, 0, 0, 'A', 0,
    0, 'A', 0, 'C', 'T', 0, 0, 'G', 0, 0, 0, 0, 0, 0, 'A', 0) ;
  return convertNoFilterAVX2 (seq, end, convert, lut) ;
}

/************ AVX2 helper from kmerhash.c ************/

bool isMatchAVX2 (U64 *u, U64 *v, int n)
{
  // For n == 3, scalar is fine
  if (n == 3) return u[0] == v[0] && u[1] == v[1] && u[2] == v[2] ;
  // For n >= 4, use AVX2 in chunks of 4
  while (n >= 4)
    { __m256i a = _mm256_loadu_si256 ((const __m256i*)u) ;
      __m256i b = _mm256_loadu_si256 ((const __m256i*)v) ;
      __m256i cmp = _mm256_cmpeq_epi64 (a, b) ;
      if (_mm256_movemask_epi8 (cmp) != (int)0xFFFFFFFF) return false ;
      u += 4 ; v += 4 ; n -= 4 ;
    }
  // Handle remainder
  while (n--) if (*u++ != *v++) return false ;
  return true ;
}

/**************** end of file ****************/
