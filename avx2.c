/*  File: avx2.c
 *  Description: all AVX2 intrinsic code consolidated into one file.
 *               Only this file is compiled with -mavx2 -march=native.
 *               Contains: seqPackAVX2, seqPackRevCompAVX2, convertFilterAVX2
 *               (from seqio.c), isMatchAVX2 (from kmerhash.c), and
 *               syncmer functions (from syncmer_avx2.c).
 */

#include <immintrin.h>
#include "avx2.h"
#include "syncmer_iter.h"
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

/* AVX2 convert + filter for FASTA: uppercase, N->A, strip newlines.
   Works in-place (t <= s always). Returns output length. */
size_t convertFilterAVX2 (char *seq, char *end, int *convert)
{
  char *s = seq, *t = seq ;
  const __m256i vmaskDF = _mm256_set1_epi8 ((char)0xDF) ;
  const __m256i vN      = _mm256_set1_epi8 ('N') ;
  const __m256i vA      = _mm256_set1_epi8 ('A') ;
  const __m256i vlo     = _mm256_set1_epi8 (0x40) ;  /* '@' */
  const __m256i vhi     = _mm256_set1_epi8 (0x5B) ;  /* '[' */

  while (s + 32 <= end)
    { __m256i d = _mm256_loadu_si256 ((const __m256i*)s) ;
      __m256i u = _mm256_and_si256 (d, vmaskDF) ;             /* uppercase */
      __m256i ok = _mm256_and_si256 (_mm256_cmpgt_epi8 (u, vlo),
                                     _mm256_cmpgt_epi8 (vhi, u)) ; /* in [A-Z]? */
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

/************ Seqhash create/read/write ************/

Seqhash *seqhashCreate (int k, int w)
{
  assert (sizeof (U64) == 8) ;
  Seqhash *sh = new0 (1, Seqhash) ;
  sh->k = k ; if (k < 1 || k >= 32) die ("seqhash k %d must be between 1 and 32\n", k) ;
  sh->w = w ; if (w < 1) die ("seqhash w %d must be positive\n", w) ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  for (i = 0 ; i < 4 ; ++i) { sh->patternRC[i] = (3-i) ; sh->patternRC[i] <<= 2*(k-1) ; }
  return sh ;
}

void seqhashWrite (Seqhash *sh, FILE *f)
{ if (fwrite ("SQHSHv3",8,1,f) != 1) die ("failed to write seqhash header") ;
  if (fwrite (sh,sizeof(Seqhash),1,f) != 1) die ("failed to write seqhash") ;
}

Seqhash *seqhashRead (FILE *f)
{ Seqhash *sh = new (1, Seqhash) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read seqhash header") ;
  if (strcmp (name, "SQHSHv3")) die ("seqhash read mismatch - expected SQHSHv3") ;
  if (fread (sh,sizeof(Seqhash),1,f) != 1) die ("failed to read seqhash") ;
  return sh ;
}

void seqhashReport (Seqhash *sh, FILE *f)
{ fprintf (f, "SH k %d  w/m %d\n", sh->k, sh->w) ; }

/************ syncmer iterator using AVX2 TWOSTACK batch ***********/

SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashIterator *si = new0 (1, SeqhashIterator) ;
  si->sh = sh ;
  si->s = s ;
  si->iter = NULL ;  /* not used in AVX2 batch path */

  int K = sh->w + sh->k ;

  if (len < K) {
    si->isDone = true ;
    return si ;
  }

  /* Allocate batch buffers: 4x expected syncmer density gives ample headroom */
  size_t max_syncmers = 4 * (size_t)len / (sh->w + 1) ;
  if (max_syncmers < 64) max_syncmers = 64 ;

  si->batch_positions = (uint32_t *)malloc (max_syncmers * sizeof(uint32_t)) ;
  si->batch_strands   = (uint8_t *)malloc (max_syncmers * sizeof(uint8_t)) ;
  if (!si->batch_positions || !si->batch_strands)
    die ("syncmerIterator: failed to allocate batch buffers for len %d", len) ;

  si->batch_count = csyncmer_twostack_simd_32_canonical_positions (
    s, (size_t)len, (size_t)K, (size_t)sh->k,
    si->batch_positions, si->batch_strands, max_syncmers) ;

  si->batch_index = 0 ;
  si->isDone = (si->batch_count == 0) ;
  return si ;
}

bool syncmerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ;

  if (si->batch_index < si->batch_count) {
    if (pos)  *pos  = (int)si->batch_positions[si->batch_index] ;
    if (kmer) *kmer = 0 ;
    if (isF)  *isF  = (si->batch_strands[si->batch_index] == 0) ;
    si->batch_index++ ;
    if (si->batch_index >= si->batch_count)
      si->isDone = true ;
    return true ;
  }

  si->isDone = true ;
  return false ;
}

/*********************************************/

char *seqString (U64 kmer, int len)
{
  static char trans[4] = { 'a', 'c', 'g', 't' } ;
  static char buf[33] ;
  assert (len <= 32) ;
  buf[len] = 0 ;
  while (len--) { buf[len] = trans[kmer & 0x3] ; kmer >>= 2 ; }
  return buf ;
}

void seqhashForCompilerHappiness (Seqhash *sh, SeqhashIterator *sit)
{ seqhashDestroy (sh) ; seqhashIteratorDestroy (sit) ; }

/**************** end of file ****************/
