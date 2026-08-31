/*  File: syncmer_iter.c
 *  Description: syncmer detection using ntHash, from the csyncmer_fast library.
 *               Drop-in replacement for the seeded hash in seqhash.c, selected by
 *               CSYNCMER=1 (which requires AVX2=1).  Expects ASCII (acgtACGT) input.
 */

#include "syncmer_iter.h"
#include <stdio.h>


/************ Seqhash create/read/write ************/

Seqhash *seqhashCreate (int k, int w, int seed)
{ (void)seed ; // unused for ntHash64 path, kept for API compatibility
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

  int K = sh->w + sh->k - 1 ; // K-mer length must match kh->len = params.w + params.k

  if (len < K) {
    si->isDone = true ;
    return si ;
  }

  /* 4*len/(w+1) assumes random-sequence density; low-complexity approaches 1/pos */
  size_t max_syncmers = (size_t)len - (size_t)K + 1 ; // true bound: 1 per k-mer start
  if (max_syncmers < 64) max_syncmers = 64 ;

  si->batch_positions = (uint32_t *)malloc (max_syncmers * sizeof(uint32_t)) ;
  si->batch_strands   = (uint8_t *)malloc (max_syncmers * sizeof(uint8_t)) ;
  if (!si->batch_positions || !si->batch_strands)
    die ("syncmerIterator: failed to allocate batch buffers for len %d", len) ;
  si->batch_capacity = max_syncmers ;

  si->batch_count = csyncmer_twostack_simd_32_canonical_positions (
    s, (size_t)len, (size_t)K, (size_t)sh->k,
    si->batch_positions, si->batch_strands, max_syncmers) ;

  si->batch_index = 0 ;
  si->isDone = (si->batch_count == 0) ;
  return si ;
}

void syncmerIteratorReinit (SeqhashIterator *si, char *s, int len)
{
  Seqhash *sh = si->sh ;
  si->s = s ;
  int K = sh->w + sh->k - 1 ;

  if (len < K) {
    si->isDone = true ;
    si->batch_count = 0 ;
    si->batch_index = 0 ;
    return ;
  }

  size_t max_syncmers = (size_t)len - (size_t)K + 1 ; // true bound: 1 per k-mer start
  if (max_syncmers < 64) max_syncmers = 64 ;

  if (max_syncmers > si->batch_capacity) {
    free (si->batch_positions) ;
    free (si->batch_strands) ;
    si->batch_positions = (uint32_t *)malloc (max_syncmers * sizeof(uint32_t)) ;
    si->batch_strands   = (uint8_t *)malloc (max_syncmers * sizeof(uint8_t)) ;
    if (!si->batch_positions || !si->batch_strands)
      die ("syncmerIteratorReinit: failed to allocate batch buffers for len %d", len) ;
    si->batch_capacity = max_syncmers ;
  }

  si->batch_count = csyncmer_twostack_simd_32_canonical_positions (
    s, (size_t)len, (size_t)K, (size_t)sh->k,
    si->batch_positions, si->batch_strands, si->batch_capacity) ;

  si->batch_index = 0 ;
  si->isDone = (si->batch_count == 0) ;
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

void syncmerMultiRead (Seqhash *sh, char *seqs[8], int lens[8], int n,
                       uint32_t *positions[8], uint8_t *strands[8],
                       size_t max_per_read, size_t counts[8],
                       uint8_t *work_buf, size_t work_buf_size)
{
  int K = sh->w + sh->k - 1 ;
  const char *mseqs[8] = {NULL} ;
  size_t mlens[8] = {0} ;
  for (int i = 0 ; i < n && i < 8 ; ++i)
    { mseqs[i] = seqs[i] ; mlens[i] = (size_t)lens[i] ; }
  csyncmer_twostack_simd_32_multi_canonical_positions (
    mseqs, mlens, (size_t)K, (size_t)sh->k,
    positions, strands, max_per_read, counts,
    work_buf, work_buf_size) ;
}

size_t syncmerMultiWorkBufSize (Seqhash *sh, size_t max_read_len)
{
  int K = sh->w + sh->k - 1 ;
  return csyncmer_multi_work_buf_size (max_read_len, (size_t)K, (size_t)sh->k) ;
}

void syncmerThreadCleanup (void)
{ csyncmer_twostack_thread_cleanup () ;
}

void seqhashForCompilerHappiness (Seqhash *sh, SeqhashIterator *sit)
{ seqhashDestroy (sh) ; seqhashIteratorDestroy (sit) ; }

/**************** end of file ****************/
