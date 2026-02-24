/*  File: syncmer_iter.c
 *  Description: csyncmer_fast-based syncmer detection using ntHash64.
 *               Replaces the seeded hash in seqhash.c (kept for reference).
 *               Expects ASCII (A/C/G/T) input sequences.
 */

#include "syncmer_iter.h"
#include <stdio.h>

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

/************ syncmer iterator using csyncmer_fast ***********/

SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashIterator *si = new0 (1, SeqhashIterator) ;
  si->sh = sh ;
  si->s = s ;

  // Check minimum length for syncmers (need at least K = w + k bases)
  if (len < sh->w + sh->k) {
    si->isDone = true ;
    si->iter = NULL ;
    return si ;
  }

  // Create csyncmer_fast canonical iterator (expects ASCII input)
  // K = syncmer length = w + k
  // S = smer size = k
  si->iter = csyncmer_iterator_create_canonical_64 (s, (size_t)len, (size_t)(sh->w + sh->k), (size_t)sh->k) ;
  if (!si->iter) {
    si->isDone = true ;
    return si ;
  }

  si->isDone = false ;
  return si ;
}

bool syncmerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ;

  size_t position ;
  int strand ;
  if (csyncmer_iterator_next_canonical_64 (si->iter, &position, &strand)) {
    if (pos) *pos = (int)position ;
    if (kmer) *kmer = 0 ;
    if (isF) *isF = (strand == 0) ;
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
