/*  File: syncmer_iter.h
 *  Description: csyncmer_fast-based syncmer detection using ntHash64.
 *               Replaces the seeded hash in seqhash.h (kept for reference).
 */

#ifndef SYNCMER_ITER_DEFINED
#define SYNCMER_ITER_DEFINED

#include "utils.h"
#include "csyncmer_fast.h"

typedef struct {
  int k ;			/* kmer/smer size */
  int w ;			/* window */
  U64 mask ;			/* 2*k bits */
  U64 patternRC[4] ;		/* one per base */
} Seqhash ;

typedef struct {
  Seqhash *sh ;
  char *s ;     		/* sequence currently being hashed */
  CsyncmerIteratorCanonical64 *iter ;    /* scalar iterator (NULL for AVX2 batch) */
  uint32_t *batch_positions ;            /* AVX2 batch: syncmer positions */
  uint8_t  *batch_strands ;             /* AVX2 batch: strand flags */
  size_t    batch_count ;               /* AVX2 batch: number of syncmers found */
  size_t    batch_index ;               /* AVX2 batch: current iteration index */
  size_t    batch_capacity ;            /* allocated capacity of batch buffers */
  bool isDone ;
} SeqhashIterator ;

Seqhash *seqhashCreate (int k, int w) ;
static void seqhashDestroy (Seqhash *sh) { free (sh) ; }

void seqhashWrite (Seqhash *sh, FILE *f) ;
Seqhash *seqhashRead (FILE *f) ;
void seqhashReport (Seqhash *sh, FILE *f) ;

// for all iterators sequence must continue to exist through the life of the iterator
// all the *next functions return any/all of kmer, pos, isF - get hash from seqhash(sh,kmer)

// (closed) syncmer extracts w-mers that end with a minimal kmer
// these provide a cover, and have good distribution properties
// Now uses ntHash64 via csyncmer_fast.h
SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len) ;
void syncmerIteratorReinit (SeqhashIterator *si, char *s, int len) ;
bool syncmerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF) ;

static void seqhashIteratorDestroy (SeqhashIterator *si)
{ if (si->iter) csyncmer_iterator_destroy_canonical_64 (si->iter) ;
  free (si->batch_positions) ;
  free (si->batch_strands) ;
  free (si) ; }

// utilities
char *seqString (U64 kmer, int len)  ;
static inline char* seqhashString (Seqhash *sh, U64 kmer) { return seqString (kmer, sh->k) ; }

#endif /* SYNCMER_ITER_DEFINED */

/******* end of file ********/
