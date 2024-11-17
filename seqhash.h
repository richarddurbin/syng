/*  File: seqhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: header file for seqhash package - minimizers and moshers
 * Exported functions: see below
 * HISTORY:
 * Last edited: Jul 24 14:38 2023 (rd109)
 * Created: Mon Mar  5 08:43:45 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"

typedef struct {
  int seed ;			/* seed */
  int k ;			/* kmer */
  int w ;			/* window */
  U64 mask ;			/* 2*k bits */
  int shift1, shift2 ;
  U64 factor1, factor2 ;
  U64 patternRC[4] ;		/* one per base */
} Seqhash ;

typedef struct {
  Seqhash *sh ;
  char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
  U64 h, hRC ;			/* current k-mer values */
  U64 *hash ;			/* buffer of length w holding hashes for current window */
  bool *isForward ;		/* buffer of length w holding isForward for current window */
  int base ;			/* start of buf in sequence */
  int iStart, iMin ;		/* position in buf of start of current window, next min */
  U64 min ;                     /* needed for syncmers */
  bool isDone ;
} SeqhashIterator ;

Seqhash *seqhashCreate (int k, int w, int seed) ;
static void seqhashDestroy (Seqhash *sh) { free (sh) ; }

void seqhashWrite (Seqhash *sh, FILE *f) ;
Seqhash *seqhashRead (FILE *f) ;
void seqhashReport (Seqhash *sh, FILE *f) ;

// for all iterators sequence must continue to exist through the life of the iterator
// all the *next functions return any/all of kmer, pos, isF - get hash from seqhash(sh,kmer)

// simple iterator to return all kmers
SeqhashIterator *seqhashIterator (Seqhash *sh, char *s, int len) ;
bool seqhashNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF) ;

static void seqhashIteratorDestroy (SeqhashIterator *si)
{ free (si->hash) ; free (si->isForward) ; free (si) ; }

// iterator to extract minimizers from a sequence
// NB sequence must continue to exist through the life of the iterator
SeqhashIterator *minimizerIterator (Seqhash *sh, char *s, int len) ;
bool minimizerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF) ;

// modimizer extracts hashes that are divisible by m->w
// this is faster and more robust to errors - same mean density without evenness guarantees
SeqhashIterator *modIterator (Seqhash *sh, char *s, int len) ;
bool modNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF) ;

// (closed) syncmer extracts w-mers that end with a minimal kmer
// these provide a cover, and have good distribution properties
SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len) ;
bool syncmerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF) ;

// utilities
static inline U64 kHash (Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1) ; }
char *seqString (U64 kmer, int len)  ;
static inline char* seqhashString (Seqhash *sh, U64 kmer) { return seqString (kmer, sh->k) ; }

/******* end of file ********/
