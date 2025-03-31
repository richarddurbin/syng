/*  File: kmerhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: fixed length DNA string hash set package (e.g. syncmers)
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 16 22:07 2025 (rd109)
 * Created: Tue Sep  3 19:39:02 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "array.h"  // includes utils.h
#include "seqio.h"  // for SeqPack

// fixed length kmer hasing, based on dict.[ch]
// does not support removal
// pack sequences to 2bit encoding as in ONElib, so can read/write ONEfiles natively
// kmers are stored in canonical orientation: kmer < reverseComplement(kmer) 
// new kmers are added from 1 .. kmerHashMax() inclusive: so 1-based numbering
// this allows to return matches in negative orientation as -index
// kmerHashAdd() and kmerHashFind() are not threadsafe, but can use kmerHasFindThreadSafe() in threads

typedef struct {
  int      len ;     // length of dna sequences stored
  int      dim ;     // dimension of table (size is 2^dim)
  I64     *table ;   // the main table indexed by the hash
  I64      mask ;    // mask limited to tsize bits
  I64      max ;     // current number of entries
  int      plen ;    // packed sequence length = (len+31) >> 5
  U64      psize ;   // max number of elements in pack before doubling
  U64     *pack ;    // the packed sequences
  char    *seqbuf ;  // size len+1, for printing out sequences
  U64      finds ;   // stats: number of finds
  U64      deltas ;  // stats: number of deltas (not in remapping)
  SeqPack *seqPack ; // for binary packing/unpacking
} KmerHash ;

KmerHash *kmerHashCreate (U64 initialSize, int len) ;
void      kmerHashDestroy (KmerHash *kh) ;

bool      kmerHashAdd (KmerHash *kh, char *dna, I64 *index) ;  // true if addded - always fill *index
bool      kmerHashFind (KmerHash *kh, char *dna, I64 *index) ; // true if found
bool      kmerHashFindThreadSafe (KmerHash *kh, char *dna, I64 *index, U64 *buf) ;
// buf must point to user memory of size kh->plen or larger
// NB: Add, Find and FindPacked increment count. FindThreadSafe does not - user must do this.
bool      kmerHashAddPacked (KmerHash *kh, U64 *u, I64 *index) ;
bool      kmerHashFindPacked (KmerHash *kh, U64 *u, I64 *index) ; // true if found
// these versions add/find already packed and correctly oriented kmers

char*     kmerHashSeq (KmerHash *kh, I64 i, char *buf) ; // retrieve i'th sequence (rev-comp if i < 0)
                                                         // buf can be 0, but then not thread-safe
#define   kmerHashMax(kh)  ((kh)->max)      // number of stored kmers

#ifdef ONE_DEFINED
bool      kmerHashWriteOneFile (KmerHash *kh, OneFile *of) ;
KmerHash *kmerHashReadOneFile (OneFile *of) ;
#endif

/*********** end of file ***********/
