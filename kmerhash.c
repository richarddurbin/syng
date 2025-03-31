/*  File: kmerhash.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: fixed length DNA string hash set package (e.g. syncmers)
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 16 22:50 2025 (rd109)
 * Created: Tue Sep  3 19:38:07 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "kmerhash.h"

static inline void Compress_DNA (int len, char *s, U64 *u) ;    // from s to u
static inline void Compress_DNA_RC (int len, char *s, U64 *u) ; // from s to u
static inline void Uncompress_DNA (int len, U64 *u, char *t) ;  // from u to t

KmerHash *kmerHashCreate (U64 initialSize, int len)
{
  if (!(len & 0x01)) die ("kmerHash len %d must be odd", len) ;
  KmerHash *kh = new0(1, KmerHash) ;
  kh->len = len ;
  kh->dim = 20 ;
  U64 size ;
  for (size = 1 << kh->dim ; size < initialSize ; ++kh->dim, size <<= 1) ;
  kh->table = new0(size, I64) ;
  kh->mask = size - 1 ;
  kh->plen = (len+31) >> 5 ;
  kh->psize = size * 0.3 ;
  kh->pack = new0(kh->plen*kh->psize, U64) ;
  kh->seqbuf = new0(len+1,char) ;
  kh->seqPack = seqPackCreate ('a') ; // 'a' means unpack into acgt
  return kh ; 
}

void kmerHashDestroy (KmerHash *kh)
{
  newFree (kh->table, (U64)1 << kh->dim, I64) ;
  newFree (kh->pack, kh->psize*kh->plen, U64) ;
  newFree (kh->seqbuf, kh->len+1, char) ;
  seqPackDestroy (kh->seqPack) ;
  newFree (kh, 1, KmerHash) ;
}

static U8 comp[] = {   /* sends N (indeed any non-CGT) to A, except 0,1,2,3 are maintained */
   3,   2, 1,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0,   0, 0,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0,   0, 0,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0,   0, 0,   0,   0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 
   0, 'T', 0, 'G',   0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0,
   0,   0, 0,   0, 'A', 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0,
   0, 't', 0, 'g',   0, 0, 0, 'c', 0, 0, 0, 0, 0, 0, 'n', 0,
   0,   0, 0,   0, 'a', 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0
} ;

static inline bool isCanonical (char *dna, int len)
{
  int x = -1, y = len ;
  while (dna[++x] == comp[(int)dna[--y]]) ;
  return (dna[x] < comp[(int)dna[y]]) ;
}

static inline bool isMatch (U64 *u, U64 *v, int n)
{
  while (n--) if (*u++ != *v++) return false ;
  return true ;
}

static inline U64 hashDelta (U64 *u, int plen, int dim)
{ // very simple hash is enough here I think
  U64 v = *u ;
  while (--plen) v = v ^ u[plen] ;
  plen = 64 / dim ;
  while (--plen) v = v ^ (v >> dim) ;
  v |= 1 ; // ensure that delta is odd so coprime with 2^dim
  return v ;
}

#define packseq(kh,i) ((kh)->pack + (i)*(kh)->plen)

static inline bool find (KmerHash *kh, U64 *u, I64 *index, bool isRC, U64 *newLoc, bool isSafe)
{
  if (isSafe) ++kh->finds ;
  U64 loc = *u & kh->mask ; // initially try the simplest thing
  U64 delta = 0 ;
  I64 x = kh->table[loc] ;
  while (x)
    { if (isMatch (u, packseq(kh,x), kh->plen)) // found
	{ if (index) *index = isRC ? -x : x ;
	  return true ;
	}
      if (!delta) delta = hashDelta (u, kh->plen, kh->dim) ;
      loc = (loc + delta) & kh->mask ;
      if (isSafe) ++kh->deltas ;
      x = kh->table[loc] ;
    }
  *newLoc = loc ;
  return false ; // not found
}

bool kmerHashFind (KmerHash *kh, char *dna, I64 *index)
{
  U64 *u = packseq(kh,kh->max+1) ; // pack into the next free slot in pack[]
  bool isRC = !isCanonical (dna,kh->len) ;
  if (isRC) seqPackRevComp (kh->seqPack, dna, (U8*)u, kh->len) ;
  else seqPack (kh->seqPack, dna, (U8*)u, kh->len) ;
  U64 loc ;
  return find (kh, u, index, isRC, &loc, true) ;
}

bool kmerHashFindThreadSafe (KmerHash *kh, char *dna, I64 *index, U64 *u)
{
  bool isRC = !isCanonical (dna,kh->len) ;
  if (isRC) seqPackRevComp (kh->seqPack, dna, (U8*)u, kh->len) ;
  else seqPack (kh->seqPack, dna, (U8*)u, kh->len) ;
  U64 loc ;
  return find (kh, u,  index, isRC, &loc, false) ;
}

bool kmerHashFindPacked (KmerHash *kh, U64 *u, I64 *index) // assume packed and correctly oriented
{
  U64 loc ;
  return find (kh, u, index, false, &loc, true) ;
}

static void doubleTable (KmerHash *kh)
{
  ++kh->dim ;
  I64 *newTable = new0 (1 << kh->dim, I64) ;
  kh->mask = (1 << kh->dim) - 1 ;
  U64  i ;
  for (i = 1 ; i <= kh->max ; ++i) // remap all the packed sequences into newTable
    { U64 loc = *packseq(kh,i) & kh->mask ;
      U64 delta = 0 ;
      while (newTable[loc])
	{ if (!delta) delta = hashDelta(packseq(kh,i), kh->plen, kh->dim) ;
	  loc = (loc + delta) & kh->mask ;
	}
      newTable[loc] = i ;
    }
  newFree (kh->table, (U64)1 << (kh->dim-1), I64) ;
  kh->table = newTable ;
  kh->pack = newResize (kh->pack, kh->plen*kh->psize, 2*kh->plen*kh->psize, U64) ;
  kh->psize *= 2 ;
  // printf ("doubled at max = %llu to %llu\n", kh->max, kh->psize) ;
}

bool kmerHashAdd (KmerHash *kh, char *dna, I64 *index)
{
  U64 *u = packseq(kh,kh->max+1) ; // pack into the next free slot in pack[]
  bool isRC = !isCanonical (dna,kh->len) ;
  if (isRC) seqPackRevComp (kh->seqPack, dna, (U8*)u, kh->len) ;
  else seqPack (kh->seqPack, dna, (U8*)u, kh->len) ;
  U64 newLoc ;
  bool isFound = find (kh, u, index, isRC, &newLoc, true) ;
  if (isFound) return false ;

  kh->table[newLoc] = ++kh->max ; // add the new location
  if (index) *index = isRC ? -kh->max : kh->max ;
  
  if (kh->max == kh->psize-1) doubleTable (kh) ;

  // printf ("add: loc %llx index %lld max %lld\n", loc, *index, kh->max) ;
  return true ;
}

bool kmerHashAddPacked (KmerHash *kh, U64 *u, I64 *index) // assume packed and correctly oriented
{
  U64 newLoc ;
  bool isFound = find (kh, u, index, false, &newLoc, true) ;
  if (isFound) return false ;

  memcpy (packseq(kh,++kh->max), u, kh->plen*sizeof(I64)) ; // have to copy it
  kh->table[newLoc] = kh->max ; // add the new location
  if (index) *index = kh->max ;
  
  if (kh->max == kh->psize-1) doubleTable (kh) ;

  // printf ("add: loc %llx index %lld max %lld\n", loc, *index, kh->max) ;
  return true ;
}

char* kmerHashSeq (KmerHash *kh, I64 i, char *buf) // user must provide buf to be threadsafe
{
  bool isRC = false ;
  if (!buf) buf = kh->seqbuf ;
  if (!i) die ("kmerHashSeq() illegally received i = 0") ;
  if (i < 0) { i = -i ; isRC = true ; }
  if (i > kh->max) die ("out of range in kmerHashSeq: %lld > %lld", i, kh->max) ;
  if (isRC) return seqUnpackRevComp (kh->seqPack, (U8*)packseq(kh,i), buf, 0, kh->len) ;
  else return seqUnpack (kh->seqPack, (U8*)packseq(kh,i), buf, 0, kh->len) ;
}

#ifdef ONEIO

#include "ONElib.h"

static char *schemaText =
  "P 5 khash                 KMER HASH\n"
  "O t 3 3 INT 3 INT 3 INT   max, len, dim for KmerHash table\n"
  "D S 1 3 DNA               packed sequences aligned to 64-bit boundaries\n" 
  "D L 1 8 INT_LIST          locations in the table\n"
  ".\n" ;

bool kmerHashWriteOneFile (KmerHash *kh, OneFile *of)
{
  if (!kh || !of || !of->isWrite || !oneFileCheckSchemaText (of, schemaText)) return false ;
  // write the kmerHash parameters
  oneInt(of,0) = kh->max ; oneInt(of,1) = kh->len ; oneInt(of,2) = kh->dim ;
  oneWriteLine (of, 't', 0, 0) ;

  // next write the DNA - may need to partition into chunks
  I64  chunk = (I64)1 << 22 ;              // chunk size in 8-byte words, giving 128 Mb chunks
  I64  total = kh->max * kh->plen ;        // number of 8-byte words in kh->pack to copy
  U64 *dna = packseq(kh,1) ;
  while (total > chunk)
    { oneWriteLineDNA2bit (of, 'S', chunk << 5, (U8*)dna) ; dna += chunk ; total -= chunk ; }
  if (total) oneWriteLineDNA2bit (of, 'S', total << 5, (U8*)dna) ;
  
  // next find the locations of all the kmers and write them
  I64 i, size = 1 << kh->dim, *loc0 = new(kh->max,I64), *loc = loc0 ;
  for (i = 0 ; i < size ; ++i) if (kh->table[i]) loc[kh->table[i]-1] = i ;
  // similarly chunk the location buffer
  total = kh->max ;
  while (total > chunk) { oneWriteLine (of, 'L', chunk, loc) ; loc += chunk ; total -= chunk ; }
  if (total) oneWriteLine (of, 'L', total, loc) ;
  newFree (loc0, kh->max, I64) ;
  
  return true ;
}

KmerHash *kmerHashReadOneFile (OneFile *of)
{
  if (!of || !oneFileCheckSchemaText (of, schemaText) || !oneGoto (of, 't', 1)) return 0 ;
  
  oneReadLine (of) ;
  KmerHash *kh = kmerHashCreate (1 << oneInt(of,2), oneInt(of,1)) ;
  kh->max = oneInt(of,0) ;

  // read the DNA, maybe in chunks - assume they are at least multiple of 4bp, i.e. full bytes
  U64 *dna = packseq(kh,1) ;
  while (oneReadLine (of) && of->lineType == 'S')
    { memcpy(dna, oneDNA2bit(of), oneLen(of)>>2) ;
      dna += oneLen(of)>>5 ;
    }
  if (dna - packseq(kh,1) != kh->max * kh->plen) die ("wrong number of bp read in kmerhash") ;

  // read the locations
  I64 x = 0 ;
  while (of->lineType == 'L')
    { I64 j, *loc = oneIntList(of) ;
      for (j = 0 ; j < oneLen(of) ; ++j) kh->table[loc[j]] = ++x ;
      oneReadLine(of) ;
    }
  if (x != kh->max) die ("wrong number of locations read in kmerhash") ;

  return kh ;
}

#endif // ONE_DEFINED

/******************************************************************************************/

#ifdef TEST

#include <stdlib.h>

int len ;

static inline void seqSim (char *s)
{ int i ;
  for (i = 0 ; i < len ; ++i) *s++ = random() & 0x3 ;
}

int main (int argc, char *argv[])
{
  U64 i, count ;
  int target ;
  timeUpdate (stdout) ;
  if (argc != 4) die ("Usage: dnahash <len> <number> <target>") ;
  len = atoi(argv[1]) ; count = atoi(argv[2]) ; target = atoi(argv[3]) ;

  KmerHash *kh = kmerHashCreate (0, len) ;

  /* test the compression/decompression code
  char *dna = "cagtatcttagatcctatgggagtacgagcgagcactagcggagggaggagcatcagcacgagcgagggcacaatcta" ;
  char buf[50] ;
  U64 u[2] ;
  printf ("    %s\n", dna) ;
  for (i = 1 ; i < 40 ; ++i)
    { Compress_DNA (i, dna, u) ; Uncompress_DNA (i, u, buf) ; printf ("%2d  %-40s", (int)i, buf) ;
      Compress_DNA_RC (i, dna, u) ; Uncompress_DNA (i, u, buf) ; printf ("  %40s\n", buf) ;
    }
  printf ("    %s\n", dna) ;
  exit (0) ;
  */
  
  char *data = new(count*len, char) ;
  char *s = data ;
  for (i = 0 ; i < count ; ++i, s += len) seqSim(s) ;
  printf ("simulated %llu sequences of length %d\n", count, len) ;
  timeUpdate (stdout) ;
  
  bool isRC ;
  U64 oriTotal[2] ; oriTotal[0] = 0 ; oriTotal[1] = 0 ;
  U64 index ;
  s = data ;
  for (i = 0 ; i < count ; ++i, s += len)
    { kmerHashAdd (kh, s, &index, &isRC) ;
      if (i < 4) printf ("orientation[%llu] is %c\n", i, isRC ? '-' : '+') ;
      ++oriTotal[isRC] ;
    }

  char *acgt = "acgt" ;
  s = data + target*len ;
  printf ("seq %d is  ", target) ; for (i = 0 ; i < len ; ++i) putchar (acgt[s[i]]) ; putchar ('\n') ;
  printf ("hash %d is %s\n", target, kmerHashSeq (kh,target)) ;
  
  printf ("loaded: max %lld finds %lld deltas %lld + %llu - %llu\n",
	  kh->max, kh->finds, kh->deltas, oriTotal[0], oriTotal[1]) ;
  index = 0 ; kmerHashFind (kh, data + target*len, &index, &isRC) ;
  printf ("found target %d at %llu orientation %d\n", target, index, isRC) ;
  timeUpdate (stdout) ;
  kmerHashDestroy (kh) ;
  newFree (data, count*len, char) ;
}

#endif // TEST

/*********** end of file ***********/
