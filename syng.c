/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 11 22:59 2025 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include <pthread.h>
#include <sys/stat.h>

#include "seqio.h"
#ifdef USE_CSYNCMER
#include "syncmer_iter.h"
#else
#include "seqhash.h"
#endif
#include "syng.h"

static OneSchema *schema ;

/************************* thread functions ****************************/

typedef struct {
  Seqhash  *sh ;
  KmerHash *kh ;        // read-only for Find, concurrent insert for Add
  OneFile  *ofIn ;
  SyngBWT  *sbwt ;
  I64       nPath ;
  bool      isAdd ;     // true: use CAS AddThreadSafe; false: use FindThreadSafe
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  	                // seq stores just the starts and ends if only those are required
  Array     seqInfo ;	// of SeqInfo: per sequence
  Array     syncPos ;	// of SyncPos: syncmers and their start positions
} ThreadInfo ;

typedef struct {
  I64       len ;	// sequence length
  I64       nSync ;	// number of syncs
  I32       iSource ;   // ordinal of original source file = sequence set, typically genome
  I32       inSource ;  // position of sequence in original file
} SeqInfo ;

typedef struct {
  I32	    sync ;	// initially 0 if not in kh; -ve if if reverse orientation
  I32       pos ;       // offset in the current sequence
} SyncPos ;

typedef enum { NONE=0, SEQ, PATH, GBWT } OutType ;

static OutType outType = NONE ;

/***** threadProcessSequences() handles input from SeqIO: maps sequences to sync,pos *****/

#define PF_DIST  16
#define PF_DIST2  8
#define PF_TOTAL (PF_DIST + PF_DIST2)

// 3-stage prefetch pipeline: pack kmer → prefetch hash table → match/add.
// Processes a pre-computed position array (from multi-read SIMD path).
static inline void processSyncmerPositions (
    ThreadInfo *ti, char *seq, uint32_t *positions, uint8_t *strands, size_t count,
    U64 *pipePack, int *pipePos, bool *pipeIsRC, I64 *pipeTableVal, I64 *batchState)
{
  KmerHash *kh = ti->kh ;
  int plen = kh->plen ;
  size_t posIdx = 0 ;
  int nFilled = 0 ;
  bool hasMore = true ;

  while (nFilled < PF_TOTAL && posIdx < count)
    { int pos = (int)positions[posIdx] ;
      pipePos[nFilled] = pos ;
      pipeIsRC[nFilled] = !isCanonical (seq+pos, kh->len) ;
      U64 *pk = pipePack + nFilled * plen ;
      if (pipeIsRC[nFilled]) seqPackRevComp (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
      else seqPack (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
      __builtin_prefetch (&kh->table[pk[0] & kh->mask], 0, 0) ;
      nFilled++ ; posIdx++ ;
    }
  if (posIdx >= count) hasMore = false ;

  int nB = nFilled < PF_DIST2 ? nFilled : PF_DIST2 ;
  for (int b = 0 ; b < nB ; ++b)
    pipeTableVal[b] = kmerHashPrefetchPack (kh, pipePack + b * plen) ;

  int nProcessed = 0 ;
  while (nProcessed < nFilled)
    { int slotC = nProcessed % PF_TOTAL ;
      U64 *pkC = pipePack + slotC * plen ;
      // Stage C: match/add — packed kmer is now in cache from earlier prefetch
      I64 sync = 0 ;
      if (kmerHashMatchAt (kh, pkC, pipeTableVal[slotC]))
	sync = pipeIsRC[slotC] ? -pipeTableVal[slotC] : pipeTableVal[slotC] ;
      else if (ti->isAdd)
	kmerHashAddPackedThreadSafe (kh, pkC, &sync, pipeIsRC[slotC], batchState) ;
      else if (pipeTableVal[slotC] > 0)
	kmerHashFindPackedThreadSafe (kh, pkC, &sync, pipeIsRC[slotC]) ;
      if (sync > 2 || sync < -2 || !sync)
	{ SyncPos *sp = arrayp(ti->syncPos, arrayMax(ti->syncPos), SyncPos) ;
	  sp->pos = pipePos[slotC] ;
	  sp->sync = sync ;
	}
      // Stage B: read table slot + prefetch packed kmer for item PF_DIST2 ahead
      if (nProcessed + PF_DIST2 < nFilled)
	{ int slotB = (nProcessed + PF_DIST2) % PF_TOTAL ;
	  pipeTableVal[slotB] = kmerHashPrefetchPack (kh, pipePack + slotB * plen) ;
	}
      // Stage A: pack kmer + prefetch table slot for next item (reuses freed slot)
      if (hasMore)
	{ int pos = (int)positions[posIdx] ;
	  pipePos[slotC] = pos ;
	  pipeIsRC[slotC] = !isCanonical (seq+pos, kh->len) ;
	  U64 *pk = pipePack + slotC * plen ;
	  if (pipeIsRC[slotC]) seqPackRevComp (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
	  else seqPack (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
	  __builtin_prefetch (&kh->table[pk[0] & kh->mask], 0, 0) ;
	  nFilled++ ; posIdx++ ;
	  if (posIdx >= count) hasMore = false ;
	}
      nProcessed++ ;
    }
}

// Same pipeline but fed by syncmerIterator (for single-read / non-CSYNCMER path).
static inline void processIteratorPipeline (
    ThreadInfo *ti, SeqhashIterator *sit, char *seq,
    U64 *pipePack, int *pipePos, bool *pipeIsRC, I64 *pipeTableVal, I64 *batchState)
{
  KmerHash *kh = ti->kh ;
  int plen = kh->plen ;
  int pos ;
  int nFilled = 0 ;
  bool hasMore = true ;

  while (nFilled < PF_TOTAL && syncmerNext (sit, 0, &pos, 0))
    { pipePos[nFilled] = pos ;
      pipeIsRC[nFilled] = !isCanonical (seq+pos, kh->len) ;
      U64 *pk = pipePack + nFilled * plen ;
      if (pipeIsRC[nFilled]) seqPackRevComp (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
      else seqPack (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
      __builtin_prefetch (&kh->table[pk[0] & kh->mask], 0, 0) ;
      nFilled++ ;
    }
  if (nFilled < PF_TOTAL) hasMore = false ;

  int nB = nFilled < PF_DIST2 ? nFilled : PF_DIST2 ;
  for (int b = 0 ; b < nB ; ++b)
    pipeTableVal[b] = kmerHashPrefetchPack (kh, pipePack + b * plen) ;

  int nProcessed = 0 ;
  while (nProcessed < nFilled)
    { int slotC = nProcessed % PF_TOTAL ;
      U64 *pkC = pipePack + slotC * plen ;
      // Stage C: match/add
      I64 sync = 0 ;
      if (kmerHashMatchAt (kh, pkC, pipeTableVal[slotC]))
	sync = pipeIsRC[slotC] ? -pipeTableVal[slotC] : pipeTableVal[slotC] ;
      else if (ti->isAdd)
	kmerHashAddPackedThreadSafe (kh, pkC, &sync, pipeIsRC[slotC], batchState) ;
      else if (pipeTableVal[slotC] > 0)
	kmerHashFindPackedThreadSafe (kh, pkC, &sync, pipeIsRC[slotC]) ;
      if (sync > 2 || sync < -2 || !sync)
	{ SyncPos *sp = arrayp(ti->syncPos, arrayMax(ti->syncPos), SyncPos) ;
	  sp->pos = pipePos[slotC] ;
	  sp->sync = sync ;
	}
      // Stage B: read table + prefetch pack for item PF_DIST2 ahead
      if (nProcessed + PF_DIST2 < nFilled)
	{ int slotB = (nProcessed + PF_DIST2) % PF_TOTAL ;
	  pipeTableVal[slotB] = kmerHashPrefetchPack (kh, pipePack + slotB * plen) ;
	}
      // Stage A: pack + prefetch table for next syncmer
      if (hasMore)
	{ if (syncmerNext (sit, 0, &pos, 0))
	    { pipePos[slotC] = pos ;
	      pipeIsRC[slotC] = !isCanonical (seq+pos, kh->len) ;
	      U64 *pk = pipePack + slotC * plen ;
	      if (pipeIsRC[slotC]) seqPackRevComp (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
	      else seqPack (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
	      __builtin_prefetch (&kh->table[pk[0] & kh->mask], 0, 0) ;
	      nFilled++ ;
	    }
	  else hasMore = false ;
	}
      nProcessed++ ;
    }
}

#ifdef USE_CSYNCMER
#define MULTI_THRESHOLD 32768
#define MULTI_BATCH_SZ 8
#define BUCKET_SHIFT 11          // 2048 bases per bucket
#define N_BUCKETS    16          // covers reads up to K + 16*2048 ≈ 33K bases
#define CHUNK_SIZE   1024        // reads per chunk — arena fits in L2/L3

typedef struct {
  char *seq[MULTI_BATCH_SZ] ;
  int   len[MULTI_BATCH_SZ] ;
  int   idx[MULTI_BATCH_SZ] ;   // index within chunk (0..CHUNK_SIZE-1)
  int   n ;
} LenBucket ;

typedef struct {
  uint32_t *positions ;          // points into arena
  uint8_t  *strands ;
  size_t    count ;
} ReadSyncmers ;

// Fire one bucket through multi-read SIMD, writing directly into arena.
// Each lane gets a reserved slot of maxPerRead entries; no intermediate copy.
static inline void fireBucket (
    LenBucket *bk, Seqhash *sh, ReadSyncmers *readSync,
    size_t maxPerRead, size_t *arenaUsed,
    uint32_t *arenaPos, uint8_t *arenaStr,
    uint8_t *workBuf, size_t workBufSize)
{
  if (bk->n == 0) return ;

  // Point each lane's output directly into arena
  uint32_t *dPos[MULTI_BATCH_SZ] ;
  uint8_t  *dStr[MULTI_BATCH_SZ] ;
  size_t base = *arenaUsed ;
  for (int b = 0 ; b < bk->n ; ++b)
    { dPos[b] = arenaPos + base + b * maxPerRead ;
      dStr[b] = arenaStr + base + b * maxPerRead ;
    }
  for (int b = bk->n ; b < MULTI_BATCH_SZ ; ++b)
    { dPos[b] = NULL ; dStr[b] = NULL ; }

  size_t bCounts[MULTI_BATCH_SZ] ;
  syncmerMultiRead (sh, bk->seq, bk->len, bk->n,
		    dPos, dStr, maxPerRead, bCounts,
		    workBuf, workBufSize) ;

  for (int b = 0 ; b < bk->n ; ++b)
    { ReadSyncmers *rs = &readSync[bk->idx[b]] ;
      rs->count = bCounts[b] ;
      rs->positions = dPos[b] ;
      rs->strands = dStr[b] ;
    }
  *arenaUsed = base + (size_t)bk->n * maxPerRead ;
  bk->n = 0 ;
}
#endif

static void *threadProcessSequences (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  KmerHash *kh = ti->kh ;
  int plen = kh->plen ;
  int K = ti->sh->w + ti->sh->k - 1 ;

  // Pipeline buffers
  U64  *pipePack = new0(PF_TOTAL * plen, U64) ;
  int   pipePos[PF_TOTAL] ;
  bool  pipeIsRC[PF_TOTAL] ;
  I64   pipeTableVal[PF_TOTAL] ;
  I64   batchState[2] = {0, 0} ; // batch ID state for thread-safe insertions

  arrayMax(ti->syncPos) = 0 ;

  int nSeqs = arrayMax(ti->seqInfo) ;
  SeqhashIterator *sit = NULL ;

#ifdef USE_CSYNCMER
  // Short reads (< MULTI_THRESHOLD) use 8-lane SIMD syncmer detection;
  // long reads (genomes, contigs) fall through to the single-read iterator.

  LenBucket buckets[N_BUCKETS] ;
  ReadSyncmers chunkSync[CHUNK_SIZE] ;
  bool chunkIsMulti[CHUNK_SIZE] ;
  char *chunkSeqPtr[CHUNK_SIZE] ;

  size_t maxPerRead = 4 * MULTI_THRESHOLD / (ti->sh->w + 1) ;
  if (maxPerRead < 64) maxPerRead = 64 ;
  size_t arenaSize = CHUNK_SIZE * maxPerRead ;
  uint32_t *arenaPos = (uint32_t *)malloc (arenaSize * sizeof(uint32_t)) ;
  uint8_t  *arenaStr = (uint8_t *)malloc (arenaSize) ;

  size_t workBufSize = syncmerMultiWorkBufSize (ti->sh, MULTI_THRESHOLD) ;
  uint8_t *workBuf = (uint8_t *)aligned_alloc (32, workBufSize) ;

  I64 ss = 0 ;
  for (int chunkStart = 0 ; chunkStart < nSeqs ; chunkStart += CHUNK_SIZE)
    { int chunkEnd = chunkStart + CHUNK_SIZE ;
      if (chunkEnd > nSeqs) chunkEnd = nSeqs ;
      int chunkN = chunkEnd - chunkStart ;

      for (int b = 0 ; b < N_BUCKETS ; ++b) buckets[b].n = 0 ;
      memset (chunkSync, 0, chunkN * sizeof(ReadSyncmers)) ;
      size_t arenaUsed = 0 ;

      for (int ci = 0 ; ci < chunkN ; ++ci)
	{ int gi = chunkStart + ci ;  // global index
	  I64 seqLen = arrp(ti->seqInfo, gi, SeqInfo)->len ;
	  chunkSeqPtr[ci] = arrp(ti->seq, ss, char) ;
	  ss += seqLen ;
	  if (seqLen < MULTI_THRESHOLD && seqLen >= K)
	    { chunkIsMulti[ci] = true ;
	      int bi = (int)((seqLen - K) >> BUCKET_SHIFT) ;
	      if (bi >= N_BUCKETS) bi = N_BUCKETS - 1 ;
	      LenBucket *bk = &buckets[bi] ;
	      bk->seq[bk->n] = chunkSeqPtr[ci] ;
	      bk->len[bk->n] = (int)seqLen ;
	      bk->idx[bk->n] = ci ;
	      bk->n++ ;
	      if (bk->n == MULTI_BATCH_SZ)
		fireBucket (bk, ti->sh, chunkSync, maxPerRead,
			    &arenaUsed, arenaPos, arenaStr,
			    workBuf, workBufSize) ;
	    }
	  else
	    chunkIsMulti[ci] = false ;
	}
      for (int bi = 0 ; bi < N_BUCKETS ; ++bi) // flush partial buckets
	fireBucket (&buckets[bi], ti->sh, chunkSync, maxPerRead,
		     &arenaUsed, arenaPos, arenaStr,
		     workBuf, workBufSize) ;

      // replay chunk in seqInfo order: multi-read results from arena, long reads via iterator
      for (int ci = 0 ; ci < chunkN ; ++ci)
	{ int gi = chunkStart + ci ;
	  I64 seqLen = arrp(ti->seqInfo, gi, SeqInfo)->len ;
	  char *seq = chunkSeqPtr[ci] ;
	  int spStart = arrayMax(ti->syncPos) ;

	  if (chunkIsMulti[ci])
	    { ReadSyncmers *rs = &chunkSync[ci] ;
	      if (rs->count > 0)
		processSyncmerPositions (ti, seq, rs->positions, rs->strands,
					rs->count, pipePack, pipePos, pipeIsRC,
					pipeTableVal, batchState) ;
	    }
	  else if (seqLen >= K)
	    { if (!sit) sit = syncmerIterator (ti->sh, seq, seqLen) ;
	      else syncmerIteratorReinit (sit, seq, seqLen) ;
	      processIteratorPipeline (ti, sit, seq, pipePack, pipePos,
				       pipeIsRC, pipeTableVal, batchState) ;
	    }
	  arrp(ti->seqInfo, gi, SeqInfo)->nSync = arrayMax(ti->syncPos) - spStart ;
	  arrp(ti->seqInfo, gi, SeqInfo)->inSource = 0 ;
	}
    } // chunks

  free (arenaPos) ; free (arenaStr) ; free (workBuf) ;

#else // !USE_CSYNCMER
  I64 seqStart = 0 ;
  for (i = 0 ; i < nSeqs ; ++i)
    { I64 seqLen = arrp(ti->seqInfo, i, SeqInfo)->len ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += seqLen ;
      int spStart = arrayMax(ti->syncPos) ;
      if (!sit) sit = syncmerIterator (ti->sh, seq, seqLen) ;
      else syncmerIteratorReinit (sit, seq, seqLen) ;
      processIteratorPipeline (ti, sit, seq, pipePack, pipePos,
			       pipeIsRC, pipeTableVal, batchState) ;
      arrp(ti->seqInfo, i, SeqInfo)->nSync = arrayMax(ti->syncPos) - spStart ;
      arrp(ti->seqInfo, i, SeqInfo)->inSource = 0 ;
    }
#endif

  if (sit) seqhashIteratorDestroy (sit) ;
  syncmerThreadCleanup () ;

  newFree (pipePack, PF_TOTAL * plen, U64) ;
  return 0 ;
}

/***** threadProcessPaths() handles input from 1path: maps sync/pos to sequence ******/

static void *threadProcessPaths (void* arg) // read in paths, make sequences if necessary
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int klen = ti->kh->len ;
  int i, j, k ;
  I64 seqStart = 0, spStart = 0 ;
  U64 *uBuf = new(ti->kh->plen,U64) ;

  arrayMax(ti->seq) = 0 ;
  arrayMax(ti->seqInfo) = 0 ;
  arrayMax(ti->syncPos) = 0 ;

  oneReadLine (ti->ofIn) ; // the oneGoto() call setting up the thread located just before the object
  for (i = 0 ; i < ti->nPath ; ++i) // need to set nPath
    { if (ti->ofIn->lineType != 'P') break ;
      SeqInfo *si = arrayp(ti->seqInfo, arrayMax(ti->seqInfo), SeqInfo) ;
      si->len = oneInt(ti->ofIn,0) ;
      si->iSource = oneInt(ti->ofIn,1) ;
      si->inSource = oneInt(ti->ofIn,2) ;
      si->nSync = 0 ; // default if no Z or z line
      array(ti->seq, seqStart+si->len, char) = 0 ;
      char *seq = arrp(ti->seq, seqStart, char) ; seqStart += si->len ;
      SyncPos *sp = 0 ;
      // could read the source and sequence position in source from oneInt(0), oneInt(1), but no need
      while (oneReadLine (ti->ofIn) && ti->ofIn->lineType != 'P')
	switch (ti->ofIn->lineType)
	  { // expect a Z line or a (z, o) line pair - X and Y are optional
	  case 'Z':
	    si->nSync = oneInt(ti->ofIn,3) ;
	    arrayp(ti->syncPos, spStart+si->nSync, SyncPos)->pos = 0 ;
	    sp = arrp(ti->syncPos, spStart, SyncPos) ; spStart += si->nSync ;
	    sp->sync = oneInt(ti->ofIn,0) ; sp->pos = oneInt(ti->ofIn,1) ; // starting sync, pos
	    SyngBWTpath *sbp = syngBWTpathStartOld (ti->sbwt, sp->sync, oneInt(ti->ofIn,2)) ;
//	    printf ("starting path %d length %d with sync %d count %d at %d\n", i, (int)si->nSync, (int)sp->sync, (int)oneInt(ti->ofIn,2), (int)sp->pos) ;
	    for (j = 1 ; j < si->nSync ; ++j)
	      if (syngBWTpathNext (sbp, &sp[j].sync, &sp[j].pos))
		{ sp[j].pos += sp[j-1].pos ;
//		  printf (" %d:%d", (int)sp[j].sync, (int)sp[j].pos) ;
		}
	      else
		die ("failed GBWT extension: seq %d count %d max %d from sync %d pos %d",
		     i, j, si->nSync, sp[j-1].sync, sp[j-1].pos) ;
	    syngBWTpathDestroy (sbp) ;
//	    printf ("\n") ;
	    break ;
	  case 'z':
	    si->nSync = oneLen(ti->ofIn) ;
	    arrayp(ti->syncPos, spStart+si->nSync, SyncPos)->pos = 0 ;
	    sp = arrp(ti->syncPos, spStart, SyncPos) ; spStart += si->nSync ;
	    I64 *iList = oneIntList(ti->ofIn) ;
	    for (j = 0 ; j < si->nSync ; ++j) sp[j].sync = iList[j] ;
	    break ;
	  case 'o': // comes with z
	    if (oneLen(ti->ofIn) != si->nSync) die ("o error in threadProcessPaths %d", i) ;
	    iList = oneIntList(ti->ofIn) ;
	    for (j = 0 ; j < si->nSync ; ++j) sp[j].pos = iList[j] ;
	    break ;
	  case 'X':
	    if (si->nSync && oneLen(ti->ofIn) != sp->pos) die ("X error in threadProcessPaths %d", i) ;
	    memcpy (seq, oneDNAchar(ti->ofIn), oneLen(ti->ofIn)) ;
	    break ;
	  case 'Y':
	    if (si->nSync && oneLen(ti->ofIn) != si->len - (sp[si->nSync-1].pos + ti->kh->len))
	      die ("Y error in threadProcessPaths %d nSync %d oneLen %d si->len %d pos %d len %d",
		   i, si->nSync, oneLen(ti->ofIn), si->len, sp[si->nSync-1].pos, ti->kh->len) ;
	    int endLen = oneLen(ti->ofIn) ;
	    memcpy (seq + si->len - endLen, oneDNAchar(ti->ofIn), endLen) ;
	    break ;
	  }
      if (outType == SEQ) // build the sequence from the syncs
	{ I64 lastMax = 0 ;
	  for (j = 0 ; j < si->nSync ; ++j)
	    { assert (sp) ;
	      if (!sp[j].sync) continue ;
	      if (j && sp[j].pos > sp[j-1].pos + klen) // a poly-X run - fill it in
		{ char c = seq[sp[j-1].pos+klen-1] ;
		  for (k = sp[j-1].pos+klen ; k < sp[j].pos ; ++k) seq[k] = c ;
		}
	      kmerHashSeq (ti->kh, sp[j].sync, seq+sp[j].pos) ;
	    }
	}
    }
  
  return 0 ;
}


/************ start of package for sorting-based approach ************/

static int KMERSIZE = 0 ;

static int kmerOrder (const void *a, const void *b)
{
  U64 *ua = (U64*) a, *ub = (U64*) b ;
  int i = KMERSIZE ;
  while (i)
    if (*ua > *ub) return 1 ;
    else if (*ua < *ub) return -1 ;
    else { ua++ ; ub++ ; i-- ; }
  return 0 ;
}

/*****************************************************/

static bool oneFileTest (char* fname)
{
  FILE *f ;
  if (!(f = fopen (fname, "r"))) return false ;
  char peek = getc (f) ;
  fclose (f) ;
  return (peek == '1') ;
}

static char **collectSources (int argc, char *argv[], int *maxSources)
{
  OneFile *of ;
  int  i, nSources = 0 ;
  *maxSources = 1024 ;
  char **sources = new0 (*maxSources, char*) ;
  while (argc--)
    { if (oneFileTest (*argv) && (of = oneFileOpenRead (*argv, 0, 0, 1)))
	{ if (nSources + oneReferenceCount(of) >= *maxSources)
	    newDouble (sources, *maxSources, char*) ;
	  for (i = 0 ; i < oneReferenceCount(of) ; ++i)
	    sources[nSources++] = strdup (of->reference[i].filename) ;
	  oneFileClose (of) ;
	}
      else
	{ if (nSources + 1 >= *maxSources)
	    newDouble (sources, *maxSources, char*) ;
	  sources[nSources++] = strdup (*argv) ;
	}
      ++argv ;
    }
  return sources ;
}

static void addSourceReferences (OneFile *of, char **sources)
{
  int n = oneReferenceCount(of) ;
  while (*sources) { oneAddReference (of, *sources, ++n) ; ++sources ; }
}

/******** estimate total input size for hash pre-sizing ********/

static U64 estimateInputSize (int argc, char **argv)
{
  U64 totalSize = 0 ;
  struct stat st ;
  int i ;
  for (i = 0 ; i < argc ; ++i)
    { if (stat (argv[i], &st) == 0)
        { U64 fileSize = st.st_size ;
          // Check for compressed files and estimate 3x expansion
          int len = strlen (argv[i]) ;
          if (len > 3 && (!strcmp (argv[i]+len-3, ".gz") || !strcmp (argv[i]+len-3, ".bz")))
            fileSize *= 3 ;
          else if (len > 4 && !strcmp (argv[i]+len-4, ".bz2"))
            fileSize *= 3 ;
          totalSize += fileSize ;
        }
    }
  return totalSize ;
}

/******** quadratic histogram package ********/

static void qhist (Array a, FILE *f)
{ I64 i = 0, j, k = 0, n = 0, ni, total = 0 ;
  if (a->size != sizeof(I64)) die ("qhist size mismatch") ;
  for (j = 1 ; k < arrayMax(a) ; ++j)
    { i += j ;
      for (ni = 0 ; k <= i && k < arrayMax(a) ; ++k)
	{ ni += arr(a,k,I64) ; total += k * arr(a,k,I64) ; }
      if (ni) { fprintf (f, "%8lld %lld\n", i, ni) ; n += ni ; }
      if (k < arrayMax(a))
	{ i += j ;
	  for (ni = 0 ; k <= i && k < arrayMax(a) ; ++k)
	    { ni += arr(a,k,I64) ; total += k * arr(a,k,I64) ; }
	  if (ni) { fprintf (f, "%8lld %lld\n", i, ni) ; n += ni ; }
	}
    }
  if (n)
    { I64 median, n50 = 0, ptotal = 0, pn = 0 ;
      for (k = arrayMax(a) ; --k ;)
	{ pn += arr(a,k,I64) ; if (pn > n/2) { median = k+1 ; break ; }
	  if (!n50) { ptotal += k * arr(a,k,I64) ; if (ptotal > total/2) n50 = k+1 ; }
	}
      fprintf (f, "n %lld total %lld mean %.2f median %lld N50 %lld\n",
	       n, total, total/(double)n, median, n50) ;
    }
}

#include "batch.h"  // fillBatch, postProcessBatch, outputBatch

/**************** main program ********************/

static char usage[] =
  "Usage: syng <operation>* <input>*\n"
  "possible operations are:\n"
  "  -w <window length>     : [55] syncmer length = w + k\n"
  "  -k <smer length>       : [8] must be under 32\n"
  "  -T <threads>           : [8] number of threads\n"
  "  -o <outfile prefix>    : [syngOut] applies to all following write* options\n"
  "  -readK <.1khash file>  : read and start from this syncmer (khash) file\n"
  "  -zeroK                 : zero the kmer counts\n"
  "  -noAddK                : do not create new syncmers - convert unmatched syncmers to 0\n"
  //  "  -limitK <min> <max>    : filter current kmer set based on counts ; max 0 to have no upper bound\n"
  "  -histK                 : output quadratic histogram of kmer counts (after sequence processsing)\n"
  "  -writeK                : write the syncmers as a .1khash file\n"
  "  -writeKfa              : write the syncmers as a fasta file, with ending .kmer.fa.gz\n"
  "  -writeNewK <file prefix> : write new syncmers as a .1khash file; implies -noAddK\n"
  "  -writePath             : write a .1path file (paths of nodes)\n"
  "  -writeGBWT             : write a .1gbwt file (nodes, edges and paths in GBWT form)\n"
  "  -writeSeq              : write a .1seq file (paths converted back to sequences)\n"
  "  -outputEnds            : write the non-syncmer ends of path sequences as X,Y lines\n"
  "  -FM                    : convert input GBWT to FM before processing\n"
  //  "  -outputNames           : write the names of path sequences as I lines\n"
  "possible inputs are:\n"
  "  <sequence file>        : any of fasta[.gz], fastq[.gz], BAM/CRAM/SAM, .1seq\n"
  "  <.1path file>          : sequences as lists of kmers, with optional non-syncmer DNA ends\n"
  "  <.1gbwt file>          : graph BWT with paths, with optional ends\n"
  "Operations are carried out in order as they are parsed, with some setting up future actions,\n"
  "e.g. changing the outfile prefix affects following lower case options for file opening\n"
  "Some output files, e.g. .1gbwt will be output at the end, after all inputs are processed,\n"
  "whereas others, e.g. .1path are written as inputs are processed.\n" ;

int main (int argc, char *argv[])
{
  char       *outPrefix = "syngOut" ;
  char       *kFaFileName = 0 ;
  int         nThread = 8 ;
  pthread_t  *threads ;
  SyncmerSet *sms = 0, *smsNew = 0 ;
  OneFile    *ofK = 0, *ofNewK = 0, *ofOut = 0 ;
  SyngBWT    *gbwtOut = 0 ;
  SyncmerParams params = syncmerParamsDefault () ;
  bool        isAddSyncmers = true, isHistK = false, isOutputEnds = false, isFM = false ;
  bool        isSort = false ;
  I64         i, j, k ; // general purpose indices
  
  timeUpdate (0) ;

  schema = oneSchemaCreateFromText (syngSchemaText) ;

  U64 nSeq = 0, totSeq = 0, totSync = 0 ;
  
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { fprintf (stderr, "%s",usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-w") && argc > 1) { params.w = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-k") && argc > 1) { params.k = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-T") && argc > 1) { nThread = atoi(argv[1]) ; argc -=2 ; argv +=2 ; }
    else if (!strcmp (*argv, "-readK") && argc > 1)
      { sms = syncmerSetRead (argv[1]) ; argc -= 2 ; argv +=2 ; }
    else if (!strcmp (*argv, "-zeroK"))
      { if (!sms) die ("-zeroK error: no khash file read in yet") ;
	memset (arrp(sms->count,0,I64), 0, arrayMax(sms->count)*sizeof(I64)) ;
	memset (arrp(sms->maxCount,0,I64), 0, arrayMax(sms->maxCount)*sizeof(I64)) ;
	argc-- ; argv++ ; 
      }
    else if (!strcmp (*argv, "-limitK") && argc > 2)
      { if (!sms) die ("-limitK error: no khash file read in yet") ;
	int min = atoi (argv[1]), max = atoi(argv[2]) ;
	I64 n = 0, *c = arrp(sms->count,1,I64) ;
	for (i = 1 ; i <= sms->kh->max ; ++i, ++c) if (*c >= min && (*c <= max || !max)) ++n ;
	SyncmerSet *sms2 = syncmerSetCreate (params, n*4) ;
	array(sms2->count,n,I64) = 0 ; array(sms2->maxCount,n,char) = 0 ;
	c = arrp(sms->count,1,I64) ;
	I64 j = 0 ;
	for (i = 1 ; i <= sms->kh->max ; ++i, ++c)
	  if (*c >= min && (*c <= max || !max))
	    { kmerHashAddPacked (sms2->kh, sms->kh->pack + i*sms->kh->plen, 0) ;
	      ++j ;
	      arr(sms2->count,j,I64) = arr(sms->count,i,I64) ;
	      arr(sms2->maxCount,j,I64) = arr(sms->count,i,I64) ;
	    }
	syncmerSetDestroy (sms) ;
	sms = sms2 ;
	argc -= 3 ; argv += 3 ;
      }
    else if (!strcmp (*argv, "-histK")) { isHistK = true ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-noAddK"))   { isAddSyncmers = false ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-writeK"))
      { char *fname = fnameTag (outPrefix,"1khash") ;
	ofK = oneFileOpenWriteNew (fname, schema, "khash", true, 1) ;
	if (!ofK) die ("failed to open kmer file %s to write", fname) ;
	oneAddProvenance (ofK, "syng", SYNG_VERSION, getCommandLine()) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writeKfa"))
      { kFaFileName = outPrefix ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writeGBWT"))
      { outType = GBWT ;
	char *fname = fnameTag (outPrefix,"1gbwt") ;
	ofOut = oneFileOpenWriteNew (fname, schema, "gbwt", true, 1) ;
	if (!ofOut) die ("failed to open gbwt file %s to write", fname) ;
	if (sms && sms->kh && sms->kh->max)
	  gbwtOut = syngBWTcreate (params.k + params.w, sms->kh->max + 1) ;
	else
	  gbwtOut = syngBWTcreate (params.k + params.w, 0) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writePath"))
      { outType = PATH ;
	char *fname = fnameTag (outPrefix,"1path") ;
	ofOut = oneFileOpenWriteNew (fname, schema, "path", true, 1) ;
	if (!ofOut) die ("failed to open path file %s to write", fname) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writeSeq"))
      { outType = SEQ ;
	char *fname = fnameTag (outPrefix,"1seq") ;
	ofOut = oneFileOpenWriteNew (fname, schema, "seq", true, 1) ;
	if (!ofOut) die ("failed to open seq file %s to write", fname) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-histK")) { isHistK = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-writeNewK") && argc > 1)
      { smsNew = syncmerSetCreate (params, 0) ;
	ofNewK = oneFileOpenWriteNew (fnameTag(argv[1],"1khash"), schema, "khash", true, 1) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-outputEnds")) { isOutputEnds = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-FM")) { isFM = true ; --argc ; ++argv ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  fprintf (stdout, "k, w are %d %d\n", params.k, params.w) ;
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly

  if (!sms)
    { // Estimate unique syncmer count from largest input file for hash pre-sizing
      // Syncmer density is ~2/(w+1) per base
      // FASTA has ~1.1 bytes/base → bytes_per_syncmer ≈ 0.55*(w+1), use 11*(w+1)/28 for headroom
      // FASTQ has ~4 bytes/base → bytes_per_syncmer ≈ 2*(w+1), use 11*(w+1)/7 for headroom
      // kmerHashResize() between chunks handles growth if estimate is too low.
      U64 maxFileSize = 0 ;
      bool isFastq = false ;
      for (i = 0 ; i < argc ; ++i)
	{ U64 fs = estimateInputSize (1, &argv[i]) ;
	  if (fs > maxFileSize) maxFileSize = fs ;
	  int len = strlen (argv[i]) ;
	  if ((len > 3 && !strcmp (argv[i]+len-3, ".fq")) ||
	      (len > 6 && (!strcmp (argv[i]+len-6, ".fq.gz") || !strcmp (argv[i]+len-6, ".fastq"))) ||
	      (len > 9 && !strcmp (argv[i]+len-9, ".fastq.gz")))
	    isFastq = true ;
	}
      U64 bytesPerSyncmer = (isFastq ? 11 : 11) * ((U64)params.w + 1) / (isFastq ? 7 : 28) ;
      if (bytesPerSyncmer < 1) bytesPerSyncmer = 1 ;
      U64 estSyncmers = maxFileSize / bytesPerSyncmer ;
      U64 initialSize = estSyncmers * 4 ; // 4x headroom for hash load factor
      // Cap pre-sizing so pack array doesn't exceed ~16 GB
      U64 plen = (params.w + params.k + 31) >> 5 ;
      U64 maxPsize = ((U64)16 << 30) / (plen * 8) ;
      if (initialSize > maxPsize * 3)   // psize = size * 0.3, so size = psize/0.3 ~= psize*3
	initialSize = maxPsize * 3 ;
      if (initialSize < (1<<20)) initialSize = 0 ; // use default minimum if small
      sms = syncmerSetCreate (params, initialSize) ;
      if (maxFileSize > 0)
        fprintf (stdout, "estimated from largest file %llu bytes, pre-sized hash for ~%llu syncmers\n",
                 maxFileSize, estSyncmers) ;
#ifdef USE_CSYNCMER
      static const char sentinels[4] = { 'A', 'C', 'G', 'T' } ;
#else
      static const char sentinels[4] = { 0, 1, 2, 3 } ;
#endif
      char *s = new0 (sms->kh->len, char) ;
      for (i = 0 ; i < 4 ; ++i)
	{ for (j = 0 ; j < sms->kh->len ; ++j) s[j] = sentinels[i] ;
	  I64 sync ; kmerHashAdd (sms->kh, s, &sync) ; // poly-A,C,G,T map to 1,2,-2,-1
	}
      newFree (s, sms->kh->len, char) ;
    }

  if (ofOut) oneAddProvenance (ofOut, "syng", SYNG_VERSION, getCommandLine()) ;

  int maxSources ; // collect sources from remaining files
  char **sources = collectSources (argc, argv, &maxSources) ;
  if (sources)
    { if (ofK) addSourceReferences (ofK, sources) ;
      if (ofOut) addSourceReferences (ofOut, sources) ;
      char **s = sources ; while (*s) { free (*s) ; ++s ; }
      newFree (sources, maxSources, char*) ;
    }
  
  U64 nSeq0 = 0, totSeq0 = 0, totSync0 = totSync, syncMax0 = kmerHashMax(sms->kh), nNew = 0 ;

  if (nThread < 1) die ("number of threads %d must be at least 1", nThread) ;
  threads = new (nThread, pthread_t) ;
  ThreadInfo *tiBuf[2] ;
  for (int s = 0 ; s < 2 ; ++s)
    { tiBuf[s] = new0 (nThread, ThreadInfo) ;
      for (i = 0 ; i < nThread ; ++i)
        { tiBuf[s][i].sh = sh ;
          tiBuf[s][i].kh = sms->kh ;
          tiBuf[s][i].isAdd = isAddSyncmers ;
          tiBuf[s][i].seq = arrayCreate (101<<20, char) ; // 101 Mb
          tiBuf[s][i].seqInfo = arrayCreate (20000, SeqInfo) ;
          tiBuf[s][i].syncPos = arrayCreate (1<<20, SyncPos) ;
        }
    }

  int nFile = 0 ;
  int nSource = 0 ;
  while (argc--)
    { ++nFile ;
      OneFile *ofIn = oneFileOpenRead (*argv, schema, 0, nThread) ;
      SeqIO   *sio = 0 ;
      I64      nPath ;
      if (ofIn && oneStats (ofIn, 'P', &nPath, 0, 0) && nPath) // a path file 
	{ if (!sms->kh->max) die ("%s is a path file but we have no kmers", *argv) ;
	  fprintf (stdout, "path file %d %s: ", nFile, *argv) ;
	  oneReadLine (ofIn) ; if (ofIn->lineType == 'h') syncmerParamsCheck (ofIn, params) ;
	  I64 z ;
	  if (oneStats (ofIn, 'Z', &z, 0, 0) && z) // must have a GBWT
	    { tiBuf[0]->sbwt = syngBWTread (ofIn) ;
	      if (isFM) syngBWTtoFM (tiBuf[0]->sbwt) ;
	    }
	  for (i = 0 ; i < nThread ; ++i)
	    { tiBuf[0][i].ofIn = ofIn+i ;
	      oneGoto (ofIn+i, 'P', 1 + (nPath * i) / nThread) ;
	      tiBuf[0][i].nPath = (nPath*(i+1))/nThread - (nPath*i)/nThread ;
	      tiBuf[0][i].sbwt = tiBuf[0]->sbwt ; // copy from the master
	    }
	}
      else
	{ if (ofIn) { oneFileClose (ofIn) ; ofIn = 0 ; }
#ifdef USE_CSYNCMER
	  sio = seqIOopenRead (*argv, dna2textN2AConv, 0) ;
#else
	  sio = seqIOopenRead (*argv, dna2index4Conv, 0) ;
#endif
	  if (!sio) die ("failed to open sequence file %s", *argv) ;
	  ++nSource ; // each input sequence file is a source
	  fprintf (stdout, "sequence file %d %s type %s: ", nSource, *argv, seqIOtypeName[sio->type]) ;
	}
      if (!ofIn && !sio)
	die ("failed to open sequence or path file %s", *argv) ;
      fflush (stdout) ;
      if (sio)
	{ // Double-buffered: fill next batch while threads process current batch
	  int cur = 0 ;
	  U64 batchFilled = fillBatch (sio, tiBuf[cur], nThread, &totSeq) ;
	  while (batchFilled > 0)
	    { for (i = 0 ; i < nThread ; ++i)
		pthread_create (&threads[i], 0, threadProcessSequences, &tiBuf[cur][i]) ;
	      int nxt = 1 - cur ;
	      U64 nxtFilled = fillBatch (sio, tiBuf[nxt], nThread, &totSeq) ;
	      for (i = 0 ; i < nThread ; ++i)
		pthread_join (threads[i], 0) ;
	      postProcessBatch (tiBuf[cur], nThread, sms, smsNew,
				isAddSyncmers, nxtFilled == 0) ;
	      outputBatch (tiBuf[cur], nThread, ofOut, gbwtOut,
			   sms->kh->len, isOutputEnds,
			   &nSeq, &nSeq0, &nSource, &totSync) ;
	      cur = nxt ;
	      batchFilled = nxtFilled ;
	    }
	} // sio
      else if (ofIn)
	{ for (i = 0 ; i < nThread ; ++i)
	    pthread_create (&threads[i], 0, threadProcessPaths, &tiBuf[0][i]) ;
	  for (i = 0 ; i < nThread ; ++i)
	    pthread_join (threads[i], 0) ;
	  outputBatch (tiBuf[0], nThread, ofOut, gbwtOut,
		       sms->kh->len, isOutputEnds,
		       &nSeq, &nSeq0, &nSource, &totSync) ;
	} // ofIn
      
      if (sio)
	{ seqIOclose (sio) ;
	  fprintf (stdout, "had %llu sequences %llu bp, ",
		   nSeq-nSeq0, totSeq-totSeq0) ;
	  fprintf (stdout, "yielding %llu syncs with %llu extra syncmers\n",
		   totSync-totSync0, kmerHashMax(sms->kh)-syncMax0) ;
	}
      else if (ofIn)
	{ oneFileClose (ofIn) ;
	  fprintf (stdout, "had %llu sequences containing %llu syncmers\n",
		   nSeq-nSeq0, totSync-totSync0) ;
	}
      timeUpdate (stdout) ;

      nSeq0 = nSeq ; totSeq0 = totSeq ; totSync0 = totSync ; syncMax0 = kmerHashMax(sms->kh) ;
      ++argv ;
    } // source file

  // compact CAS holes if any were created by concurrent inserts
  // skip for single-file runs unless writing khash (holes corrupt the output)
  if (isAddSyncmers && (nFile > 1 || ofK))
    { I64 nHoles = syncmerSetCompact (sms) ;
      if (nHoles)
	fprintf (stdout, "compacted %lld CAS holes from hash table\n", nHoles) ;
    }

  fprintf (stdout, "Total for this run %llu sequences, total length %llu\n", nSeq, totSeq) ;
  fprintf (stdout, "Overall total %llu instances of %llu syncmers, average %.2f coverage\n",
	   totSync, kmerHashMax(sms->kh), totSync / (double)kmerHashMax(sms->kh)) ;

  // destroy the thread objects
  for (int s = 0 ; s < 2 ; ++s)
    for (i = 0 ; i < nThread ; ++i)
      { arrayDestroy (tiBuf[s][i].seq) ;
        arrayDestroy (tiBuf[s][i].seqInfo) ;
        arrayDestroy (tiBuf[s][i].syncPos) ;
      }
  for (int s = 0 ; s < 2 ; ++s)
    newFree (tiBuf[s], nThread, ThreadInfo) ;
  newFree (threads, nThread, pthread_t) ;

  if (isHistK) // histogram of kmers
    { Array hist = arrayCreate (1024, I64) ;
      for (i = 1 ; i < arrayMax(sms->count) ; ++i) ++array(hist,arr(sms->count,i,I64),I64) ;
      qhist(hist,stdout) ;
      arrayDestroy (hist) ;
#ifdef COUNT_POLY_G
      I64 *cgCount = new0 (sms->kh->len, I64) ;
      int nBytes = (sms->kh->len+3) / 4 ;
      int byteCountC[256], byteCountG[256] ;
      for (i = 0 ; i < 256 ; ++i)
	{ byteCountC[i] = ((i&3) == 1) + ((i&12) == 4) + ((i&48) == 16) + ((i&192) == 64) ;
	  byteCountG[i] = ((i&3) == 2) + ((i&12) == 8) + ((i&48) == 32) + ((i&192) == 128) ;
	}
      for (i = 1 ; i <= sms->kh->max ; ++i)
	{ I64 nC = 0, nG = 0 ;
	  U8 *x ;
	  for (j = 0, x = (U8*)(sms->kh->pack + i*sms->kh->plen) ; j < nBytes ; ++j, ++x)
	    { nC += byteCountC[(int)*x] ; nG += byteCountG[(int)*x] ; }
	  if (nC > nG) ++cgCount[nC] ; else ++cgCount[nG] ;
	}
      for (i = 0 ; i < sms->kh->len ; ++i) printf ("CG %2lld %lld\n", i, cgCount[i]) ;
      newFree (cgCount, sms->kh->len, I64) ;
#endif // COUNT_POLY_G
      timeUpdate (stdout) ;
    }

  // now write all the files requested

  if (outType == GBWT)
    { syngBWTwrite (ofOut, gbwtOut) ;
      fprintf (stdout, "wrote gbwt to file %s\n", oneFileName(ofOut)) ;
      syngBWTdestroy (gbwtOut) ;
      timeUpdate (stdout) ;
    }
  if (ofOut) oneFileClose (ofOut) ;
  
  if (ofK) { syncmerSetWrite (sms, ofK) ; oneFileClose (ofK) ; } // write the kmerhash as ONEcode
  if (sms) syncmerSetDestroy (sms) ;

  if (ofNewK) syncmerSetWrite (smsNew, ofNewK) ; // write the new kmerhash
  if (smsNew) syncmerSetDestroy (smsNew) ;

  if (kFaFileName)
    { char  *fname = fnameTag (kFaFileName,"1kmer.fa.gz") ;
      SeqIO *sio = seqIOopenWrite (fname, FASTA, 0, 0) ;
      if (!sio) die ("failed to open fasta file %s to write", fname) ;
      char idBuf[24], *sBuf = new0 (sms->kh->len+1, char) ;
      for (i = 1 ; i <= sms->kh->max ; ++i)
	{ sprintf (idBuf, "%lld", i) ;
	  seqIOwrite (sio, idBuf, 0, sms->kh->len, kmerHashSeq (sms->kh, i, 0), sBuf) ;
	}
      newFree (sBuf, sms->kh->len+1, char) ;
      seqIOclose (sio) ;
      fprintf (stdout, "wrote %llu syncmers to file %s\n", kmerHashMax(sms->kh), fname) ;
      timeUpdate (stdout) ;
    }

  seqhashDestroy (sh) ;
	  
  fprintf (stdout, "total: ") ; timeTotal (stdout) ;
  return 0 ;
}

/*********************** end of file **********************/
