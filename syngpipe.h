/*  File: syngpipe.h
 *  Description: everything the AVX2/csyncmer optimisation adds to syng.c, kept here so
 *    that syng.c itself stays close to the original.  Included from syng.c immediately
 *    after the ThreadInfo definition, which these functions work on.
 *    Nothing here is used by any other file.
 *-------------------------------------------------------------------
 */

#include <sys/stat.h>

/***** helpers over KmerHash internals, used by the prefetch pipeline below *****/

#define pipePack(kh,i) ((kh)->pack + (i)*(kh)->plen) // as packseq() in kmerhash.c

// Two-level prefetch: stage B reads the table and prefetches pack; stage C checks the match
static inline I64 kmerHashPrefetchPack (KmerHash *kh, U64 *packed)
{ U64 loc = packed[0] & kh->mask ;
  I64 x = kh->table[loc] ;
  if (x > 0) __builtin_prefetch (pipePack(kh, x), 0, 0) ;
  return x ;
}

// prefetch the first table slot a packed kmer hashes to, ahead of a later Add or Find
static inline void kmerHashPrefetchTable (KmerHash *kh, U64 *packed)
{ __builtin_prefetch (&kh->table[packed[0] & kh->mask], 1, 1) ; }

static inline bool kmerHashMatchAt (KmerHash *kh, U64 *packed, I64 tableVal)
{ if (tableVal <= 0) return false ;
  U64 *v = pipePack(kh, tableVal) ;
  int n = kh->plen ;
  switch (n)
    { case 1: return packed[0] == v[0] ;
      case 2: return packed[0] == v[0] && packed[1] == v[1] ;
      default: while (n--) if (*packed++ != *v++) return false ; return true ;
    }
}

// as syncmerAdd(), but for a kmer already packed in canonical orientation by a worker thread
static inline void syncmerAddPacked (SyncmerSet *sms, U64 *u, bool isRC, I64 *index)
{ I64 i ;
  bool added = kmerHashAddPacked (sms->kh, u, &i) ;
  if (added)
    { array(sms->count, i, I64)++ ;
      array(sms->thisCount, i, char) = 1 ;
    }
  else
    { arr(sms->count, i, I64)++ ;
      if (++arr(sms->thisCount, i, I64) & 0x80) arr(sms->thisCount,i,I64) = 0x7f ;
    }
  if (index) *index = isRC ? -i : i ;
}


// The hash table is far bigger than L3, so every lookup is a DRAM miss.  To hide that
// latency we run a 3-stage software pipeline over the syncmer positions of a sequence:
//   A: pack the kmer, and prefetch the table slot it hashes to
//   B: PF_DIST2 items later, read that slot and prefetch the packed kmer it points at
//   C: PF_DIST items after A, do the match, by when both lines have arrived
// Threads never insert: ids must be handed out in input order for the run to be
// reproducible, so a miss just stages its packed kmer for postProcessBatch().
#define PF_DIST  16
#define PF_DIST2  8
#define PF_TOTAL (PF_DIST + PF_DIST2)

typedef struct {          // the in-flight window, indexed by slot % PF_TOTAL
  U64  *pack ;		// PF_TOTAL packed kmers, each kh->plen U64s
  int   pos[PF_TOTAL] ;
  bool  isRC[PF_TOTAL] ;
  I64   tableVal[PF_TOTAL] ;
} Pipe ;

static inline void pipeStageA (ThreadInfo *ti, Pipe *pp, char *seq, int slot, int pos)
{ KmerHash *kh = ti->kh ;
  U64 *pk = pp->pack + slot * kh->plen ;
  pp->pos[slot] = pos ;
  pp->isRC[slot] = !isCanonical (seq+pos, kh->len) ;
  if (pp->isRC[slot]) seqPackRevComp (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
  else seqPack (kh->seqPack, seq+pos, (U8*)pk, kh->len) ;
  __builtin_prefetch (&kh->table[pk[0] & kh->mask], 0, 0) ;
}

static inline void pipeStageB (ThreadInfo *ti, Pipe *pp, int nProcessed, int nFilled)
{ if (nProcessed + PF_DIST2 >= nFilled) return ;
  int slot = (nProcessed + PF_DIST2) % PF_TOTAL ;
  pp->tableVal[slot] = kmerHashPrefetchPack (ti->kh, pp->pack + slot * ti->kh->plen) ;
}

static inline void pipeStageC (ThreadInfo *ti, Pipe *pp, int slot)
{ KmerHash *kh = ti->kh ;
  U64 *pk = pp->pack + slot * kh->plen ;
  I64 tv = pp->tableVal[slot], sync = 0 ;
  if (kmerHashMatchAt (kh, pk, tv)) sync = pp->isRC[slot] ? -tv : tv ;
  else if (tv > 0)
    kmerHashFindPackedThreadSafe (kh, pk, &sync, pp->isRC[slot]) ;
  if (sync > 2 || sync < -2 || !sync) // +-1, +-2 are the homopolymers, which are not recorded
    { SyncPos *sp = arrayp(ti->syncPos, arrayMax(ti->syncPos), SyncPos) ;
      sp->pos = pp->pos[slot] ;
      sp->sync = sync ;
      if (!sync && ti->isAdd) // stage it, so the serial add need not find and pack it again
	{ I64 n = arrayMax(ti->newIsRC), plen = kh->plen ;
	  array(ti->newIsRC, n, char) = pp->isRC[slot] ;
	  array(ti->newPack, (n+1)*plen - 1, U64) = 0 ; // ensure the whole block exists
	  memcpy (arrp(ti->newPack, n*plen, U64), pk, plen * sizeof(U64)) ;
	}
    }
}

static inline void pipePrimeB (ThreadInfo *ti, Pipe *pp, int nFilled) // stage B for the first items
{ int b, n = nFilled < PF_DIST2 ? nFilled : PF_DIST2 ;
  for (b = 0 ; b < n ; ++b)
    pp->tableVal[b] = kmerHashPrefetchPack (ti->kh, pp->pack + b * ti->kh->plen) ;
}

// run the pipeline over a pre-computed position array, from the multi-read SIMD path
static inline void processSyncmerPositions (ThreadInfo *ti, Pipe *pp, char *seq,
					    uint32_t *positions, size_t count)
{
  size_t posIdx = 0 ;
  int nFilled = 0, nProcessed ;

  while (nFilled < PF_TOTAL && posIdx < count)
    pipeStageA (ti, pp, seq, nFilled++, (int)positions[posIdx++]) ;
  pipePrimeB (ti, pp, nFilled) ;
  for (nProcessed = 0 ; nProcessed < nFilled ; ++nProcessed)
    { int slot = nProcessed % PF_TOTAL ;
      pipeStageC (ti, pp, slot) ;
      pipeStageB (ti, pp, nProcessed, nFilled) ;
      if (posIdx < count) // refill the slot we have just freed
	{ pipeStageA (ti, pp, seq, slot, (int)positions[posIdx++]) ; ++nFilled ; }
    }
}

// the same pipeline, fed one position at a time by syncmerIterator
static inline void processIteratorPipeline (ThreadInfo *ti, Pipe *pp, char *seq,
					    SeqhashIterator *sit)
{
  int pos, nFilled = 0, nProcessed ;
  bool hasMore = true ;

  while (nFilled < PF_TOTAL && syncmerNext (sit, 0, &pos, 0))
    pipeStageA (ti, pp, seq, nFilled++, pos) ;
  if (nFilled < PF_TOTAL) hasMore = false ;
  pipePrimeB (ti, pp, nFilled) ;
  for (nProcessed = 0 ; nProcessed < nFilled ; ++nProcessed)
    { int slot = nProcessed % PF_TOTAL ;
      pipeStageC (ti, pp, slot) ;
      pipeStageB (ti, pp, nProcessed, nFilled) ;
      if (hasMore) // refill the slot we have just freed
	{ if (syncmerNext (sit, 0, &pos, 0)) { pipeStageA (ti, pp, seq, slot, pos) ; ++nFilled ; }
	  else hasMore = false ;
	}
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

  Pipe pipe ;
  pipe.pack = new0 (PF_TOTAL * plen, U64) ;

  arrayMax(ti->syncPos) = 0 ;
  arrayMax(ti->newPack) = 0 ;
  arrayMax(ti->newIsRC) = 0 ;
  
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

	  // count > maxPerRead means only the first maxPerRead positions were
	  // written, so redo the read through the correctly-bounded iterator
	  if (chunkIsMulti[ci] && chunkSync[ci].count <= maxPerRead)
	    { ReadSyncmers *rs = &chunkSync[ci] ;
	      if (rs->count > 0)
		processSyncmerPositions (ti, &pipe, seq, rs->positions, rs->count) ;
	    }
	  else if (seqLen >= K)
	    { if (!sit) sit = syncmerIterator (ti->sh, seq, seqLen) ;
	      else syncmerIteratorReinit (sit, seq, seqLen) ;
	      processIteratorPipeline (ti, &pipe, seq, sit) ;
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
      processIteratorPipeline (ti, &pipe, seq, sit) ;
      arrp(ti->seqInfo, i, SeqInfo)->nSync = arrayMax(ti->syncPos) - spStart ;
      arrp(ti->seqInfo, i, SeqInfo)->inSource = 0 ; // ensure not used in output loop
    }
#endif

  if (sit) seqhashIteratorDestroy (sit) ;
  syncmerThreadCleanup () ;

  newFree (pipe.pack, PF_TOTAL * plen, U64) ;
  return 0 ;
}

/************ batch helpers, so the main loop can double-buffer ************/

// if the file is mmap'ed, point seqIOread() straight at ti->seq so it need not copy
static void setDirectBuf (SeqIO *sio, Array seq, int seqStart)
{ U64 room = sio->maxSeqLen ? (U64)sio->maxSeqLen * 2 + 64 : 65536 ;
  array(seq, seqStart + room, char) = 0 ; // grow now; arrayMax is reset at the end
  sio->directBuf = arrp(seq, seqStart, char) ;
  sio->directBufSize = room ;
}

// read 100Mb DNA per thread, returning total bytes read (0 means end of file)
static U64 fillBatch (SeqIO *sio, ThreadInfo *threadInfo, int nThread, U64 *totSeq)
{
  U64 filled = 0 ;
  int i ;
  for (i = 0 ; i < nThread ; ++i)
    { ThreadInfo *ti = &threadInfo[i] ;
      arrayMax(ti->seq) = 0 ;
      arrayMax(ti->seqInfo) = 0 ;
      int seqStart = 0 ;
      if (sio->mmapSize) setDirectBuf (sio, ti->seq, seqStart) ;
      while (arrayMax(ti->seq) < 100<<20 && seqIOread (sio))
	{ arrayp(ti->seqInfo, arrayMax(ti->seqInfo), SeqInfo)->len = sio->seqLen ;
	  array(ti->seq, seqStart+sio->seqLen, char) = 0 ;
	  if (!sio->directBufUsed)
	    memcpy (arrp(ti->seq, seqStart, char), sqioSeq(sio), sio->seqLen) ;
	  seqStart += sio->seqLen ;
	  *totSeq += sio->seqLen ;
	  if (sio->mmapSize) setDirectBuf (sio, ti->seq, seqStart) ;
	}
      arrayMax(ti->seq) = seqStart ; // undo the inflation from setting up directBuf
      filled += arrayMax(ti->seq) ;
    } // loading thread
  sio->directBuf = 0 ; sio->directBufSize = 0 ;
  seqIOReleaseRead (sio) ;
  return filled ;
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

// pre-size the kmer hash from the largest input file, so it need not double so often.
// Syncmer density is ~2/(w+1) per base.
// FASTA has ~1.1 bytes/base -> bytes_per_syncmer ~ 0.55*(w+1), use 11*(w+1)/28 for headroom
// FASTQ has ~4 bytes/base   -> bytes_per_syncmer ~ 2*(w+1),    use 11*(w+1)/7 for headroom
// kmerHashAdd() doubles the table anyway if the estimate is too low.
static U64 estimateHashSize (int argc, char **argv, SyncmerParams params)
{
  U64 maxFileSize = 0 ;
  bool isFastq = false ;
  int i ;
  for (i = 0 ; i < argc ; ++i)
    { U64 fs = estimateInputSize (1, &argv[i]) ;
      if (fs > maxFileSize) maxFileSize = fs ;
      int len = strlen (argv[i]) ;
      if ((len > 3 && !strcmp (argv[i]+len-3, ".fq")) ||
	  (len > 6 && (!strcmp (argv[i]+len-6, ".fq.gz") || !strcmp (argv[i]+len-6, ".fastq"))) ||
	  (len > 9 && !strcmp (argv[i]+len-9, ".fastq.gz")))
	isFastq = true ;
    }
  U64 bytesPerSyncmer = 11 * ((U64)params.w + 1) / (isFastq ? 7 : 28) ;
  if (bytesPerSyncmer < 1) bytesPerSyncmer = 1 ;
  U64 estSyncmers = maxFileSize / bytesPerSyncmer ;
  U64 initialSize = estSyncmers * 4 ; // 4x headroom for hash load factor
  U64 plen = (params.w + params.k + 31) >> 5 ; // cap so pack[] does not exceed ~16 GB
  U64 maxPsize = ((U64)16 << 30) / (plen * 8) ;
  if (initialSize > maxPsize * 3)   // psize = size * 0.3, so size = psize/0.3 ~= psize*3
    initialSize = maxPsize * 3 ;
  if (initialSize < (1<<20)) return 0 ; // too small to bother: use the default minimum
  if (maxFileSize > 0)
    fprintf (stdout, "estimated from largest file %llu bytes, pre-sized hash for ~%llu syncmers\n",
	     maxFileSize, estSyncmers) ;
  return initialSize ;
}

/*********** end of file ***********/
