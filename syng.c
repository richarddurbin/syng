/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  6 11:47 2025 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

// AGENDA
// need to add in edges, with counts
// need to consolidate edges to form unitigs
// when do I convert back to non-hoco?
// need some ability to remove nodes or reads, or error-correct

#include <pthread.h>

#include "seqio.h"
#include "seqhash.h"
#include "ONElib.h"

#include "syng.h"

/************************* thread functions ****************************/

typedef struct {
  Seqhash  *sh ;
  KmerHash *kh ;        // read-only here
  OneFile  *ofIn ;
  SyngBWT  *sbwt ;
  I64       nPath ;
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  	                // seq stores just the starts and ends if only those are required
  Array     seqInfo ;	// of SeqInfo: per sequence
  Array     syncPos ;	// of SyncPos: syncmers and their start positions
} ThreadInfo ;

typedef struct {
  I64       len ;	// sequence length
  I64       nSync ;	// number of syncs
} SeqInfo ;

typedef struct {
  I32	    sync ;	// initially 0 if not in kh; -ve if if reverse orientation
  I32       pos ;       // offset in the current sequence
} SyncPos ;

typedef enum { NONE=0, SEQ, PATH, GBWT } OutType ;

static OutType outType = NONE ;

/***** threadProcessRead() handles input from SeqIO: maps sequences to sync,pos *****/

static void *threadProcessRead (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;
  I64 sync ;
  U64 *uBuf = new(ti->kh->plen,U64) ; // working buffer for threadsafe kmerHashFind

  arrayMax(ti->syncPos) = 0 ;
  
  for (i = 0 ; i < arrayMax(ti->seqInfo) ; ++i)
    { I64 seqLen = arrp(ti->seqInfo, i, SeqInfo)->len ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += seqLen ;
      int pos, spStart = arrayMax(ti->syncPos) ;
      SeqhashIterator *sit = syncmerIterator (ti->sh, seq, seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ sync = 0 ; // default if not found
	  kmerHashFindThreadSafe (ti->kh, seq+pos, &sync, uBuf) ;
	  if (sync > 2 || sync < -2 || !sync) // don't record poly-A, poly-C, poly-G, poly-T
	    { SyncPos *sp = arrayp(ti->syncPos, arrayMax(ti->syncPos), SyncPos) ;
	      sp->pos = pos ;
	      sp->sync = sync ;
	    }
	}
      seqhashIteratorDestroy (sit) ;
      arrp(ti->seqInfo, i, SeqInfo)->nSync = arrayMax(ti->syncPos) - spStart ;
    }

  newFree (uBuf, ti->kh->plen, U64) ;
  return 0 ;
}

/***** threadProcessPath() handles input from 1path: maps sync/pos to sequence ******/

static void *threadProcessPath (void* arg) // read in paths, make sequences if necessary
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i, j ;
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
	    for (j = 1 ; j < si->nSync ; ++j)
	      if (!syngBWTpathNext (sbp, &sp[j].sync, &sp[j].pos)) die ("failed GBWT extension") ;
	    syngBWTpathDestroy (sbp) ;
	    break ;
	  case 'z':
	    si->nSync = oneLen(ti->ofIn) ;
	    arrayp(ti->syncPos, spStart+si->nSync, SyncPos)->pos = 0 ;
	    sp = arrp(ti->syncPos, spStart, SyncPos) ; spStart += si->nSync ;
	    I64 *iList = oneIntList(ti->ofIn) ;
	    for (j = 0 ; j < si->nSync ; ++j) sp[j].sync = iList[j] ;
	    break ;
	  case 'o': // comes with z
	    if (oneLen(ti->ofIn) != si->nSync) die ("o error in threadProcessPath %d", i) ;
	    iList = oneIntList(ti->ofIn) ;
	    for (j = 0 ; j < si->nSync ; ++j) sp[j].pos = iList[j] ;
	    break ;
	  case 'X':
	    if (si->nSync && oneLen(ti->ofIn) != sp->pos) die ("X error in threadProcessPath %d", i) ;
	    char *dna = oneDNAchar (ti->ofIn) ;
	    for (j = 0 ; j < oneLen(ti->ofIn) ; ++j) seq[j] = dna2indexConv[dna[j]] ;
	    break ;
	  case 'Y':
	    if (si->nSync && oneLen(ti->ofIn) != si->len - (sp[si->nSync-1].pos + ti->kh->len))
	      die ("Y error in threadProcessPath %d", i) ;
	    dna = oneDNAchar (ti->ofIn) ;
	    int endLen = oneLen(ti->ofIn) ;
	    for (j = 0 ; j < endLen ; ++j) seq[si->len-endLen+j] = dna2indexConv[dna[j]] ;
	    break ;
	  }
      if (outType == SEQ) // build the sequence from the syncs
	{ I64 lastMax = 0 ;
	  for (j = 0 ; j < si->nSync ; ++j)
	    { assert (sp) ;
	      if (!sp[j].sync) continue ;
	      kmerHashSeq (ti->kh, sp[j].sync, seq+sp[j].pos) ;
	    }
	}
    }
  
  return 0 ;
}

/******************* package to handle syncmer parameters *************/

typedef struct {
  int w, k, seed ;
} Params ;

static int PARAMS_K_DEFAULT = 8 ;
static int PARAMS_W_DEFAULT = 55 ;
static int PARAMS_SEED_DEFAULT = 7 ;

static void writeParams (OneFile *of, Params *p)
{
  oneInt(of,0) = p->k ;
  oneInt(of,1) = p->w ;
  oneInt(of,2) = p->seed ;
  oneWriteLine (of, 'h', 0, 0) ;
}

static void readParams (OneFile *of, Params *p)
{ while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType != 'h') die ("sync file %s has no 'h' parameters record", oneFileName(of)) ;
  p->k = oneInt(of,0) ; p->w = oneInt(of,1) ; p->seed = oneInt(of,2) ;
  fprintf (stdout, "read syncmer parameters k %d w %d (size %d) seed %d\n",
	   p->k, p->w, p->k + p->w, p->seed) ;
}

static void checkParams (OneFile *of, Params *p)
{
  if (oneInt(of,0) != p->k ||
      oneInt(of,1) != p->w ||
      oneInt(of,2) != p->seed)
    die ("hash parameters mismatch: (k,w,s) file (%d,%d,%d) != code (%d,%d,%d)",
	 oneInt(of,0), oneInt(of,1), oneInt(of,2), p->k, p->w, p->seed) ;
}

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

/**************** main program ********************/

static char usage[] =
  "Usage: syng <operation>* <input>*\n"
  "possible operations are:\n"
  "  -w <window length>     : [1023] syncmer length = w + k\n"
  "  -k <smer length>       : [16] must be under 32\n"
  "  -seed <seed>           : [7] for the hashing function\n"
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
  char       *kFileName = 0, *kFaFileName = 0 ;
  int         nThread = 8 ;
  pthread_t  *threads ;
  ThreadInfo *threadInfo ;
  KmerHash   *kh = 0, *khNew = 0 ;
  OneSchema  *schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile    *ofK = 0, *ofNewK = 0, *ofOut = 0 ;
  SyngBWT    *gbwtOut = 0 ;
  Params      params ;
  bool        isAddSyncmers = true, isHistK = false, isOutputEnds = false, isOutputNames = false ;
  I64         i, j, k ; // general purpose indices
  
  timeUpdate (0) ;

  params.k = PARAMS_K_DEFAULT ;
  params.w = PARAMS_W_DEFAULT ;
  params.seed = PARAMS_SEED_DEFAULT ;

  U64 nSeq = 0, totSeq = 0, totSync = 0 ;
  
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { fprintf (stderr, "%s",usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-w") && argc > 1) { params.w = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-k") && argc > 1) { params.k = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-seed") && argc > 1) { params.seed = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-T") && argc > 1) { nThread = atoi(argv[1]) ; argc -=2 ; argv +=2 ; }
    else if (!strcmp (*argv, "-readK") && argc > 1)
      { kFileName = argv[1] ; argc -= 2 ; argv +=2 ;
	OneFile *of = oneFileOpenRead (kFileName, schema, "khash", 1) ;
	if (!of) die ("failed to open syncmer file %s to read", kFileName) ;

	readParams (of, &params) ;             // read the syncmer hash parameters
	kh = kmerHashReadOneFile (of) ;        // read in the kmerhash
	if (kh->len != params.w + params.k)
	  die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;
	for (i = 1 ; i <= kmerHashMax(kh) ; ++i) totSync += arr(kh->count,i,I64) ;
	fprintf (stdout, "read %llu syncmers from %s with total count %llu\n",
		kmerHashMax(kh), kFileName, totSync) ;
	oneFileClose (of) ; // close the old file
	timeUpdate (stdout) ;
      }
    else if (!strcmp (*argv, "-zeroK"))
      { if (!kh) die ("-zeroK error: no khash file read in yet") ;
	memset (arrp(kh->count,0,I64), 0, arrayMax(kh->count)*sizeof(I64)) ;
	argc-- ; argv++ ; 
      }
    else if (!strcmp (*argv, "-limitK") && argc > 2)
      { if (!kh) die ("-limitK error: no khash file read in yet") ;
	int min = atoi (argv[1]), max = atoi(argv[2]) ;
	I64 n = 0, *c = arrp(kh->count,1,I64) ;
	for (i = 1 ; i <= kh->max ; ++i, ++c) if (*c >= min && (*c <= max || !max)) ++n ;
	KmerHash *kh2 = kmerHashCreate (n*4, kh->len) ;
	c = arrp(kh->count,1,I64) ;
	for (i = 1 ; i <= kh->max ; ++i, ++c)
	  if (*c >= min && (*c <= max || !max)) kmerHashAddPacked (kh2, kh->pack + i*kh->plen, 0) ;
	kmerHashDestroy (kh) ;
	kh = kh2 ;
	argc -= 3 ; argv += 3 ;
      }
    else if (!strcmp (*argv, "-histK")) { isHistK = true ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-noAddK"))   { isAddSyncmers = false ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-writeK"))
      { char *fname = fnameTag (outPrefix,"1khash") ;
	if (kFileName)
	  { OneFile *ofIn = oneFileOpenRead (kFileName, 0, "khash", 1) ;
	    ofK = oneFileOpenWriteFrom (fname, ofIn, true, 1) ;
	    oneFileClose (ofIn) ;
	  }
	else
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
	if (kh && kh->max)
	  gbwtOut = syngBWTcreate (params.k + params.w, kh->max + 1) ;
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
      { khNew = kmerHashCreate (0, params.w + params.k) ;
	ofNewK = oneFileOpenWriteNew (fnameTag(argv[1],"1khash"), schema, "khash", true, 1) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-outputEnds")) { isOutputEnds = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-outputNames")) { isOutputNames = true ; --argc ; ++argv ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  fprintf (stdout, "k, w, seed are %d %d %d\n", params.k, params.w, params.seed) ;
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly

  if (!kh)
    { kh = kmerHashCreate (0, params.w + params.k) ;
      char *s = new0 (kh->len, char) ;
      for (i = 0 ; i < 4 ; ++i)
	{ for (j = 0 ; j < kh->len ; ++j) s[j] = i ;
	  I64 sync ; kmerHashAdd (kh, s, &sync) ; // 0,1,2,3 map to 1,2,-2,-1
	}
      newFree (s, kh->len, char) ;
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
  
  U64 nSeq0 = 0, totSeq0 = 0, totSync0 = totSync, syncMax0 = kmerHashMax(kh), nNew = 0 ;
  
  if (nThread < 1) die ("number of threads %d must be at least 1", nThread) ;
  threads = new (nThread, pthread_t) ;
  threadInfo = new0 (nThread, ThreadInfo) ;
  for (i = 0 ; i < nThread ; ++i)
    { threadInfo[i].sh = sh ;
      threadInfo[i].kh = kh ;
      threadInfo[i].seq = arrayCreate (101<<20, char) ; // 101 Mb
      threadInfo[i].seqInfo = arrayCreate (20000, SeqInfo) ;
      threadInfo[i].syncPos = arrayCreate (1<<20, SyncPos) ;
    }

  int nSource = 0 ;
  while (argc--)
    { ++nSource ;
      OneFile *ofIn = oneFileOpenRead (*argv, schema, 0, nThread) ;
      SeqIO   *sio = 0 ;
      I64      nPath ;
      if (ofIn && oneStats (ofIn, 'P', &nPath, 0, 0) && nPath) // a path file 
	{ if (!kh->max) die ("%s is a path file but we have no kmers", *argv) ; // must have syncmers
	  fprintf (stdout, "path file %d %s: ", nSource, *argv) ;
	  oneReadLine (ofIn) ; if (ofIn->lineType == 'h') checkParams (ofIn, &params) ;
	  I64 z ;
	  if (oneStats (ofIn, 'Z', &z, 0, 0) && z) // must have a GBWT
	    threadInfo->sbwt = syngBWTread (ofIn) ;
	  for (i = 0 ; i < nThread ; ++i)
	    { threadInfo[i].ofIn = ofIn+i ;
	      oneGoto (ofIn+i, 'P', 1 + (nPath * i) / nThread) ;
	      threadInfo[i].nPath = (nPath*(i+1))/nThread - (nPath*i)/nThread ;
	      threadInfo[i].sbwt = threadInfo->sbwt ; // copy from the master
	    }
	}
      else
	{ if (ofIn) { oneFileClose (ofIn) ; ofIn = 0 ; }
	  sio = seqIOopenRead (*argv, dna2indexConv, 0) ;
	  fprintf (stdout, "sequence file %d %s type %s: ", nSource, *argv, seqIOtypeName[sio->type]) ;
	}
      if (!ofIn && !sio)
	die ("failed to open sequence or path file %s", *argv) ;
      fflush (stdout) ;
      bool isDone = false ;
      while (!isDone) // read this file
	{ if (sio)
	    { for (i = 0 ; i < nThread ; ++i)  // read 100Mb DNA per thread and find kmers in parallel
		{ ThreadInfo *ti = &threadInfo[i] ;
		  arrayMax(ti->seq) = 0 ;
		  arrayMax(ti->seqInfo) = 0 ;
		  int seqStart = 0 ;
		  while (arrayMax(ti->seq) < 100<<20 && seqIOread (sio))
		    { arrayp(ti->seqInfo, arrayMax(ti->seqInfo), SeqInfo)->len = sio->seqLen ;
		      array(ti->seq, seqStart+sio->seqLen, char) = 0 ;
		      memcpy (arrp(ti->seq, seqStart, char), sqioSeq(sio), sio->seqLen) ;
		      seqStart += sio->seqLen ;
		      totSeq += sio->seqLen ;
		    }
		  if (!arrayMax(ti->seq)) isDone = true ; // we are done after processing this lot
		} // loading thread
	      for (i = 0 ; i < nThread ; ++i) // create threads
		pthread_create (&threads[i], 0, threadProcessRead, &threadInfo[i]) ;
	      for (i = 0 ; i < nThread ; ++i)
		pthread_join (threads[i], 0) ; // wait for threads to complete
	      for (i = 0 ; i < nThread ; ++i) // must go through threads linearly to add missing syncs
		{ ThreadInfo *ti = threadInfo + i ;
		  SyncPos *sp = arrp(ti->syncPos, 0, SyncPos) ;
		  char *seq = arrp(ti->seq, 0, char) ;
		  totSync += arrayMax(ti->syncPos) ;
		  for (j = 0 ; j < arrayMax(ti->seqInfo) ; ++j)
		    { for (k = 0 ; k < arrp(ti->seqInfo, j, SeqInfo)->nSync ; ++k, ++sp)
			if (sp->sync) // increment kh->count, because FindThreadSafe could not
			  ++array(kh->count,(sp->sync >= 0) ? sp->sync : -sp->sync,I64) ;
			else if (isAddSyncmers)
			  { I64 sync ;
			    kmerHashAdd (kh, seq + sp->pos, &sync) ;
			    sp->sync = sync ;
			  }
			else if (khNew)
			  kmerHashAdd (khNew, seq + sp->pos, 0) ;
		      seq += arrp(ti->seqInfo, j, SeqInfo)->len ;
		    } // read
		} // thread i
	    } // sio - reading a sequence file
	  else if (ofIn)
	    { for (i = 0 ; i < nThread ; ++i) // create threads
		pthread_create (&threads[i], 0, threadProcessPath, &threadInfo[i]) ;
	      for (i = 0 ; i < nThread ; ++i)
		pthread_join (threads[i], 0) ; // wait for threads to complete
	      for (i = 0 ; i < nThread ; ++i) // must go through threads linearly to count syncs
		{ ThreadInfo *ti = threadInfo + i ;
		  totSync += arrayMax(ti->syncPos) ;
		  SyncPos *sp = arrp(ti->syncPos, 0, SyncPos) ;
		  for (j = 0 ; j < arrayMax(ti->seqInfo) ; ++j)
		    for (k = 0 ; k < arrp(ti->seqInfo, j, SeqInfo)->nSync ; ++k, ++sp)
		      if (sp->sync) // need to increment kh->count, because FindThreadSafe could not
			++array(kh->count,((sp->sync >= 0) ? sp->sync : -sp->sync),I64) ;
		} // thread i
	      isDone = true ; // do each file in one go
	    }

	  // now the output
	  
	  for (i = 0 ; i < nThread ; ++i) // go through threads and sequences within threads
	    { ThreadInfo *ti = threadInfo + i ;
	      char *seq = arrp(ti->seq, 0, char) ;
	      SyncPos *sp = arrp(ti->syncPos, 0, SyncPos) ;
	      for (j = 0 ; j < arrayMax (ti->seqInfo) ; ++j, ++nSeq)
		{ if (outType == SEQ)
		    oneWriteLine (ofOut, 'S', arrp(ti->seqInfo, j, SeqInfo)->len, seq) ;
		  else if (outType == PATH || outType == GBWT)
		    { I64 nSync = arrp(ti->seqInfo, j, SeqInfo)->nSync ; // number of syncs
		      oneInt(ofOut, 0) = arrp(ti->seqInfo, j, SeqInfo)->len ;
		      oneInt(ofOut, 1) = nSource ; // first write the path number
		      oneInt(ofOut, 2) = nSeq-nSeq0+1 ;
		      oneWriteLine (ofOut, 'P', 0, 0) ;
		      if (nSync && outType == GBWT) // add paths to the GBWT and write the start nodes
			{ SyngBWTpath *sbp = syngBWTpathStartNew (gbwtOut, sp->sync) ;
			  oneInt(ofOut, 0) = sp->sync ;
			  oneInt(ofOut, 1) = sp->pos ;
			  oneInt(ofOut, 2) = sbp->inCount ;
			  oneInt(ofOut, 3) = nSync ;
			  oneWriteLine (ofOut, 'Z', 0, 0) ;
			  for (k = 1 ; k < nSync ; ++k)
			    syngBWTpathAdd (sbp, sp[k].sync, sp[k].pos - sp[k-1].pos) ;
			  syngBWTpathFinish (sbp) ;
			  // now add the reverse path
			  sbp = syngBWTpathStartNew (gbwtOut, -sp[nSync-1].sync) ;
			  for (k = nSync-2 ; k >= 0 ; --k)
			    syngBWTpathAdd (sbp, -sp[k].sync, sp[k+1].pos - sp[k].pos) ;
			  syngBWTpathFinish (sbp) ;
			}
		      else if (nSync && outType == PATH)
			{ static I64 *x = 0 ; // memory leak here...
			  static size_t xSize = 0 ;
			  if (!x)
			    { xSize = nSync ; x = new (xSize, I64) ; }
			  else if (xSize < nSync)
			    { newFree (x,xSize,I64) ; xSize = nSync ; x = new (xSize, I64) ; }
			  for (k = 0 ; k < nSync ; ++k) x[k] = sp[k].sync ;
			  oneWriteLine (ofOut, 'z', nSync, x) ;
			  for (k = 0 ; k < nSync ; ++k) x[k] = sp[k].pos ;
			  oneWriteLine (ofOut, 'o', nSync, x) ;
			}
		      if (isOutputEnds)
			{ I64 len = arrp(ti->seqInfo,j,SeqInfo)->len ;
			  if (nSync)
			    { oneWriteLine (ofOut, 'X', sp->pos, seq) ;
			      I64 endOff = sp[nSync-1].pos + kh->len ;
			      oneWriteLine (ofOut, 'Y', len - endOff, seq + endOff) ;
			    }
			  else // split sequence into two - both parts must be less than kh->len
			    { oneWriteLine (ofOut, 'X', len/2, seq) ;
			      oneWriteLine (ofOut, 'Y', len - len/2, seq + len/2) ;
			    }
			}
		      sp += nSync ;
		    } // PATH or GBWT
		  seq +=  arrp(ti->seqInfo, j, SeqInfo)->len ;
		} // j sequences
	    } // i threads
	} // isDone: end of file
      
      if (sio)
	{ seqIOclose (sio) ;
	  fprintf (stdout, "had %llu sequences %llu bp, ",
		   nSeq-nSeq0, totSeq-totSeq0) ;
	  fprintf (stdout, "yielding %llu syncs with %llu extra syncmers\n",
		   totSync-totSync0, kmerHashMax(kh)-syncMax0) ;
	}
      else if (ofIn)
	{ oneFileClose (ofIn) ;
	  fprintf (stdout, "had %llu sequences containing %llu syncmers\n",
		   nSeq-nSeq0, totSync-totSync0) ;
	}
      timeUpdate (stdout) ;

      nSeq0 = nSeq ; totSeq0 = totSeq ; totSync0 = totSync ; syncMax0 = kmerHashMax(kh) ;
      ++argv ;
    } // source file

  fprintf (stdout, "Total for this run %llu sequences, total length %llu\n", nSeq, totSeq) ;
  fprintf (stdout, "Overall total %llu instances of %llu syncmers, average %.2f coverage\n", 
	   totSync, kmerHashMax(kh), totSync / (double)kmerHashMax(kh)) ;

  // destroy the thread objects
  for (i = 0 ; i < nThread ; ++i)
    { arrayDestroy (threadInfo[i].seq) ;
      arrayDestroy (threadInfo[i].seqInfo) ;
      arrayDestroy (threadInfo[i].syncPos) ;
    }
  newFree (threads, nThread, pthread_t) ;
  newFree (threadInfo, nThread, ThreadInfo) ;

  if (isHistK) // histogram of kmers
    { Array hist = arrayCreate (1024, I64) ;
      for (i = 1 ; i < arrayMax(kh->count) ; ++i) ++array(hist,arr(kh->count,i,I64),I64) ;
      qhist(hist,stdout) ;
      arrayDestroy (hist) ;
#ifdef COUNT_POLY_G
      I64 *cgCount = new0 (kh->len, I64) ;
      int nBytes = (kh->len+3) / 4 ;
      int byteCountC[256], byteCountG[256] ;
      for (i = 0 ; i < 256 ; ++i)
	{ byteCountC[i] = ((i&3) == 1) + ((i&12) == 4) + ((i&48) == 16) + ((i&192) == 64) ;
	  byteCountG[i] = ((i&3) == 2) + ((i&12) == 8) + ((i&48) == 32) + ((i&192) == 128) ;
	}
      for (i = 1 ; i <= kh->max ; ++i)
	{ I64 nC = 0, nG = 0 ;
	  U8 *x ;
	  for (j = 0, x = (U8*)(kh->pack + i*kh->plen) ; j < nBytes ; ++j, ++x)
	    { nC += byteCountC[(int)*x] ; nG += byteCountG[(int)*x] ; }
	  if (nC > nG) ++cgCount[nC] ; else ++cgCount[nG] ;
	}
      for (i = 0 ; i < kh->len ; ++i)	printf ("CG %2lld %lld\n", i, cgCount[i]) ;
      newFree (cgCount, kh->len, I64) ;
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
  
  if (ofK) // write the kmerhash as ONEcode
    { writeParams (ofK, &params) ;
      kmerHashWriteOneFile (kh, ofK) ;
      fprintf (stdout, "wrote %llu syncmers to file %s\n", kmerHashMax(kh), oneFileName(ofK)) ;
      oneFileClose (ofK) ;
      timeUpdate (stdout) ;
    }
  if (ofNewK) // write the new kmerhash
    { writeParams (ofNewK, &params) ;
      kmerHashWriteOneFile (khNew, ofNewK) ;
      fprintf (stdout, "wrote %llu syncmers to file %s\n", kmerHashMax(khNew), oneFileName(ofNewK)) ;
      oneFileClose (ofNewK) ;
      timeUpdate (stdout) ;
    }
  if (kFaFileName)
    { char  *fname = fnameTag (kFaFileName,"1kmer.fa.gz") ;
      SeqIO *sio = seqIOopenWrite (fname, FASTA, 0, 0) ;
      if (!sio) die ("failed to open fasta file %s to write", fname) ;
      char idBuf[24], *sBuf = new0 (kh->len+1, char) ;
      for (i = 1 ; i <= kh->max ; ++i)
	{ sprintf (idBuf, "%lld", i) ;
	  seqIOwrite (sio, idBuf, 0, kh->len, kmerHashSeq (kh, i, 0), sBuf) ;
	}
      newFree (sBuf, kh->len+1, char) ;
      seqIOclose (sio) ;
      fprintf (stdout, "wrote %llu syncmers to file %s\n", kmerHashMax(kh), fname) ;
      timeUpdate (stdout) ;
    }

  if (kh) kmerHashDestroy (kh) ;
  if (khNew) kmerHashDestroy (khNew) ;
  seqhashDestroy (sh) ;
	  
  fprintf (stdout, "total: ") ; timeTotal (stdout) ;
  return 0 ;
}

/*********************** end of file **********************/
