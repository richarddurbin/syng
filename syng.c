/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  9 20:37 2024 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

#define SYNG_VERSION "1.0"

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
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  Array     seqLen ;	// of I64, input: per sequence
  Array     pos ;	// of I64, output: concatenated start positions of syncmers
  Array     posLen ;	// of I64, output: number of syncmers per sequence
  Array     sync ;      // of I64, if non-zero, found in kh; -ve if reverse direction
} ThreadInfo ;

/***** threadProcessRead() handles input from SeqIO: maps sequences to sync,pos *****/

static void *threadProcessRead (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;
  I64 iSync ;
  U64 *uBuf = new(ti->kh->plen,U64) ;

  arrayMax(ti->pos) = 0 ;
  arrayMax(ti->posLen) = 0 ;
  
  for (i = 0 ; i < arrayMax(ti->seqLen) ; ++i)
    { I64 seqLen = arr(ti->seqLen, i, I64) ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += seqLen ;
      int pos, posStart = arrayMax(ti->pos) ;
      SeqhashIterator *sit = syncmerIterator (ti->sh, seq, seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ I64 iPos = arrayMax(ti->pos) ;
	  array(ti->pos, iPos, I64) = pos ;
	  iSync = 0 ;
	  kmerHashFindThreadSafe (ti->kh, seq+pos, &iSync, uBuf) ;
	  array(ti->sync,iPos,I64) = iSync ; // will be 0 if not found
	}
      seqhashIteratorDestroy (sit) ;
      array(ti->posLen, i, I64) = arrayMax(ti->pos) - posStart ;
    }

  newFree (uBuf, ti->kh->plen, U64) ;
  return 0 ;
}

/***** threadProcessPath() handles input from 1path: maps sync/pos to sequence ******/


static void *threadProcessPath (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;
  I64 iSync ;
  U64 *uBuf = new(ti->kh->plen,U64) ;
  return 0 ;
}

/******************* package to handle syncmer parameters *************/

typedef struct {
  int w, k, seed ;
} Params ;

static int PARAMS_K_DEFAULT = 16 ;
static int PARAMS_W_DEFAULT = 1023 ;
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
  fprintf (stderr, "read syncmer parameters k %d w %d (size %d) seed %d\n",
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
  "  -limitK <min> <max>    : filter current kmer set based on counts ; max 0 to have no upper bound\n"
  "  -histK                 : output quadratic histogram of kmer counts (after sequence processsing)\n"
  "  -writeK                : write the syncmers as a .1khash file\n"
  "  -writeKfa              : write the syncmers as a fasta file, with ending .kmer.fa.gz\n"
  "  -writeNewK <file prefix> : write new syncmers as a .1khash file; implies -N\n"
  "  -noNewK                : do not create new syncmers - convert unmatched syncmers to 0\n"
  "  -writePath             : write a .1path file (paths of nodes)\n"
  "  -writeGBWT             : write a .1gbwt file (nodes, edges and paths in GBWT form)\n"
  "  -writeEnds             : write the non-syncmer ends of sequences to 1path or 1gbwt files\n"
  "  -writeSeq              : write a .1seq file (paths converted back to sequences)\n"
  "  -newSummary            : summary of the new data set\n"
  "  -filterG <nG>          : remove input sequences with runs of >nG consecutive Gs (bad Illumina reads)\n"
  "  -filterQ <QT>          : remove input sequences with average quality score below QT\n"
  "  -filterIllumina        : equivalent to -filterPolyG 60 -filterQual 30\n"
  "possible inputs are:\n"
  "  <sequence file>        : any of fasta[.gz], fastq[.gz], BAM/CRAM/SAM, .1seq\n"
  "  <.1path file>          : sequences as lists of kmers, with optional non-syncmer DNA ends\n"
  "  <.1gbwt file>          : graph BWT with paths, with optional ends\n"
  "Operations are carried out in order as they are parsed, with some setting up future actions,\n"
  "e.g. changing the outfile prefix affects following lower case options for file opening\n"
  "Some output files, e.g. .1gbwt will be output at the end, after all inputs are processed,\n"
  "whereas others, e.g. .1path are written as inputs are processed.\n" ;

typedef enum { NONE=0, SEQ, PATH, GBWT } OutType ;

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
  OutType     outType = NONE ;
  SyngBWT    *gbwtOut = 0 ;
  Params      params ;
  bool        isAddSyncmers = true, isHistK = false, isWriteEnds = false ;
  int         filterG = 0, filterQ = 0 ;
  I64         i, j, k ; // general purpose indices
  
  timeUpdate (0) ;

  params.k = PARAMS_K_DEFAULT ;
  params.w = PARAMS_W_DEFAULT ;
  params.seed = PARAMS_SEED_DEFAULT ;

  U64 nSeq = 0, nFilt = 0, totSeq = 0, totSync = 0 ;
  
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
	fprintf (stderr, "read %llu syncmers from %s with total count %llu\n",
		kmerHashMax(kh), kFileName, totSync) ;
	oneFileClose (of) ; // close the old file
	timeUpdate (stderr) ;
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
    else if (!strcmp (*argv, "-noNewK"))   { isAddSyncmers = false ; argc-- ; argv++ ; }
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
	gbwtOut = syngBWTcreate (params.k + params.w) ;
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
    else if (!strcmp (*argv, "-writeEnds")) { isWriteEnds = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-filterG") && argc > 1)
      { filterG = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-filterQ") && argc > 1)
      { filterQ = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-filterIllumina")) { filterG = 60 ; filterQ = 20 ; --argc ; ++argv ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  fprintf (stderr, "k, w, seed are %d %d %d\n", params.k, params.w, params.seed) ;
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly

  if (!kh) kh = kmerHashCreate (0, params.w + params.k) ; 

  if (ofOut) oneAddProvenance (ofOut, "syng", SYNG_VERSION, getCommandLine()) ;

  int maxSources ; // collect sources from remaining files
  char **sources = collectSources (argc, argv, &maxSources) ;
  if (sources)
    { if (ofK) addSourceReferences (ofK, sources) ;
      if (ofOut) addSourceReferences (ofOut, sources) ;
      char **s = sources ; while (*s) { free (*s) ; ++s ; }
      newFree (sources, maxSources, char*) ;
    }
  
  U64 nSeq0 = 0, nFilt0 = 0, totSeq0 = 0, totSync0 = totSync, syncMax0 = kmerHashMax(kh), nNew = 0 ;

  Array seqDir = 0 ;
  if (ofOut) seqDir = arrayCreate(1024,char) ; // might be needed if we will flip
  
  if (nThread < 1) die ("number of threads %d must be at least 1", nThread) ;
  threads = new (nThread, pthread_t) ;
  threadInfo = new0 (nThread, ThreadInfo) ;
  for (i = 0 ; i < nThread ; ++i)
    { threadInfo[i].sh = sh ;
      threadInfo[i].kh = kh ;
      threadInfo[i].seq = arrayCreate (101<<20, char) ; // 101 Mb
      threadInfo[i].seqLen = arrayCreate (20000, I64) ;
      threadInfo[i].pos = arrayCreate (1<<20, I64) ;
      threadInfo[i].sync = arrayCreate (1<<20, I64) ;
      threadInfo[i].posLen = arrayCreate (20000, I64) ;
    }

  int nSource = 0 ;
  while (argc--)
    { ++nSource ;
      OneFile *ofIn = oneFileOpenRead (*argv, schema, 0, nThread) ;
      SeqIO   *sio = 0 ;
      I64      nPath ;
      if (ofIn && oneStats (ofIn, 'P', &nPath, 0, 0) && nPath) // a path file 
	{ if (!kh->max) die ("%s is a path file but we have no kmers", *argv) ; // must have syncmers
	  fprintf (stderr, "path file %d %s: ", nSource, *argv) ;
	  for (i = 0 ; i < nThread ; ++i)
	    { threadInfo[i].ofIn = ofIn+i ;
	      oneGoto (ofIn+i, 'P', 1 + nPath * i / nThread) ;
	      
	    }
	}
      else
	{ if (ofIn) { oneFileClose (ofIn) ; ofIn = 0 ; }
	  sio = seqIOopenRead (*argv, dna2indexConv, (filterQ > 0)) ;
	  fprintf (stderr, "sequence file %d %s type %s: ", nSource, *argv, seqIOtypeName[sio->type]) ;
	}
      if (!ofIn && !sio)
	die ("failed to open sequence or path file %s", *argv) ;
      if (filterQ > 0 && (!sio || !sio->isQual))
	die ("\nfile %s has no qualities, but filterQ is set", *argv) ;
      fflush (stderr) ;
      bool isDone = false ;
      while (!isDone) // read this file
	{ if (sio)
	    { for (i = 0 ; i < nThread ; ++i)  // read 100Mb DNA per thread and find kmers in parallel
		{ ThreadInfo *ti = &threadInfo[i] ;
		  arrayMax(ti->seq) = 0 ;
		  arrayMax(ti->seqLen) = 0 ;
		  int seqStart = 0 ;
		  while (arrayMax(ti->seq) < 100<<20 && seqIOread (sio))
		    { if (filterG)
			{ int nG = 0 ;
			  char *s = sqioSeq(sio) ;
			  for (i = 0 ; i < sio->seqLen ; ++i)
			    if (*s++ != 2) nG = 0 ;
			    else if (++nG > filterG) { sio->seqLen = 0 ; ++nFilt ; break ; }
			}
		      if (filterQ)
			{ I64 totQ = 0 ;
			  char *q = sqioQual(sio) ;
			  for (i = 0 ; i < sio->seqLen ; ++i) totQ += *q++ ;
			  if (totQ < filterQ * sio->seqLen) { sio->seqLen = 0 ; ++nFilt ; }
			}
		      array(ti->seqLen, arrayMax(ti->seqLen), I64) = sio->seqLen ;
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
		  I64 *pos = arrp(ti->pos, 0, I64) ;
		  I64 *sync = arrp(ti->sync, 0, I64) ;
		  char *seq = arrp(ti->seq, 0, char) ;
		  totSync += arrayMax(ti->sync) ;
		  for (j = 0 ; j < arrayMax(ti->posLen) ; ++j)
		    { for (k = 0 ; k < arr(ti->posLen, j, I64) ; ++k, ++sync, ++pos)
			if (*sync) // need to increment kh->count, because FindThreadSafe could not
			  ++array(kh->count,((*sync >= 0) ? *sync : -*sync),I64) ;
			else if (isAddSyncmers)
			  kmerHashAdd (kh, seq + *pos, sync) ;
			else if (khNew)
			  kmerHashAdd (khNew, seq + *pos, 0) ;
		      seq += arr(ti->seqLen, j, I64) ;
		    } // read
		} // thread i
	    } // sio - reading a sequence file
	  else if (ofIn)
	    { for (i = 0 ; i < nThread ; ++i)  // reconstruct 100Mb-worth of DNA or paths per thread
		{ ThreadInfo *ti = &threadInfo[i] ;
		  arrayMax(ti->seq) = 0 ;
		  arrayMax(ti->seqLen) = 0 ;
		} // loading thread
	      for (i = 0 ; i < nThread ; ++i) // create threads
		pthread_create (&threads[i], 0, threadProcessPath, &threadInfo[i]) ;
	      for (i = 0 ; i < nThread ; ++i)
		pthread_join (threads[i], 0) ; // wait for threads to complete
	      for (i = 0 ; i < nThread ; ++i) // must go through threads linearly to count syncs
		{ ThreadInfo *ti = threadInfo + i ;
		  totSync += arrayMax(ti->sync) ;
		  I64 *sync = arrp(ti->sync, 0, I64) ;
		  for (j = 0 ; j < arrayMax(ti->posLen) ; ++j)
		    for (k = 0 ; k < arr(ti->posLen, j, I64) ; ++k, ++sync)
		      if (*sync) // need to increment kh->count, because FindThreadSafe could not
			++array(kh->count,((*sync >= 0) ? *sync : -*sync),I64) ;
		} // thread i
	    }

	  // now the output
	  
	  for (i = 0 ; i < nThread ; ++i) // go through threads and sequences within threads
	    { ThreadInfo *ti = threadInfo + i ;
	      char *seq = arrp(ti->seq, 0, char) ;
	      I64  *posList = arrp(ti->pos, 0, I64) ;
	      I64  *syncList = arrp(ti->sync, 0 , I64) ;
	      for (j = 0 ; j < arrayMax (ti->seqLen) ; ++j, ++nSeq)
		if (outType == SEQ)
		  { oneWriteLine (ofOut, 'S', arr(ti->seqLen, j, I64), seq) ;
		    seq += arr(ti->seqLen, j, I64) ;
		  }
		else if (outType == PATH || outType == GBWT)
		  { I64 posLen = arr(ti->posLen, j, I64) ; // number of syncs
		    oneInt(ofOut, 0) = nSource ; // first write the path number
		    oneInt(ofOut, 1) = nSeq-nSeq0+1 ;
		    oneWriteLine (ofOut, 'P', 0, 0) ;
		    if (isWriteEnds) // and its ends if they were requested
		      { oneWriteLine (ofOut, 'X', posList[0], seq) ;
			if (posLen)
			  oneWriteLine (ofOut, 'Y',
					arr(ti->seqLen,j,I64) - posList[posLen-1] - kh->len,
					seq+posList[posLen-1]) ;
		      }
		    if (posLen && outType == GBWT) // add paths to the GBWT and write the start nodes
		      { SyngBWTpath *sbp = syngBWTpathStart (gbwtOut, syncList[0]) ;
			oneInt(ofOut, 0) = syncList[0] ;
			oneInt(ofOut, 1) = sbp->startCount ;
			oneInt(ofOut, 2) = posLen ;
			oneWriteLine (ofOut, 'Z', 0, 0) ;
			for (k = 1 ; k < posLen ; ++k)
			  syngBWTpathAdd (sbp, syncList[k], posList[k] - posList[k-1]) ;
			syngBWTpathFinish (sbp) ;
			// now add the reverse path
			sbp = syngBWTpathStart (gbwtOut, -syncList[posLen-1]) ;
			for (k = posLen-2 ; k >= 0 ; --k)
			  syngBWTpathAdd (sbp, -syncList[k], posList[k+1] - posList[k]) ;
			syngBWTpathFinish (sbp) ;
		      }
		    else if (posLen && outType == PATH) // second because it can corrupt syncList[]
		      { // test if worth flipping sites
			bool isFlip = false ;
			I64 *sync = syncList, sync0 = *sync++ ;
			for (k = 1 ; k < posLen ; ++k, ++sync)
			  if (*sync > sync0-64 && *sync < sync0+64) sync0 = *sync ;
			  else if (-*sync > sync0-64 && -*sync < sync0+64)
			    { isFlip = true ; sync0 = -*sync ; }
			  else
			    { isFlip = false ; break ; }
			if (isFlip)
			  { arrayMax(seqDir) = 0 ;
			    array(seqDir,posLen,char) = 0 ;
			    for (k = 0 ; k < posLen ; ++k)
			      if (syncList[k] >= 0)
				arr(seqDir,k,char) = '+' ;
			      else
				{ arr(seqDir,k,char) = '-' ;
				  syncList[k] = -syncList[k] ;
				}
			  }
			oneWriteLine (ofOut, 'z', posLen, syncList) ;
			oneWriteLine (ofOut, 'o', posLen, posList) ;
			if (isFlip) oneWriteLine (ofOut, 'd', posLen, arrp(seqDir,0,char)) ;
		      }
		    posList += posLen ;
		    syncList += posLen ;
		  } // posLen
	    } // thread
	} // isDone: end of file
      
      if (sio)
	{ seqIOclose (sio) ;
	  fprintf (stderr, "had %llu sequences (%lld filtered) %llu bp, ",
		   nSeq-nSeq0, nFilt-nFilt0, totSeq-totSeq0) ;
	  fprintf (stderr, "yielding %llu syncs with %llu extra syncmers\n",
		   totSync-totSync0, kmerHashMax(kh)-syncMax0) ;
	}
      else if (ofIn)
	{ oneFileClose (ofIn) ;
	  fprintf (stderr, "had %llu sequences containing %llu syncmers\n",
		   nSeq-nSeq0, totSync-totSync0) ;
	}
      timeUpdate (stderr) ;

      nSeq0 = nSeq ; nFilt0 = nFilt ; totSeq0 = totSeq ; totSync0 = totSync ;
      syncMax0 = kmerHashMax(kh) ;
      ++argv ;
    } // source file

  fprintf (stderr, "Total for this run %llu sequences, total length %llu\n", nSeq, totSeq) ;
  fprintf (stderr, "Overall total %llu instances of %llu syncmers, average %.2f coverage\n", 
	   totSync, kmerHashMax(kh), totSync / (double)kmerHashMax(kh)) ; 

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
      timeUpdate (stderr) ;
    }

  // now write all the files requested

  if (outType == GBWT)
    { syngBWTwrite (ofOut, gbwtOut) ;
      fprintf (stderr, "wrote gbwt to file %s\n", oneFileName(ofOut)) ;
      syngBWTdestroy (gbwtOut) ;
      timeUpdate (stderr) ;
    }
  if (ofOut) oneFileClose (ofOut) ;
  
  if (ofK) // write the kmerhash as ONEcode
    { writeParams (ofK, &params) ;
      kmerHashWriteOneFile (kh, ofK) ;
      fprintf (stderr, "wrote %llu syncmers to file %s\n", kmerHashMax(kh), oneFileName(ofK)) ;
      oneFileClose (ofK) ;
      timeUpdate (stderr) ;
    }
  if (ofNewK) // write the new kmerhash
    { writeParams (ofNewK, &params) ;
      kmerHashWriteOneFile (khNew, ofNewK) ;
      fprintf (stderr, "wrote %llu syncmers to file %s\n", kmerHashMax(khNew), oneFileName(ofNewK)) ;
      oneFileClose (ofNewK) ;
      timeUpdate (stderr) ;
    }
  if (kFaFileName)
    { char  *fname = fnameTag (kFaFileName,"1kmer.fa.gz") ;
      SeqIO *sio = seqIOopenWrite (fname, FASTA, 0, 0) ;
      if (!sio) die ("failed to open fasta file %s to write", fname) ;
      char   idBuf[24] ;
      for (i = 1 ; i <= kh->max ; ++i)
	{ sprintf (idBuf, "%lld", i) ;
	  seqIOwrite (sio, idBuf, 0, kh->len, kmerHashSeq (kh, i), 0) ;
	}
      seqIOclose (sio) ;
      fprintf (stderr, "wrote %llu syncmers to file %s\n", kmerHashMax(kh), fname) ;
      timeUpdate (stderr) ;
    }

  if (kh) kmerHashDestroy (kh) ;
  if (khNew) kmerHashDestroy (khNew) ;
  seqhashDestroy (sh) ;
  if (seqDir) arrayDestroy (seqDir) ;
	  
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
