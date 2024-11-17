/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 16 11:33 2024 (rd109)
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

typedef struct {
  Seqhash  *sh ;
  KmerHash *kh ;        // read-only here
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  Array     seqLen ;	// of I64, input: per sequence
  Array     pos ;	// of I64, output: concatenated start positions of syncmers
  Array     posLen ;	// of I64, output: number of syncmers per sequence
  Array     sync ;      // of I64, if non-zero, found in kh; -ve if reverse direction
} ThreadInfo ;

static void *threadProcess (void* arg) // find the start positions of all the syncmers
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
  "  -readK <.1khash file>  : read and start from this syncmer (khash) file\n"
  "  -zeroK                 : zero the kmer counts\n"
  "  -noNewK                : do not create new syncmers - convert unmatched syncmers to 0\n"
  "  -o <outfile prefix>    : [syngOut] applies to all write* options\n"
  "  -writeK                : write the syncmers as a .1khash file\n"
  "  -writeKfa              : write the syncmers as a fasta file, with ending .kmer.fa.gz\n"
  "  -histK                 : write a quadratic histogram of kmer counts\n"
  "  -writeNewK             : write a .nk binary file of unmatched syncmers for sorting; implies -N\n"
  "  -checkNewK             : deep dive into newK\n"
  "  -histNewK              : write a quadratic histogram of new kmers; implies -N\n"
  "  -writePath             : write a .1path file (paths of nodes)\n"
  "  -writeGBWT             : write a .1gbwt file (nodes, edges and paths in GBWT form)\n"
  "  -fullPath              : write X and Y tags with non-syncmer sequence ends to .1path, .1gbwt\n"
  "  -writeSeq              : write a .1seq file (paths converted back to sequences)\n"
  "possible inputs are:\n"
  "  <sequence file>        : any of fasta[.gz], fastq[.gz], BAM/CRAM/SAM, .1seq\n"
  "  <.1path file>          : sequences as lists of kmers, with optional prefix and suffix DNA\n"
  "  <.1gbwt file>          : graph BWT with paths\n"
  "Operations are carried out in order as they are parsed, with some setting up future actions,\n"
  "e.g. changing the outfile prefix affects following lower case options for file opening\n"
  "Some output files, e.g. .1gbwt will be output at the end, after all inputs are processed,\n"
  "whereas others, e.g. .1path are written as inputs are processed.\n" ;

int main (int argc, char *argv[])
{
  char       *outPrefix = "syngOut" ;
  char       *syncFileIn = 0, *syncFaOut = 0 ;
  int         nThread = 8 ;
  pthread_t  *threads ;
  ThreadInfo *threadInfo ;
  SeqIO      *sio = 0 ;
  KmerHash   *kh = 0, *khNew = 0 ;
  OneSchema  *schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile    *ofSync = 0, *ofGBWT = 0, *ofPath = 0, *ofSeq = 0 ;
  SyngBWT    *gbwtOut = 0 ;
  Params      params ;
  bool        isAddSyncmers = true, isHistK = false, isHistNewK = false ;
  I64         i, j ; // general purpose indices
  
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
      { syncFileIn = argv[1] ; argc -= 2 ; argv +=2 ;
	OneFile *of = oneFileOpenRead (syncFileIn, schema, "khash", 1) ;
	if (!of) die ("failed to open syncmer file %s to read", syncFileIn) ;

	readParams (of, &params) ;             // read the syncmer hash parameters
	kh = kmerHashReadOneFile (of) ;        // read in the kmerhash
	if (kh->len != params.w + params.k)
	  die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;
	for (i = 1 ; i <= kmerHashMax(kh) ; ++i) totSync += arr(kh->count,i,I64) ;
	fprintf (stderr, "read %llu syncmers from %s with total count %llu\n",
		kmerHashMax(kh), syncFileIn, totSync) ;
	oneFileClose (of) ; // close the old file
	timeUpdate (stderr) ;
      }
    else if (!strcmp (*argv, "-zeroK"))
      { if (!kh) die ("-zeroK error: no khash file read in yet") ;
	memset (arrp(kh->count,0,I64), 0, arrayMax(kh->count)*sizeof(I64)) ;
	argc-- ; argv++ ; 
      }
    else if (!strcmp (*argv, "-histK")) { isHistK = true ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-noNewK"))   { isAddSyncmers = false ; argc-- ; argv++ ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-writeK"))
      { char *fname = fnameTag (outPrefix,"1khash") ;
	if (syncFileIn)
	  { OneFile *ofIn = oneFileOpenRead (syncFileIn, 0, "khash", 1) ;
	    ofSync = oneFileOpenWriteFrom (fname, ofIn, true, 1) ;
	    oneFileClose (ofIn) ;
	  }
	else
	  ofSync = oneFileOpenWriteNew (fname, schema, "khash", true, 1) ;
	if (!ofSync) die ("failed to open syncmer file %s to write", fname) ;
	oneAddProvenance (ofSync, "syng", SYNG_VERSION, getCommandLine()) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writeKfa"))
      { syncFaOut = outPrefix ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writeGBWT"))
      { char *fname = fnameTag (outPrefix,"1gbwt") ;
	ofGBWT = oneFileOpenWriteNew (fname, schema, "gbwt", true, 1) ;
	if (!ofGBWT) die ("failed to open gbwt file %s to write", fname) ;
	oneAddProvenance (ofGBWT, "syng", SYNG_VERSION, getCommandLine()) ;
	gbwtOut = syngBWTcreate (params.k + params.w) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writePath"))
      { char *fname = fnameTag (outPrefix,"1path") ;
	ofPath = oneFileOpenWriteNew (fname, schema, "path", true, 1) ;
	if (!ofPath) die ("failed to open path file %s to write", fname) ;
	oneAddProvenance (ofPath, "syng", SYNG_VERSION, getCommandLine()) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-writeSeq"))
      { char *fname = fnameTag (outPrefix,"1seq") ;
	ofSeq = oneFileOpenWriteNew (fname, schema, "seq", true, 1) ;
	if (!ofSeq) die ("failed to open seq file %s to write", fname) ;
	oneAddProvenance (ofSeq, "syng", SYNG_VERSION, getCommandLine()) ;
	argc -= 1 ; argv += 1 ;
      }
    else if (!strcmp (*argv, "-histK")) { isHistK = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-histNewK")) { isHistNewK = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-checkNewK"))
      { khNew = kmerHashCreate (0, params.w + params.k) ;
	--argc ; ++argv ;
      }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  fprintf (stderr, "k, w, seed are %d %d %d\n", params.k, params.w, params.seed) ;
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly

  if (!kh) kh = kmerHashCreate (0, params.w + params.k) ; 

  int maxSources ; // collect sources from remaining files
  char **sources = collectSources (argc, argv, &maxSources) ;
  if (sources)
    { if (ofSync) addSourceReferences (ofSync, sources) ;
      if (ofGBWT) addSourceReferences (ofGBWT, sources) ;
      if (ofPath) addSourceReferences (ofPath, sources) ;
      if (ofSeq)  addSourceReferences (ofSeq, sources) ;
      char **s = sources ; while (*s) { free (*s) ; ++s ; }
      newFree (sources, maxSources, char*) ;
    }
  
  U64 nSeq0 = 0, totSeq0 = 0, totSync0 = totSync, syncMax0 = kmerHashMax(kh), nNew = 0 ;

  Array newKmers = 0 ;
  if (!isAddSyncmers)
    { KMERSIZE = ((params.w + params.k) + 31) >> 5 ;  // number of U64 needed to pack a kmer
      newKmers = uArrayCreate (1 << 20, KMERSIZE*8) ; // use uArrayCreate because not a simple type
    }
  Array seqDir = 0 ;
  if (ofPath) seqDir = arrayCreate(1024,char) ; // might be needed if we will flip
  
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
      if (!(sio = seqIOopenRead (*argv, dna2indexConv, false)))
	die ("failed to open read sequence file %s", *argv) ;
      fprintf (stderr, "sequence file %d %s type %s ", nSource, *argv, seqIOtypeName[sio->type]) ;
      fflush (stderr) ;
      bool isDone = false ;
      while (!isDone) 
	{ for (i = 0 ; i < nThread ; ++i)  // read 100Mb DNA per thread and find nodes in parallel
	    { ThreadInfo *ti = &threadInfo[i] ;
	      arrayMax(ti->seq) = 0 ;
	      arrayMax(ti->seqLen) = 0 ;
	      int seqStart = 0 ;
	      while (arrayMax(ti->seq) < 100<<20 && seqIOread (sio))
		{ array(ti->seqLen, arrayMax(ti->seqLen), I64) = sio->seqLen ;
		  array(ti->seq, seqStart+sio->seqLen, char) = 0 ;
		  memcpy (arrp(ti->seq, seqStart, char), sqioSeq(sio), sio->seqLen) ;
		  seqStart += sio->seqLen ;
		  totSeq += sio->seqLen ; 
		}
	      if (!arrayMax(ti->seq)) // we are done
		{ // nThread = i ;
		  isDone = true ;
		}
	    }
	  for (i = 0 ; i < nThread ; ++i) // create threads
	    pthread_create (&threads[i], 0, threadProcess, &threadInfo[i]) ;
	  for (i = 0 ; i < nThread ; ++i)
	    pthread_join (threads[i], 0) ; // wait for threads to complete
	  for (i = 0 ; i < nThread ; ++i) //
	    { ThreadInfo *ti = threadInfo + i ;
	      I64 *posList = arrp(ti->pos, 0, I64) ;
	      I64 *syncList = arrp(ti->sync, 0, I64) ;
	      char *seq = arrp(ti->seq, 0, char) ;
	      int iRead ;
	      for (iRead = 0 ; iRead < arrayMax(ti->posLen) ; ++iRead)
		{ I64 posLen = arr(ti->posLen, iRead, I64) ;
		  I64 *sync = syncList ;
		  for (j = 0 ; j < posLen ; ++j, ++sync)
		    if (*sync) // need to increment kh->count, because FindThreadSafe did not do that
		      ++array(kh->count,((*sync >= 0) ? *sync : -*sync),I64) ;
		    else if (isAddSyncmers)
		      kmerHashAdd (kh, seq + posList[j], sync) ;
		    else // write the canonicalised kmer into newKmers array
		      { kmerHashFindThreadSafe (kh, seq + posList[j], 0,
						arrayp(newKmers,arrayMax(newKmers),U64)) ;
			if (khNew && !(++nNew % 100)) // print out high coverage new kmers
			  { I64 syncNew ;
			    kmerHashAdd (khNew, seq + posList[j], &syncNew) ;
			    if (syncNew < 0) syncNew = -syncNew ;
			    if (arr(khNew->count,syncNew,I64) >= 10)
			      printf ("read %8lld syncNew %6lld %s pos %2lld in %s\n",
				      nSeq, syncNew, kmerHashSeq(khNew,syncNew), posList[j], seq) ;
			  }
		      }
		  if (ofGBWT && posLen)
		    { oneInt(ofGBWT, 0) = nSource ;
		      oneInt(ofGBWT, 1) = nSeq-nSeq0+1 ;
		      oneWriteLine (ofGBWT, 'P', 0, 0) ;
		      if (posLen > 1)
			{ SyngBWTpath *sbp = syngBWTpathStart (gbwtOut, syncList[0]) ;
			  oneInt(ofGBWT, 0) = syncList[0] ;
			  oneInt(ofGBWT, 1) = sbp->startCount ;
			  oneInt(ofGBWT, 2) = posLen ;
			  oneWriteLine (ofGBWT, 'Z', 0, 0) ;
			  for (j = 1 ; j < posLen ; ++j)
			    syngBWTpathAdd (sbp, syncList[j], posList[j] - posList[j-1]) ;
			  syngBWTpathFinish (sbp) ;
			  // now add the reverse path
			  sbp = syngBWTpathStart (gbwtOut, -syncList[posLen-1]) ;
			  for (j = posLen-2 ; j >= 0 ; --j)
			    syngBWTpathAdd (sbp, -syncList[j], posList[j+1] - posList[j]) ;
			  syngBWTpathFinish (sbp) ;
			}
		    }
		  if (ofPath && posLen) // do this last because it can corrupt syncList[]
		    { // test if worth flipping sites
		      bool isFlip = false ;
		      I64 *sync = syncList, sync0 = *sync++ ;
		      for (j = 1 ; j < posLen ; ++j, ++sync)
			if (*sync > sync0-64 && *sync < sync0+64) sync0 = *sync ;
			else if (-*sync > sync0-64 && -*sync < sync0+64)
			  { isFlip = true ; sync0 = -*sync ; }
			else
			  { isFlip = false ; break ; }
		      if (isFlip)
			{ arrayMax(seqDir) = 0 ;
			  array(seqDir,posLen,char) = 0 ;
			  for (j = 0 ; j < posLen ; ++j)
			    if (syncList[j] >= 0)
			      arr(seqDir,j,char) = '+' ;
			    else
			      { arr(seqDir,j,char) = '-' ;
				syncList[j] = -syncList[j] ;
			      }
			}
		      oneInt(ofPath, 0) = nSource ;
		      oneInt(ofPath, 1) = nSeq-nSeq0+1 ;
		      oneWriteLine (ofPath, 'P', 0, 0) ;
		      oneWriteLine (ofPath, 'z', posLen, syncList) ;
		      oneWriteLine (ofPath, 'o', posLen, posList) ;
		      if (isFlip) oneWriteLine (ofPath, 'd', posLen, arrp(seqDir,0,char)) ;
		    }
		  nSeq++ ; totSync += posLen ;
		  seq += arr(ti->seqLen, iRead, I64) ;
		  posList += posLen ;
		  syncList += posLen ;
		} // read
	    } // thread
	} // isDone: end of file
      seqIOclose (sio) ;
      fprintf (stderr, "had %llu sequences %llu bp, yielding %llu syncs with %llu extra syncmers\n",
	      nSeq-nSeq0, totSeq-totSeq0, totSync-totSync0, kmerHashMax(kh)-syncMax0) ;
      timeUpdate (stderr) ;
      nSeq0 = nSeq ; totSeq0 = totSeq ; totSync0 = totSync ; syncMax0 = kmerHashMax(kh) ;
      ++argv ;
    } // source file

  fprintf (stderr, "Total for this run %llu sequences, total length %llu\n", nSeq, totSeq) ;
  fprintf (stderr, "Overall total %llu instances of %llu syncmers, average %.2f coverage\n", 
	   totSync, kmerHashMax(kh), totSync / (double)kmerHashMax(kh)) ; 
  
  if (ofSeq) oneFileClose (ofSeq) ;
  if (ofPath) oneFileClose (ofPath) ;
  if (seqDir) arrayDestroy (seqDir) ;

  if (newKmers)
    { if (isHistNewK) // histogram of new kmers
	{ arraySort (newKmers, kmerOrder) ;
	  Array hist = arrayCreate (1024, I64) ;
	  I64 last = 0 ;
	  for (i = 1 ; i < arrayMax(newKmers) ; ++i)
	    if (kmerOrder (arrp(newKmers,i,I64), arrp(newKmers,i-1,I64)))
	      { ++array(hist,i-last,I64) ; last = i ; }
	  ++array(hist,i-last,I64) ;
	  qhist(hist,stdout) ;
	  arrayDestroy (hist) ;
	  timeUpdate (stderr) ;
	}
      arrayDestroy (newKmers) ;
    }

  if (isHistK) // histogram of kmers
    { Array hist = arrayCreate (1024, I64) ;
      for (i = 1 ; i < arrayMax(kh->count) ; ++i) ++array(hist,arr(kh->count,i,I64),I64) ;
      qhist(hist,stdout) ;
      arrayDestroy (hist) ;
      timeUpdate (stderr) ;
    }

  // now write all the files requested

  if (ofSync) // write the kmerhash as ONEcode
    { writeParams (ofSync, &params) ;
      kmerHashWriteOneFile (kh, ofSync) ;
      fprintf (stderr, "wrote %llu syncmers to file %s\n", kmerHashMax(kh), oneFileName(ofSync)) ;
      oneFileClose (ofSync) ;
      timeUpdate (stderr) ;
    }
  if (syncFaOut)
    { char  *fname = fnameTag (syncFaOut,"1kmer.fa.gz") ;
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
  if (ofGBWT)
    { syngBWTwrite (ofGBWT, gbwtOut) ;
      fprintf (stderr, "wrote gbwt to file %s\n", oneFileName(ofGBWT)) ;
      oneFileClose (ofGBWT) ;
      syngBWTdestroy (gbwtOut) ;
      timeUpdate (stderr) ;
    }

  if (kh) kmerHashDestroy (kh) ;
  if (khNew) kmerHashDestroy (khNew) ;
  seqhashDestroy (sh) ;
	  
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
