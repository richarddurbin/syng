/*  File: syngmap.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: maps sequences to a <syng>.1gbwt with its <syng>.1khash
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 14 09:28 2025 (rd109)
 * Created: Mon Dec  2 16:18:27 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "syng.h"
#include "seqhash.h"

typedef struct {
  I32 seq, pos ;
} Loc ;

/****************************************************/

typedef struct {
  int w, k, seed ;
} Params ;

static void readParams (OneFile *of, Params *p)
{ while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType != 'h') die ("sync file %s has no 'h' parameters record", oneFileName(of)) ;
  p->k = oneInt(of,0) ; p->w = oneInt(of,1) ; p->seed = oneInt(of,2) ;
  fprintf (stderr, "read syncmer parameters k %d w %d (size %d) seed %d\n",
	   p->k, p->w, p->k + p->w, p->seed) ;
}

static void checkParams (OneFile *of, Params *p)
{
  if (oneInt(of,0) != p->k || oneInt(of,1) != p->w || oneInt(of,2) != p->seed)
    die ("hash parameters mismatch: (k,w,s) file (%d,%d,%d) != code (%d,%d,%d)",
	 oneInt(of,0), oneInt(of,1), oneInt(of,2), p->k, p->w, p->seed) ;
}

/****************************************************/

static int filterG = 0, filterQ = 0 ;

typedef struct {
  Seqhash  *sh ;
  KmerHash *kh ;
  SyngBWT  *gbwt ;
  OneFile  *of ;        // output file
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  Array     seqInfo ;   // of seqInfo - see below
  Array     pos ;
  Array     posLen ;
  Array     sync ;
} ThreadInfo ;

typedef struct {
  I64       seqLen ;
  char*     id ;        // if requested
  double    avQ ;
  I32       inFile ;    // position in filep
} SeqInfo ;

static void *threadProcessRead (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;
  I64 iSync ;
  U64 *uBuf = new(ti->kh->plen,U64) ;

  arrayMax(ti->pos) = 0 ;
  arrayMax(ti->posLen) = 0 ;
  
  for (i = 0 ; i < arrayMax(ti->seqInfo) ; ++i)
    { SeqInfo *si = arrp(ti->seqInfo, i, SeqInfo) ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += si->seqLen ;
      if (filterG)
	{ int j, nG = 0 ;
	  char *s = seq ;
	  for (j = 0 ; j < si->seqLen ; ++j)
	    if (*s++ != 2) nG = 0 ;
	    else if (++nG > filterG) { si->seqLen = 0 ; break ; }
	}
      if (filterQ && si->avQ < filterQ) si->seqLen = 0 ;
      int pos, posStart = arrayMax(ti->pos) ;
      SeqhashIterator *sit = syncmerIterator (ti->sh, seq, si->seqLen) ;
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

/****************************************************/

static char *usage = 
  "Usage: syngmap [options]* <.1khash file> <.1gbwt file> <sequence file>\n"
  "possible options are:\n"
  "  -T <threads>           : [8] number of threads\n"
  "  -o <outfile prefix>    : default is syngmap, file type is .1map\n"
  "  -writeNewK <file prefix> : write new syncmers as a .1khash file; implies -N\n"
  "  -outputNames           : write the names of path sequences I lines\n"
  "  -filterG <nG>          : remove input sequences with runs of >nG consecutive Gs (bad Illumina reads)\n"
  "  -filterQ <QT>          : remove input sequences with average quality score below QT\n"
  "  -filterIllumina        : equivalent to -filterPolyG 60 -filterQual 30\n" ;

int main (int argc, char *argv[])
{
  storeCommandLine (argc--, argv++) ;
  timeUpdate (0) ;

  if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }

  int   nThread = 8 ;
  char *outPrefix = "syngmap" ;

  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-T") && argc > 1) { nThread = atoi(argv[1]) ; argc -=2 ; argv +=2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-filterG") && argc > 1)
      { filterG = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-filterQ") && argc > 1)
      { filterQ = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-filterIllumina")) { filterG = 60 ; filterQ = 20 ; --argc ; ++argv ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  if (argc != 3) die ("missing the three required arguments\n%s", usage) ;

  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
  if (!schema) die ("problem creating schema - needs debugging/recompilation") ;

  // open all the files here, so we can die quickly if any file opens fail
  OneFile *ofK = oneFileOpenRead (argv[0], schema, "khash", 1) ;
  if (!ofK) die ("failed to open .1khash file %s", argv[0]) ;
  OneFile *ofGBWT = oneFileOpenRead (argv[1], schema, "gbwt", 1) ;
  if (!ofGBWT) die ("failed to open .1gbwt file %s", argv[1]) ;
  SeqIO   *sio = seqIOopenRead (argv[2], dna2index4Conv, false) ;
  if (!sio) die ("failed to open sequence file %s", argv[2]) ;

  // this will be our output file
  OneFile *ofOut = oneFileOpenWriteNew (fnameTag(outPrefix,"1map"), schema, "map", true, nThread) ;
  if (!ofOut) die ("failed to open output file %s", argv[2]) ;
  oneAddProvenance (ofOut, "syngmap", SYNG_VERSION, getCommandLine()) ;
  oneAddReference (ofOut, argv[0], 1) ;
  oneAddReference (ofOut, argv[1], 2) ;
  oneAddReference (ofOut, argv[2], 3) ;

  // read the khash 
  Params     params ;
  readParams (ofK, &params) ;                    // read the syncmer hash parameters
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly
  KmerHash  *kh = kmerHashReadOneFile (ofK) ;    // read in the kmerhash
  if (kh->len != params.w + params.k)
    die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;
  printf ("read %llu kmers from %s\n", kh->max, oneFileName (ofK)) ;
  oneFileClose (ofK) ;
  timeUpdate (stdout) ;

  // read the BWT
  SyngBWT *sgb = syngBWTread (ofGBWT) ;
  if (!sgb) die ("failed to read GBWT from %s", argv[1]) ;
  printf ("read GBWT with %llu nodes from %s\n", arrayMax(sgb->node)-1, argv[1]) ;
  oneFileClose (ofGBWT) ;
  timeUpdate (stdout) ;

  int i, j ; // general index variables
  
  // set up the threads
  if (nThread < 1) die ("number of threads %d must be at least 1", nThread) ;
  pthread_t  *threads = new (nThread, pthread_t) ;
  ThreadInfo *threadInfo = new0 (nThread, ThreadInfo) ;
  for (i = 0 ; i < nThread ; ++i)
    { threadInfo[i].sh     = sh ;
      threadInfo[i].kh     = kh ;
      threadInfo[i].gbwt   = sgb ;
      threadInfo[i].of     = ofOut + i ; // i'th thread
      threadInfo[i].seq    = arrayCreate (101<<20, char) ; // 101 Mb
      threadInfo[i].seqInfo = arrayCreate (20000, SeqInfo) ;
    }

  // process the sequences
  I64 totSeq = 0 ;
  bool isDone = false ;
  while (seqIOread (sio))
    { for (i = 0 ; i < nThread ; ++i)  // read 100Mb DNA per thread and find kmers in parallel
	{ ThreadInfo *ti = &threadInfo[i] ;
	  arrayMax(ti->seq) = 0 ;
	  arrayMax(ti->seqInfo) = 0 ;
	  int seqStart = 0 ;
	  while (arrayMax(ti->seq) < 100<<20 && seqIOread (sio))
	    { int kSeq = arrayMax(ti->seqInfo) ;
	      I64 seqLen = sio->seqLen ;
	      SeqInfo *si = arrayp(ti->seqInfo, kSeq, SeqInfo) ;
	      si->seqLen = seqLen ;
	      array(ti->seq, seqStart+seqLen, char) = 0 ;
	      memcpy (arrp(ti->seq, seqStart, char), sqioSeq(sio), seqLen) ;
	      seqStart += seqLen ;
	      totSeq += seqLen ;
	      if (sio->isQual && seqLen)
		{ I64 totQ = 0 ;
		  char *q = sqioQual(sio) ;
		  for (i = 0 ; i < seqLen ; ++i) totQ += *q++ ;
		  si->avQ = totQ / (double)seqLen ;
		}
	    }
	  if (!arrayMax(ti->seq)) isDone = true ; // we are done after processing this lot
	} // loading thread
      for (i = 0 ; i < nThread ; ++i) // create threads
	pthread_create (&threads[i], 0, threadProcessRead, &threadInfo[i]) ;
      for (i = 0 ; i < nThread ; ++i)
	pthread_join (threads[i], 0) ; // wait for threads to complete
    }
  seqIOclose (sio) ;

  oneFileClose (ofOut) ;

  syngBWTdestroy (sgb) ;
  kmerHashDestroy (kh) ;
  oneSchemaDestroy (schema) ;				   
				   
  timeUpdate (stderr) ;
  timeTotal (stderr) ;
}
