/*  File: syngmap.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: maps sequences to a <syng>.1gbwt with its <syng>.1khash
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 14 00:49 2026 (rd109)
 * Created: Mon Dec  2 16:18:27 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "syng.h"
#include "seqhash.h"

extern int pathCount ;

static const int DEBUG = 0 ;

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

static int filterG = 0, filterQ = 0 ; // global across threads

typedef struct {
  Seqhash  *sh ;
  KmerHash *kh ;
  SyngBWT  *gbwt ;
  OneFile  *of ;        // output file
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  Array     seqInfo ;   // of seqInfo - see below
  U64       nMatch, nNoMatch, matchLen, matchCount ;
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
  U64 *uBuf = new(ti->kh->plen,U64) ;

  ti->nMatch = ti->matchLen = ti->matchCount = ti->nNoMatch = 0 ;
  
  for (i = 0 ; i < arrayMax(ti->seqInfo) ; ++i)
    { SeqInfo *si = arrp(ti->seqInfo, i, SeqInfo) ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += si->seqLen ;
      ++pathCount ;
      if (filterG)
	{ int j, nG = 0 ;
	  char *s = seq ;
	  for (j = 0 ; j < si->seqLen ; ++j)
	    if (*s++ != 2) nG = 0 ;
	    else if (++nG > filterG) { si->seqLen = 0 ; break ; }
	}
      if (filterQ && si->avQ < filterQ) si->seqLen = 0 ;
      if (si->seqLen)
	{ SeqhashIterator *sit = syncmerIterator (ti->sh, seq, si->seqLen) ;
	  SyngBWTpath *sbp = 0 ;
	  int pos, firstPos, lastPos ;
	  U32 low = 0, high ;
	  U64 nMatchIn = ti->nMatch ;
	  while (syncmerNext (sit, 0, &pos, 0))
	    { I64 iSync = 0 ;
	      kmerHashFindThreadSafe (ti->kh, seq+pos, &iSync, uBuf) ;
	      if (iSync && (iSync > 2 || iSync < -2)) // ignore poly-N
		{ I32 sync = iSync ;
		  if (!sbp)
		    { sbp = syngBWTmatchStart (ti->gbwt, sync, &high) ;
		      ++ti->nMatch ; firstPos = pos ;
		    }
		  else if (!syngBWTmatchNext (sbp, sync, pos-lastPos, &low, &high))
		    { syngBWTpathDestroy (sbp) ;
		      ti->matchCount += high - low ;
		      ti->matchLen += lastPos - firstPos + ti->kh->len ;
		      low = 0 ;
		      sbp = syngBWTmatchStart (ti->gbwt, sync, &high) ; // start a new match
		      ++ti->nMatch ;
		    }
		  lastPos = pos ;
		}
	      else if (sbp)
		{ syngBWTpathDestroy (sbp) ;
		  ti->matchCount += high - low ;
		  ti->matchLen += lastPos - firstPos + ti->kh->len ;
		  sbp = 0 ; low = 0 ;
		} // end of match
	    }
	  if (sbp)
	    { syngBWTpathDestroy (sbp) ; // clean up at end
	      ti->matchCount += high - low ;
	      ti->matchLen += lastPos - firstPos + ti->kh->len ;
	    }
	  seqhashIteratorDestroy (sit) ;
	  if (ti->nMatch == nMatchIn) ++ti->nNoMatch ;
	}
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
  "  -filterIllumina        : equivalent to -filterG 60 -filterQ 20\n" ;

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
  OneFile *ofGBWT = oneFileOpenRead (argv[1], schema, "gbwt", nThread) ;
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
  oneFileClose (ofGBWT) ;

  int i, j ; // general index variables

  if (DEBUG) nThread = 1 ;
  
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
  I64 nSeq = 0, totSeq = 0, totMatch = 0, totCount = 0, totLen = 0, totNoMatch = 0 ;
  bool isDone = false ;
  while (!isDone)
    { for (i = 0 ; !isDone && i < nThread ; ++i)  // read 100Mb DNA per thread and process in parallel
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
	      ++nSeq ; 
	      totSeq += seqLen ;
	      if (sio->isQual && seqLen)
		{ I64 totQ = 0 ;
		  char *q = sqioQual(sio) ;
		  for (i = 0 ; i < seqLen ; ++i) totQ += *q++ ;
		  si->avQ = totQ / (double)seqLen ;
		}
	    }
	  if (arrayMax(ti->seq))
	    pthread_create (&threads[i], 0, threadProcessRead, ti) ;
	  else
	    { isDone = true ; // we are done after processing this lot
	      nThread = i ;   // need to only wait for threads less than this one
	    }
	} // loading thread
      for (i = 0 ; i < nThread ; ++i)
	{ pthread_join (threads[i], 0) ; // wait for threads to complete
	  ThreadInfo *ti = &threadInfo[i] ;
	  totMatch += ti->nMatch ;
	  totNoMatch += ti->nNoMatch ;
	  totCount += ti->matchCount ;
	  totLen   += ti->matchLen ;
	}
      printf ("done %'lld seqs, totLen %'lld, matches %'lld, noMatch %'lld ",
	      nSeq, totSeq, totMatch, totNoMatch) ;
      if (totMatch) printf ("avLen %.1f avCnt %.1f ", totLen/(1.*totMatch), totCount/(1.*totMatch)) ;
      timeUpdate (stdout) ;
    }
  seqIOclose (sio) ;
  if (nSeq)
    { printf ("processed %'lld sequences total %'lld bp (av %.1f) yielding %'lld matches",
	      nSeq, totSeq, totSeq/(1.0*nSeq), totMatch) ;
      if (totMatch) printf (" avLen %.1f avCnt %.1f", totLen/(1.*totMatch), totCount/(1.*totMatch)) ;
      putchar ('\n') ;
      printf ("%'lld without matches\n", totNoMatch) ;
      timeUpdate (stdout) ;
    }
  
  oneFileClose (ofOut) ;

  syngBWTdestroy (sgb) ;
  kmerHashDestroy (kh) ;
  oneSchemaDestroy (schema) ;				   
				   
  printf ("Total usage: ") ; timeTotal (stdout) ;
}
