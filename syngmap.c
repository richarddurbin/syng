/*  File: syngmap.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: maps sequences to a <syng>.1gbwt with its <syng>.1khash
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 24 21:53 2026 (rd109)
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
  if (of->lineType != 'h') die ("syncmer file %s has no 'h' parameters record", oneFileName(of)) ;
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
  Seqhash    *sh ;
  SyncmerSet *sms ;
  SyngBWT    *gbwt ;
  OneFile    *ofIn ;      // input file
  OneFile    *ofOut ;     // output file
  Array       seq ;	// of char, input: concatenated sequences in index 0..3
  Array       seqInfo ;   // of seqInfo - see below
  U64         nMatch, nNoMatch, matchLen, matchCount ;
} ThreadInfo ;

typedef struct {
  I64         seqLen ;
  char*       id ;        // if requested
  double      avQ ;
  I32         index ;    // index in file - 1-based
} SeqInfo ;

static inline void writeMem (ThreadInfo *ti, int start, int end, int count)
{
  end += ti->sms->kh->len ;
  ++ti->nMatch ; ti->matchCount += count ; ti->matchLen += end - start ;
  oneInt(ti->ofOut, 0) = start ; oneInt(ti->ofOut, 1) = end ;
  oneInt(ti->ofOut, 2) = count ;
  oneWriteLine (ti->ofOut, 'M', 0, 0) ;
}

static inline void writeLoc (ThreadInfo *ti, I32 uniqueSync, int uniquePos, int firstPos)
{
  I64 loc = syncmerLoc (ti->sms, uniqueSync) ;
  I64 memLoc = loc > 0 ? loc - (uniquePos - firstPos) : loc + (uniquePos - firstPos) ;
  I64 file, path, offset ;
  if (syngBWTlocFind (ti->gbwt, memLoc, &file, &path, &offset))
    { oneInt(ti->ofOut, 0) = file ; oneInt(ti->ofOut, 1) = path ; oneInt(ti->ofOut, 2) = offset ;
      oneWriteLine (ti->ofOut, 'U', 0, 0) ;
    }
}

static void *threadProcessRead (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;
  U64 *uBuf = new(ti->sms->kh->plen,U64) ;

  ti->nMatch = ti->matchLen = ti->matchCount = ti->nNoMatch = 0 ;

  // printf ("THREAD nSeq %'lld pathCount %'d\n", arrayMax(ti->seqInfo), pathCount) ;

  Array syncStack = arrayCreate (4096, I32) ;
  Array posStack  = arrayCreate (4096, int) ;
  
  for (i = 0 ; i < arrayMax(ti->seqInfo) ; ++i)
    { SeqInfo *si = arrp(ti->seqInfo, i, SeqInfo) ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += si->seqLen ;
      ++pathCount ;
      
      oneInt(ti->ofOut, 0) = si->index ; oneInt(ti->ofOut, 1) = si->seqLen ;
      oneWriteLine (ti->ofOut, 'S', 0, 0) ;
      if (si->id)
	{ oneWriteLine (ti->ofOut, 'I', strlen(si->id), si->id) ;
	  free (si->id) ; // allocated by strdup()
	}
      char filter = 0 ;
      if (filterG)
	{ int j, nG = 0 ;
	  char *s = seq ;
	  for (j = 0 ; j < si->seqLen ; ++j)
	    if (*s++ != 2) nG = 0 ;
	    else if (++nG > filterG) { filter = 'G' ; break ; }
	}
      if (filterQ && si->avQ && si->avQ < filterQ) filter = 'Q' ;
      if (filter)
	{ oneChar (ti->ofOut, 0) = filter ;
	  oneWriteLine (ti->ofOut, 'F', 0, 0) ;
	}
      else
	{ SeqhashIterator *sit = syncmerIterator (ti->sh, seq, si->seqLen) ;
	  SyngBWTpath *sbp = 0 ;
	  int pos, firstPos, lastPos ;
	  U32 low, high ;
	  I32 uniqueSync = 0 ; int uniquePos ;
	  U64 nMatchIn = ti->nMatch ;
	  while (syncmerNext (sit, 0, &pos, 0))
	    { I64 iSync = 0 ;
	      kmerHashFindThreadSafe (ti->sms->kh, seq+pos, &iSync, uBuf) ;
	      // printf ("path %d pos %d sync %d\n", pathCount, pos, (int)iSync) ;
	      if (iSync && (iSync > 2 || iSync < -2)) // ignore poly-N
		{ I32 sync = iSync ;
		  array(syncStack,arrayMax(syncStack),I32) = sync ;
		  array(posStack, arrayMax(posStack), int) = pos ;
		  // printf ("  added to stack %d\n", (int)arrayMax(posStack)-1) ;
		  if (!sbp)
		    { low = 0 ; sbp = syngBWTmatchStart (ti->gbwt, sync, &high) ;
		      firstPos = pos ;
		      // printf ("  new sbp low %u high %u\n", low, high) ;
		    }
		  else
		    { // printf ("    matching sync %d off %d\n", sync, pos-lastPos) ;
		      if (!syngBWTmatchNext (sbp, sync, pos-lastPos, &low, &high))
			{ syngBWTpathDestroy (sbp) ;
			  // printf ("    no-match: writeMem start %d end %d count %d\n", firstPos, lastPos+ti->sms->kh->len, high-low) ;
			  writeMem (ti, firstPos, lastPos, high-low) ;
			  if (uniqueSync)
			    { writeLoc (ti, uniqueSync, uniquePos, firstPos) ; uniqueSync = 0 ; }
			  // run back to find the start of the next mem
			  low = 0 ; sbp = syngBWTmatchStart (ti->gbwt, -sync, &high) ;
			  // printf ("      reverse from %d low %u high %u\n", -sync, low, high) ;
			  int k = arrayMax(syncStack) - 2 ; // previous entry
			  while (true) // could be while (k), but want the check below
			    { // printf ("        matching k %d sync %d offset %d\n", k, -arr(syncStack,k,I32), arr(posStack,k+1,int)-arr(posStack,k,int)) ;
			      if (!syngBWTmatchNext (sbp, -arr(syncStack,k,I32),
						     arr(posStack,k+1,int)-arr(posStack,k,int),
						     &low, &high))
				break ;
			      // printf ("         match low %u high %u\n", low, high) ;
			      if (!k--) die ("problem reversing in memfinding") ;
			    }
			  ++k ; // increment to the last one that worked
			  firstPos = arr(posStack,k,int) ;
			  // printf ("      restart at pos %d\n", firstPos) ;
			  low = 0 ; sbp = syngBWTmatchStart (ti->gbwt, arr(syncStack,k,I32), &high) ;
			  while (++k < arrayMax(syncStack))
			    if (!syngBWTmatchNext (sbp, arr(syncStack,k,I32),
						   arr(posStack,k,int)-arr(posStack,k-1,int),
						   &low, &high))
			      die ("error reverting in memfinding") ;
			}
		      // else printf ("    match: low %u high %u\n", low, high) ;
		    }
		  if (syncmerCount (ti->sms, sync) == 1 && !uniqueSync)
		    { uniqueSync = sync ; uniquePos = pos ; }
		  lastPos = pos ;
		}
	      else
		{ if (sbp)
		    { syngBWTpathDestroy (sbp) ;
		      writeMem (ti, firstPos, lastPos, high-low) ;
		      if (uniqueSync)
		{ writeLoc (ti, uniqueSync, uniquePos, firstPos) ; uniqueSync = 0 ; }
		      sbp = 0 ; low = 0 ;
		      arrayMax(syncStack) = 0 ; // restart the syncStack
		      arrayMax(posStack)  = 0 ; // and the posStack
		    } // end of match
		  if (!iSync) // missing syncmer - write its sequence
		    { oneInt(ti->ofOut, 0) = pos ;
		      oneWriteLine (ti->ofOut, 'X', ti->sms->kh->len, seq+pos) ;
		    }
		}
	    }
	  if (sbp)
	    { syngBWTpathDestroy (sbp) ; // clean up at end
	      writeMem (ti, firstPos, lastPos, high-low) ;
	      if (uniqueSync) writeLoc (ti, uniqueSync, uniquePos, firstPos) ;
	    }
	  seqhashIteratorDestroy (sit) ;
	  if (ti->nMatch == nMatchIn) ++ti->nNoMatch ;
	}
    }

  newFree (uBuf, ti->sms->kh->plen, U64) ;
  arrayDestroy (syncStack) ;
  arrayDestroy (posStack) ;
  return 0 ;
}

/****************************************************/

static char *usage = 
  "Usage: syngmap [options]* <.1khash file> <.1gbwt file> <query sequence file>\n"
  "possible options are:\n"
  "  -T <threads>           : [8] number of threads\n"
  "  -o <outfile prefix>    : default is syngmap, file type is .1map\n"
  "  -writeNewK <file prefix> : write new syncmers as a .1khash file; implies -N\n"
  "  -outputIds             : write the identifiers of query sequences as I lines\n"
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
  bool  isOutputIds = false ;

  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-T") && argc > 1) { nThread = atoi(argv[1]) ; argc -=2 ; argv +=2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-outputIds")) { isOutputIds = true ; --argc ; ++argv ; }
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
  oneFileClose (ofK) ; // we will read it using syncmerSetRead() below
  OneFile *ofGBWT = oneFileOpenRead (argv[1], schema, "gbwt", nThread) ;
  if (!ofGBWT) die ("failed to open .1gbwt file %s", argv[1]) ;
  SeqIO   *sio = seqIOopenRead (argv[2], dna2index4Conv, true) ;
  if (!sio) die ("failed to open sequence file %s", argv[2]) ;

  // this will be our output file
  OneFile *ofOut = oneFileOpenWriteNew (fnameTag(outPrefix,"1map"), schema, "map", true, nThread) ;
  if (!ofOut) die ("failed to open output file %s", argv[2]) ;
  oneAddProvenance (ofOut, "syngmap", SYNG_VERSION, getCommandLine()) ;
  oneAddReference (ofOut, argv[0], 1) ;
  oneAddReference (ofOut, argv[1], 2) ;
  oneAddReference (ofOut, argv[2], 3) ;

  // read the SyncmerSet (KmerHash)
  SyncmerSet *sms = syncmerSetRead (argv[0]) ;
  Seqhash *sh = seqhashCreate (sms->params.k, sms->params.w+1, sms->params.seed) ; // need +1 for sync

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
      threadInfo[i].sms    = sms ;
      threadInfo[i].gbwt   = sgb ;
      threadInfo[i].ofOut  = ofOut + i ; // i'th thread
      threadInfo[i].seq    = arrayCreate (101<<20, char) ; // 101 Mb
      threadInfo[i].seqInfo = arrayCreate (20000, SeqInfo) ;
    }

  // process the sequences
  I64 nSeq = 0, totSeq = 0, totMatch = 0, totCount = 0, totLen = 0, totNoMatch = 0 ;
  I64 nLast = 0, totLast = 0 ;
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
	      totSeq += seqLen ;
	      si->index = ++nSeq ; 
	      if (isOutputIds) si->id = strdup(sqioId(sio)) ;
	      if (sio->isQual && seqLen)
		{ I64 totQ = 0 ;
		  char *q = sqioQual(sio) ;
		  for (int k = 0 ; k < seqLen ; ++k) totQ += *q++ ;
		  si->avQ = totQ / (double)seqLen ;
		  if (nSeq < 10) printf ("%d avQ %.1f\n", (int)nSeq, si->avQ) ;
		}
	    }
	  if (arrayMax(ti->seq))
	    pthread_create (&threads[i], 0, threadProcessRead, ti) ;
	  else
	    { isDone = true ; // we are done after processing this lot
	      nThread = i ;   // need to only wait for threads less than this one
	    }
	} // loading threads
      if (nThread)
	{ printf ("threads loaded with %'lld seqs %'lld bp: ", nSeq-nLast, totSeq-totLast) ;
	  nLast = nSeq ; totLast = totSeq ;
	  timeUpdate (stdout) ;
	  for (i = 0 ; i < nThread ; ++i)
	    { pthread_join (threads[i], 0) ; // wait for threads to complete
	      ThreadInfo *ti = &threadInfo[i] ;
	      totMatch += ti->nMatch ;
	      totNoMatch += ti->nNoMatch ;
	      totCount += ti->matchCount ;
	      totLen   += ti->matchLen ;
	    }
	  printf ("processed: ") ; timeUpdate (stdout) ;
	}
    }
  seqIOclose (sio) ;
  if (nSeq)
    { printf ("completed %'lld sequences total %'lld bp (av %'.1f) yielding %'lld matches",
	      nSeq, totSeq, totSeq/(1.0*nSeq), totMatch) ;
      if (totMatch) printf (" avLen %'.1f avCnt %.1f", totLen/(1.*totMatch), totCount/(1.*totMatch)) ;
      putchar ('\n') ;
      printf ("%'lld without matches\n", totNoMatch) ;
    }
  
  oneFileClose (ofOut) ;

  syngBWTdestroy    (sgb) ;
  syncmerSetDestroy (sms) ;
  oneSchemaDestroy  (schema) ;				   
				   
  printf ("Total usage: ") ; timeTotal (stdout) ;
}
