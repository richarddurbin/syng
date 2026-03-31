/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 24 17:29 2026 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include <pthread.h>

#include "seqio.h"
#include "seqhash.h"
#include "syng.h"

static OneSchema *schema ;

//#define PATH_CHECK
extern int pathCount ;

/**************** main program ********************/

typedef struct {
  I32	    sync ;	// initially 0 if not in kh; -ve if if reverse orientation
  U32       pos ;       // offset in the current sequence
} SyncPos ;

static char usage[] =
  "Usage: syngpath2gbwt XX.1path YY.1gbwt\n" ;

int main (int argc, char *argv[])
{
  I64 i, j, k ; // general purpose indices
  
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  if (argc != 3) { fprintf (stderr, "%s",usage) ; exit (0) ; }

  schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile *ofOut = oneFileOpenWriteNew (argv[2], schema, "gbwt", true, 1) ;
  if (!ofOut) die ("failed to open gbwt file %s to write", argv[2]) ;
  if (ofOut) oneAddProvenance (ofOut, "syngpath2gbwt", SYNG_VERSION, getCommandLine()) ;
  bool isOutputEnds = true ;

  OneFile *ofIn = oneFileOpenRead (argv[1], schema, "path", 1) ;
  if (!ofIn) die ("failed to open path file %s to read", argv[1]) ;
  oneAddReference (ofOut, argv[1], 1) ;

  // process the syncmer hash info line, if there is one
  int syncLen = 63 ; // default
  if (oneReadLine (ofIn) && ofIn->lineType == 'h') 
    { I64 kSync = oneInt(ofIn,0), wSync = oneInt(ofIn,1), seed = oneInt(ofIn,2) ;
      oneInt(ofOut,0) = kSync ; oneInt(ofOut,1) = wSync ; oneInt(ofOut,2) = seed ;
      oneWriteLine (ofOut, 'h', 0, 0) ;
      syncLen = kSync + wSync ;
      printf ("read h line k %lld w %lld len %d\n", kSync, wSync, syncLen) ;
    }

  // now process the paths
  pathCount = 0 ; // this is a debugging global in syngbwt3.c
  I64   totSync = 0 ;
  Array syncPos = arrayCreate (1024, SyncPos) ;
  SyngBWT *gbwt = syngBWTcreate (syncLen, 0) ;
  char  *dnaX, *dnaY ;
  I64   dnaXlen, dnaYlen ;
  I64   lastSource = 1 ;
  
  if (!oneGoto (ofIn, 'P', 1)) die ("failed to locate to a path (P) record in %s", argv[1]) ;
  oneReadLine (ofIn) ;
  while (ofIn->lineType == 'P') // a path line
    { ++pathCount ;
      I64 len = oneInt(ofIn,0) ;
      I64 source = oneInt(ofIn,1) ;
      if (source != lastSource)
	{ printf ("finished source %lld pathCount = %d : ", lastSource, pathCount) ;
	  lastSource = source ;
	  timeUpdate (stdout) ;
	}
      I64 inSource = oneInt(ofIn,2) ;
      I64 nSync = 0 ;
      SyncPos *sp = 0 ;
      dnaXlen = dnaYlen = 0 ;
      while (oneReadLine (ofIn) && ofIn->lineType != 'P')
	switch (ofIn->lineType)
	  { // expect a (z, o) line pair - X and Y are optional
	  case 'z':
	    nSync = oneLen(ofIn) ; totSync += nSync ;
	    arrayp(syncPos, nSync, SyncPos)->pos = 0 ;
	    sp = arrp(syncPos, 0, SyncPos) ;
	    I64 *iList = oneIntList(ofIn) ;
	    for (j = 0 ; j < nSync ; ++j) sp[j].sync = iList[j] ;
	    break ;
	  case 'o': // comes with z
	    if (oneLen(ofIn) != nSync)
	      die ("file error path %d o length %d != z length %d", pathCount, oneLen(ofIn), nSync) ;
	    iList = oneIntList(ofIn) ;
	    for (j = 0 ; j < nSync ; ++j) sp[j].pos = iList[j] ;
	    break ;
	  case 'X':
	    if (nSync && oneLen(ofIn) != sp->pos) die ("X error path %d", pathCount) ;
	    dnaXlen = oneLen(ofIn) ; dnaX = oneDNAchar (ofIn) ;
	    break ;
	  case 'Y':
	    if (nSync && oneLen(ofIn) != len - (sp[nSync-1].pos + syncLen))
	      die ("Y error path %d nSync %d oneLen %d si->len %d pos %d len %d",
		   pathCount, nSync, oneLen(ofIn), len, sp[nSync-1].pos, syncLen) ;
	    dnaYlen = oneLen(ofIn) ; dnaY = oneDNAchar (ofIn) ;
	    break ;
	  }

      if (nSync) // we got something
	{ oneInt(ofOut,0) = len ; // first write out the P line
	  oneInt(ofOut,1) = source ;
	  oneInt(ofOut,2) = inSource ;
	  oneWriteLine (ofOut, 'P', 0, 0) ;

	  SyngBWTpath *sbp = syngBWTpathStartNew (gbwt, sp->sync) ;
	  U32 j0 = sbp->jLast ; // used in PATH_CHECK
	  oneInt(ofOut, 0) = sp->sync ;
	  oneInt(ofOut, 1) = sp->pos ;
	  oneInt(ofOut, 2) = j0 ;
	  oneInt(ofOut, 3) = nSync ;
	  oneWriteLine (ofOut, 'Z', 0, 0) ;
	  if (dnaXlen) oneWriteLine (ofOut, 'X', dnaXlen, dnaX) ;
	  if (dnaYlen) oneWriteLine (ofOut, 'Y', dnaYlen, dnaY) ;
	  
	  for (k = 1 ; k < nSync ; ++k)
	    syngBWTpathAdd (sbp, sp[k].sync, sp[k].pos - sp[k-1].pos) ;
	  syngBWTpathFinish (sbp) ;
#ifdef PATH_CHECK			  
	  sbp = syngBWTpathStartOld (gbwt, sp->sync, j0) ;
	  I32 sync ;
	  U32 pos ;
	  for (k = 1 ; k < nSync ; ++k)
	    if (syngBWTpathNext (sbp, &sync, &pos))
	      { if (sync != sp[k].sync)
		  die ("path %d sync %d mismatch %d != %d", pathCount, k, sync, sp[k].sync) ;
		if (pos != sp[k].pos - sp[k-1].pos)
		  die ("path %d offset %d mismatch %d != %d - %d sync %d nSync %d",
		       pathCount, k, pos, sp[k].pos, sp[k-1].pos, sync, nSync) ;
	      }
	    else
	      die ("failed pathNext: path %d k %d", pathCount, k) ;
	  syngBWTpathDestroy (sbp) ;
#endif
	  // now add the reverse path
	  sbp = syngBWTpathStartNew (gbwt, -sp[nSync-1].sync) ;
	  for (k = nSync-2 ; k >= 0 ; --k)
	    syngBWTpathAdd (sbp, -sp[k].sync, sp[k+1].pos - sp[k].pos) ;
	  syngBWTpathFinish (sbp) ;
	  // printf (" path %d processed\n", pathCount) ;
	}
    } // end of file loop
  oneFileClose (ofIn) ;
  arrayDestroy (syncPos) ;
  printf ("path file %s had %'d paths containing %'llu syncmers\n", argv[1], pathCount, totSync) ;
  timeUpdate (stdout) ;

  syngBWTwrite (ofOut, gbwt) ;
  oneFileClose (ofOut) ;
  syngBWTdestroy (gbwt) ;
  printf ("wrote gbwt to file %s\n", oneFileName(ofOut)) ;
  timeUpdate (stdout) ;

  printf ("total: ") ; timeTotal (stdout) ;
  return 0 ;
}

/*********************** end of file **********************/
