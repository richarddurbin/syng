/*  File: syngmap.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  2 19:03 2024 (rd109)
 * Created: Mon Dec  2 16:18:27 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include "seqhash.h"
#include "ONElib.h"

#include "syng.h"

typedef struct {
  I32 seq, pos ;
} Loc ;

int main (int argc, char *argv[])
{
  storeCommandLine (argc--, argv++) ;
  timeUpdate (0) ;
  
  if (argc != 2) die ("Usage: syngmap <genome> <reads>") ;

  SeqIO    *gio = seqIOopenRead (argv[0], dna2indexConv, false) ;
  if (!gio) die ("failed to open genome sequence file %s", argv[0]) ;
  SeqIO    *rio = seqIOopenRead (argv[1], dna2indexConv, false) ;
  if (!rio) die ("failed to open read sequence file %s", argv[1]) ;

  // first build the kmerhash index
  
  Seqhash  *sh = seqhashCreate (8, 56, 7) ;
  KmerHash *kh = kmerHashCreate (0, 63) ;
  Array    loc = arrayCreate (1 << 24, Loc) ;
  I32      nSeq = 0 ;
  I64      totSeq = 0, totSync = 0, uniqueSync = 0, sync ;
  while (seqIOread (gio))
    { ++nSeq ; totSeq += gio->seqLen ;
      SeqhashIterator *sit = syncmerIterator (sh, sqioSeq(gio), gio->seqLen) ;
      int pos ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ ++totSync ;
	  if (kmerHashAdd (kh, sqioSeq(gio)+pos, &sync))
	    { Loc *l = arrayp(loc, (sync > 0) ? sync : -sync, Loc) ;
	      l->seq = nSeq ; l->pos = pos ;
	      ++uniqueSync ;
	    }
	  else
	    { Loc *l = arrp(loc, (sync > 0) ? sync : -sync, Loc) ;
	      if (l->seq) { --uniqueSync ; l->seq = 0 ; }
	    }
	}
    }
  seqIOclose (gio) ;
  fprintf (stderr, "%s read %d sequences %lld bp %lld syncmers %lld distinct %lld unique\n",
	   argv[0], nSeq, totSeq, totSync, arrayMax(loc)-1, uniqueSync) ;
  timeUpdate (stderr) ;

  nSeq = 0 ;
  totSeq = 0, totSync = 0, uniqueSync = 0 ;
  I64 foundSeq1 ;
  while (seqIOread (rio))
    { I32 pos, s1 = 0, p1 = 0, s2 = 0, p2 = 0 ;
      ++nSeq ; totSeq += rio->seqLen ;
      SeqhashIterator *sit = syncmerIterator (sh, sqioSeq(rio), rio->seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ ++totSync ;
	  if (kmerHashFind (kh, sqioSeq(rio)+pos, &sync))
	    { Loc *l = arrp(loc, (sync > 0) ? sync : -sync, Loc) ;
	      if (!s1) { s1 = l->seq ; p1 = l-> pos ; }
	      else if (s1 != l->seq) s1 = -1 ;
	      else if (p1 < l->pos - rio->seqLen || p1 > l->pos + rio->seqLen) s1 = -rio->seqLen ;
	    }
	}
      if (!seqIOread (rio)) break ;
      ++nSeq ; totSeq += rio->seqLen ;
      sit = syncmerIterator (sh, sqioSeq(rio), rio->seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ ++totSync ;
	  if (kmerHashFind (kh, sqioSeq(rio)+pos, &sync))
	    { Loc *l = arrp(loc, (sync > 0) ? sync : -sync, Loc) ;
	      if (!s2) { s2 = l->seq ; p2 = l->pos ; }
	      else if (s2 != l->seq) s2 = -1 ;
	      else if (p2 < l->pos - rio->seqLen || p2 > l->pos + rio->seqLen) s2 = -rio->seqLen ;
	    }
	}
      if (s1 > 0 && s2 > 0 && p1 >= 0 && p2 >= 0)
	printf ("%d %d %d %d\n", s1, s2, p1, p2) ;
    }
  seqIOclose (rio) ;
  fprintf (stderr, "%s read %d sequences %lld bp %lld syncmers\n",
	   argv[1], nSeq, totSeq, totSync) ;
  timeUpdate (stderr) ;
  timeTotal (stderr) ;
}
