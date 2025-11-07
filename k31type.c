/*  File: k32type.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov  7 22:14 2025 (rd109)
 * Created: Fri Nov  7 16:00:02 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "hash.h"
#include "seqio.h"
#include <ctype.h>

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;
  if (argc <= 1) die ("usage: syngtype <31merFile> <seqFile>+\n  31merFile is one 31mer per line") ;

  FILE *f = fopen (*argv, "r") ;
  if (!f) die ("failed to open kmerFile %s", *argv) ;
  Hash hash = hashCreate (8192) ; // of 2-bit packed 31-mers (62 bits) - this large for efficiency
  U64 forward = 0, reverse = 0 ;
  int nBases = 0, nKmers = 1 ;
  printf ("#file\t%s", fgetword(f)) ;
  while (!feof (f))
    { int c = fgetc(f) ;
      if (feof(f)) break ;
      I64 u = dna2indexConv[c] ;
      if (u < 0 || u >3)
	die ("non acgtACGT character %c %d at pos %d line %d of kmer file", c, c, nBases, nKmers) ;
      forward = (forward << 2) | u ;
      reverse = (reverse >> 2) | ((3-u) << 60) ;
      if (++nBases == 31)
	{ if (!feof(f) && fgetc(f) != '\n')
	    die ("no end of line at position 32 line %d of kmer file", nKmers) ;
	  if (!hashAdd (hash, hashInt((I64)forward), 0)) die ("duplicate 32mer line %d", nKmers) ;
	  if (!hashAdd (hash, hashInt((I64)reverse), 0)) die ("IMPOSSIBLE") ;
	  nBases = 0 ; forward = 0 ; reverse = 0 ;
	  printf ("\t%s", fgetword(f)) ;
	  ++nKmers ;
	}
    }
  putchar ('\n') ;
  if (nBases) die ("incomplete 31-mer length %d at end of kmer file", nBases) ;
  --nKmers ; // we were one ahead
  fclose (f) ;
  --argc ; ++argv ;

  U64 mask = ((U64)1 << 62) - 1 ;
  while (argc--)
    { SeqIO *sio = seqIOopenRead (*argv, dna2index4Conv, false) ;
      if (!sio) die ("failed to open %s to read", *argv) ;

      U64 N = 0, *counts = new0 (nKmers, U64) ;
      int k ;
      while (seqIOread (sio))
	{ if (sio->seqLen < 31) continue ;
	  forward = 0 ;
	  char *s = sqioSeq(sio) ;
	  U64 n ;
	  for (n = 0 ; n < 30 ; ++n) forward = (forward << 2) | *s++ ;
	  while (n++ < sio->seqLen)
	    { forward = ((forward << 2) | *s++) & mask ;
	      if (hashFind (hash, hashInt((I64)forward), &k)) ++counts[k/2] ;
	      ++N ;
	    }
	}
      seqIOclose (sio) ;

      printf ("%s\t%llu", *argv, (long long)N) ;
      int i ;
      for (i = 0 ; i < nKmers ; ++i) printf ("\t%lld", (long long)counts[i]) ;
      putchar ('\n') ;

      newFree (counts, nKmers, U64) ;
      ++argv ;
    }

  hashDestroy (hash) ;
}
