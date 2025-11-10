/*  File: k31type.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 10 12:22 2025 (rd109)
 * Created: Fri Nov  7 16:00:02 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "array.h"
#include "hash.h"
#include "seqio.h"
#include <ctype.h>

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  FILE *match = 0 ;
  if (argc > 1 && !strcmp (*argv, "-m"))
    { if (!(match = fopen (argv[1], "w"))) die ("failed to open match file %s to write", argv[1]) ;
      argc -= 2 ; argv += 2 ;
    }
  else if (argc && **argv == '-')
    die ("unknown option %s - run without arguments for usage", *argv) ;
  
  if (argc <= 1)
    { fprintf (stderr, "usage: k31type [-m <matchfile>] <31merFile> <seqFile>+\n"
	       "  <31merFile> is one 31mer per line preceeded by an id\n"
	       "  <matchfile> will contain match coordinates and sequences +/- 50bp\n") ;
      exit (1) ;
    }

  FILE *f = fopen (*argv, "r") ;
  if (!f) die ("failed to open kmerFile %s", *argv) ;
  Hash  hash = hashCreate (8192) ; // of 2-bit packed 31-mers (62 bits) - this large for efficiency
  Array ids  = arrayCreate (16, char*) ;
  U64 forward = 0, reverse = 0 ;
  int nBases = 0, nKmers = 1 ;
  char *id = fgetword(f) ;
  array(ids,arrayMax(ids),char*) = strdup(id) ;
  printf ("#file\t%s", id) ;
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
	  id = fgetword(f) ;
	  array(ids,arrayMax(ids),char*) = strdup(id) ;
	  printf ("\t%s", id) ;
	  ++nKmers ;
	}
    }
  putchar ('\n') ;
  if (nBases) die ("incomplete 31-mer length %d at end of kmer file", nBases) ;
  --nKmers ; // we were one ahead
  fclose (f) ;
  --argc ; ++argv ;

  U64  mask = ((U64)1 << 62) - 1 ;
  char buf[132] ; buf[131] = 0 ; // 31bp plus 50bp each side plus terminating 0
  while (argc--)
    { SeqIO *sio = seqIOopenRead (*argv, dna2index4Conv, false) ;
      if (!sio) die ("failed to open %s to read", *argv) ;

      U64 N = 0, *counts = new0 (nKmers, U64) ;
      int k ;
      while (seqIOread (sio))
	{ if (sio->seqLen < 31) continue ;
	  forward = 0 ;
	  char *seq = sqioSeq(sio), *s = seq ;
	  U64 n ;
	  for (n = 0 ; n < 30 ; ++n) forward = (forward << 2) | *s++ ;
	  while (n++ < sio->seqLen)
	    { forward = ((forward << 2) | *s++) & mask ;
	      if (hashFind (hash, hashInt((I64)forward), &k))
		{ ++counts[k/2] ;
		  if (match)
		    { fprintf (match, "%s\t%s\t%s", arr(ids,k/2,char*), *argv, sqioId(sio)) ;
		      if (k & 0x01) // reversed
			{ int pos = s - seq - 1 ;
			  fprintf (match, "\t-%d\t%d", k, pos) ;
			  int i = 0, j ;
			  for (j = pos+50 ; j > pos ; --j)
			    if (j >= sio->seqLen) buf[i++] = '.' ;
			    else buf[i++] = "tgca"[(int)seq[j]] ;
			  while (j > pos - 31) buf[i++] = "TGCA"[(int)seq[j--]] ;
			  while (j > pos - 81)
			    if (j < 0) { buf[i++] = '.' ; --j ; }
			    else buf[i++] = "tgca"[(int)seq[j--]] ;
			}
		      else
			{ int pos = s - seq - 31 ;
			  fprintf (match, "\t+%d\t%d", k, pos) ;
			  int i = 0, j ;
			  for (j = pos-50 ; j < pos ; ++j)
			    if (j <= 0) buf[i++] = '.' ;
			    else buf[i++] = "acgt"[(int)seq[j]] ;
			  while (j < pos + 31) buf[i++] = "ACGT"[(int)seq[j++]] ;
			  while (j < pos + 81)
			    if (j >= sio->seqLen) { buf[i++] = '.' ; ++j ; }
			    else buf[i++] = "acgt"[(int)seq[j++]] ;
			}
		      fprintf (match, "\t%s\n", buf) ;
		    }
		}
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

  if (match) fclose (match) ;
  int i ; for (i = 0 ; i < arrayMax(ids) ; ++i) free(arr(ids,i,char*)) ;
  arrayDestroy (ids) ;
  hashDestroy (hash) ;
}
