/*  File: seqhash.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: seqhash package - uses random bit patterns for bases and rotate/XOR
 	compile with -DTEST to test, and -DDEBUG to debug, -DNDEBUG turns off asserts
	see test main() at end for standard usage pattern
 * Exported functions: see seqhash.h
 * HISTORY:
 * Last edited: Sep 28 01:20 2024 (rd109)
 * Created: Sat Feb 24 19:20:18 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "seqhash.h"

Seqhash *seqhashCreate (int k, int w, int seed)
{
  assert (sizeof (U64) == 8) ;
  Seqhash *sh = new0 (1, Seqhash) ;
  sh->k = k ; if (k < 1 || k >= 32) die ("seqhash k %d must be between 1 and 32\n", k) ;
  sh->w = w ; if (w < 1) die ("seqhash w %d must be positive\n", w) ;
  sh->seed = seed ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  srandom (seed) ;
  sh->factor1 = (random() << 32) | random() | 0x01 ;
  sh->shift1 = 64 - 2*k ;
  sh->factor2 = (random() << 32) | random() | 0x01 ;
  sh->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { sh->patternRC[i] = (3-i) ; sh->patternRC[i] <<= 2*(k-1) ; }
  return sh ;
}

#include <stdio.h>

void seqhashWrite (Seqhash *sh, FILE *f)
{ if (fwrite ("SQHSHv2",8,1,f) != 1) die ("failed to write seqhash header") ;
  if (fwrite (sh,sizeof(Seqhash),1,f) != 1) die ("failed to write seqhash") ;
}

Seqhash *seqhashRead (FILE *f)
{ Seqhash *sh = new (1, Seqhash) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read seqhash header") ;
  if (strcmp (name, "SQHSHv2")) die ("seqhash read mismatch") ;
  if (fread (sh,sizeof(Seqhash),1,f) != 1) die ("failed to read seqhash") ;
  return sh ;
}

void seqhashReport (Seqhash *sh, FILE *f)
{ fprintf (f, "SH k %d  w/m %d  s %d\n", sh->k, sh->w, sh->seed) ; }

/************** basic hash functions *************/

static inline U64 hashRC (SeqhashIterator *si, bool *isForward)
{ U64 hashF = kHash (si->sh, si->h) ;
  U64 hashR = kHash (si->sh, si->hRC) ;
#ifdef DEBUG
  printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", si->h, si->hRC, hashF, hashR) ;
#endif
  if (hashF < hashR) { *isForward = true ; return hashF ; }
  else { *isForward = false ; return hashR ; }
}

static inline U64 advanceHashRC (SeqhashIterator *si, bool *isForward)
{ Seqhash *sh = si->sh ;
  if (si->s < si->sEnd)
    { si->h = ((si->h << 2) & sh->mask) | *(si->s) ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*(si->s)] ;
      ++si->s ;
      return hashRC (si, isForward) ;
    }
  else
    return U64MAX ;
}

/************ iterators to run across a sequence, returning (a subset of) hashes ***********/

/*************** this basic one returns all the hashes *********************/
/*************** and its creator is used as a base by the others ***********/

SeqhashIterator *seqhashIterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashIterator *si = new0 (1, SeqhashIterator) ;
  si->sh = sh ;
  si->s = s ; si->sEnd = s + len ;
  si->hash = new0 (sh->w, U64) ;
  si->isForward = new0 (sh->w, bool) ;
  if (len < sh->k)
    si->isDone = true ; // edge case
  else
    { int i ;			/* preinitialise the hashes for the first kmer */
      for (i = 0 ; i < sh->k ; ++i, ++si->s)
	{ si->h = (si->h << 2) | *si->s ;
	  si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*(si->s)] ;
	}
      *si->hash = hashRC (si, si->isForward) ;
    }
  return si ;
}

bool seqhashNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

  if (kmer) *kmer = si->h ;
  if (pos) *pos = si->iStart ;
  if (isF) *isF = *si->isForward ;

  if (si->s >= si->sEnd)
    si->isDone = true ;
  else
    { *si->hash = advanceHashRC (si, si->isForward) ;
      ++si->iStart ;
    }
  
  return true ;
}

/*************** this one returns all the hashes *********************/

SeqhashIterator *minimizerIterator (Seqhash *sh, char *s, int len)
{
  SeqhashIterator *si = seqhashIterator (sh, s, len) ;
  if (si->isDone) return si ;
  
  /* store first w hashes in hash and set ->iMin */
  si->min = si->hash[0] ;
  si->iMin = 0 ;
  int i ;
  for (i = 1 ; i < sh->w ; ++i, ++si->s)
    { si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      if (si->hash[i] < si->min) { si->min = si->hash[i] ; si->iMin = i ; }
    }

  return si ;
}

bool minimizerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF) /* returns u,pos,isF */
{
  if (si->isDone) return false ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, iMin %d\n", si->base, si->iStart, si->iMin) ;
  int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
  printf ("\n") ;
#endif

  if (kmer) *kmer = si->hash[si->iMin] ;
  if (pos) { *pos = si->base + si->iMin ; if (si->iMin < si->iStart) *pos += si->sh->w ; }
  if (isF) *isF = si->isForward[si->iMin] ;
  if (si->s >= si->sEnd) { si->isDone = true ; return true ; }

  int i ;	       	/* next update hash splitting into two cases */
  if (si->iMin >= si->iStart)
    for (i = si->iStart ; i <= si->iMin ; ++i)
      si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
  else
    { for (i = si->iStart ; i < si->sh->w ; ++i)
	si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      si->base += si->sh->w ;
      for (i = 0 ; i <= si->iMin ; ++i)
	si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
    }
  si->iStart = si->iMin + 1 ;
  if (si->iStart == si->sh->w) { si->iStart = 0 ; si->base += si->sh->w ; }

  /* finally find new min to set up for next call */
  if (si->hash[si->iMin] != U64MAX) /* there was a full new window */
    si->min = U64MAX ;
  else				/* otherwise, keep the last min */
    si->iMin = -1 ;
  for (i = 0 ; i < si->sh->w ; ++i)
    if (si->hash[i] < si->min) { si->min = si->hash[i] ; si->iMin = i ; }
  if (si->iMin == -1)		/* our old min was not beaten - we are done */
    si->isDone = true ;
  
  return true ;
}

/************ same for closed syncmer ***********/

SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len)
{
  SeqhashIterator *si = seqhashIterator (sh, s, len) ;
  if (len < sh->w + sh->k) si->isDone = true ; // because we are looking for w-mers not k-mers here
  if (si->isDone) return si ;
    
  /* store first w hashes in hash and set ->min */
  si->min = si->hash[0] ;
  int i ;
  for (i = 1 ; i < sh->w ; ++i)
    { si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }

  // si->iStart = 0 ; // from initialisation
  if (si->hash[0] == si->min || si->hash[sh->w-1] == si->min) return si ; // we are done
  while (true)
    { U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return si ; }
      si->hash[si->iStart++] = x ;
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ; return si ; }
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	return si ;
    }
  die ("syncmer initialisation failure") ;
}

bool syncmerNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, min %" PRIx64 "\n", si->base, si->iStart, si->min) ;
  int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
  printf ("\n") ;
#endif

  if (kmer) *kmer = si->hash[si->iStart] ;
  if (pos) *pos = si->base + si->iStart ;
  if (isF) *isF = si->isForward[si->iStart] ;

  if (si->hash[si->iStart] == si->min) // need to find new min - could use a heap, but not so bad to search here
    { int i ;
      si->min = si->hash[si->iStart] = U64MAX ;
      for (i = 0 ; i < si->sh->w ; ++i) if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }
  
  while (true) // move forwards to the next minimum
    { U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
      si->hash[si->iStart++] = x ;
      if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ;
	  //	  printf (" syncmerNext   at_end %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
	  return true ;
	}
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	{
	  // printf (" syncmerNext at_start %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
	  return true ;
	}
    }
}

/************** and for modimizers **********************/

SeqhashIterator *modIterator (Seqhash *sh, char *s, int len)
{
  SeqhashIterator *si = seqhashIterator (sh, s, len) ;
  if (si->isDone) return si ;
  
  U64 hash = si->hash[0] ;
  while ((hash % sh->w) && si->s < si->sEnd)
    { hash = advanceHashRC (si, si->isForward) ; ++si->iMin ; ++si->s ; }
  if (!(hash % sh->w)) *si->hash = hash ;
  else si->isDone = true ;

  return si ;
}

bool modNext (SeqhashIterator *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

  if (kmer) { if (*si->isForward) *kmer = si->h ; else *kmer = si->hRC ; }
  if (pos) *pos = si->iMin ; if (isF) *isF = *si->isForward ;

  if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
    
  U64 u = advanceHashRC (si, si->isForward) ; ++si->iMin ; ++si->s ;
  int w = si->sh->w ;
  while ((u % w) && si->s < si->sEnd)
    { u = advanceHashRC (si, si->isForward) ; ++si->iMin ; ++si->s ; }
  if (!(u % w)) *si->hash = u ;
  else si->isDone = true ;

  return true ;
}

/*********************************************/

char *seqString (U64 kmer, int len)
{
  static char trans[4] = { 'a', 'c', 'g', 't' } ;
  static char buf[33] ;
  assert (len <= 32) ;
  buf[len] = 0 ;
  while (len--) { buf[len] = trans[kmer & 0x3] ; kmer >>= 2 ; }
  return buf ;
}

/************** short test program, illustrating standard usage *************/

#ifdef TEST

#include "seqio.h"

int main (int argc, char *argv[])
{
  char *seq, *id ;
  int len ;
  U64 u ; int pos ; bool isF ;
  int K = 16 ;
  int W = 1023 ;
  SeqIO *sio ;

  if (argc > 1) sio = seqIOopenRead (argv[1], dna2indexConv, 0) ;
  else sio = seqIOopenRead ("-", dna2indexConv, 0) ;
  
  Seqhash *sh = seqhashCreate (K, W, 7) ;
  while (seqIOread (sio))
    { printf ("\nread sequence %s length %" PRIu64 "\n", sqioId(sio), sio->seqLen) ;
      SeqhashIterator *si = syncmerIterator (sh, sqioSeq(sio), sio->seqLen) ;
      while (syncmerNext (si, &u, &pos, &isF))
	{ printf ("\t%08llx\t%s\t%d\t%c",
		  u, seqhashString (si->sh, u), pos, isF?'F':'R') ;
	  printf ("\t%s", seqhashString (si->sh, si->h)) ;
	  printf ("\t%s\n", seqhashString (si->sh, si->hRC)) ;
	}
      seqhashIteratorDestroy (si) ;
    }
  seqIOclose (sio) ;
}

#endif

void seqhashForCompilerHappiness (Seqhash *sh, SeqhashIterator *sit)
{ seqhashDestroy (sh) ; seqhashIteratorDestroy (sit) ; }

/**************** end of file ****************/
