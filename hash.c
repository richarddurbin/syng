/*  File: hash.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
 * Description: based on acedb code from Jean Thierry-Mieg and Richard Durbin 1999-2004
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 28 00:55 2024 (rd109)
 * Created: Fri Jan  7 09:20:25 2011 (rd)
 *-------------------------------------------------------------------
 */


#include "hash.h"
#include "array.h"

/* Originally grabbed from Steve Omohundro's sather code by Richard Durbin.

   Entirely rewritten in early 1990s by Jean Thierry-Mieg as acedb ass* package, 
   using bouncing by relative primes and deletion flagging.  Assumed non-0, non-minus-one
   pointers as keys and values.

   Again rewritten by RD in 2011 as stand alone hash* package, accepting a union 
   int/float/pointer as input, and generating consecutive positive integer outputs to
   be used as array indices, i.e. similar to DICT package, generating extra indirection 
   if true type is an integer, but flexible.
   Empty and removed keys are now INT_MAX and INT_MAX-1 for integers and floats.
*/

typedef struct {
  int nbits ;			/* power of 2 = size of arrays - 1 */
  unsigned int mask ;		/* 2**nbits-1 */
  int n ;			/* number of items stored */
  int guard ;			/* number of slots to fill before doubling */
  I64 *keys ;		/* array of keys stored */
  int *values ;			/* array of indices */
  Array freeList ;     		/* list of removed integers that can be reused */
  int nFree ;			/* number in free list */
  int iter ;			/* used in iterator */
} TrueHash ;

static const int IS5 = (sizeof(I64)*8)/5 ;
static const int IS7 = (sizeof(I64)*8)/7 ;

#define HASH_FUNC(_key) { register int z = IS5, x = _key.i ; \
		  for (hash = x, x >>= 5 ; z-- ; x >>= 5) hash ^= x ; \
		  hash &= h->mask ; \
		}

#define DELTA(_key)   { register int z = IS7, x = _key.i ; \
		       for (delta = x, x >>= 7 ; z-- ; x >>= 7) delta ^= x ; \
		       delta = (delta & h->mask) | 0x01 ; \
		   }  /* delta odd is prime relative to  2^m */

static int REMOVED ;

static int nCreated = 0 ;
static int nDestroyed = 0 ;
static I64 nAdded = 0 ;
static I64 nBounced = 0 ;
static I64 nFound = 0 ;
static I64 nNotFound = 0 ;

/**************** create/destroy/clear ******************/

Hash hashCreate (int n)
{
  TrueHash *h = (TrueHash *) myalloc (sizeof (TrueHash)) ;

  if (sizeof(I64) != sizeof(HashKey)) die ("type size mismatch in hashCreate") ;
  REMOVED = (I64MAX-1)^I64MAX ; // was in initialiser, but Linux does not like that??
  
  if (n < 64) n = 64 ;
  --n ;
  h->nbits = 1 ;	       /* make room, be twice as big as needed */
  while (n >>= 1) ++h->nbits ; /* number of left most bit + 1 */
  h->mask = (1 << h->nbits) - 1 ;
  h->guard = (1 << (h->nbits - 1)) ;
  h->keys = (I64*) myalloc (sizeof(I64)*(1 << h->nbits)) ;
  memset (h->keys, 0, sizeof(I64)*(1 << h->nbits)) ;
  h->values = (int*) myalloc (sizeof(int)*(1 << h->nbits)) ;
  h->n = 0 ;
  h->freeList = arrayCreate (32, int) ;
  h->nFree = 0 ;
  ++nCreated ;
  return (Hash) h ;
}

void hashDestroy (Hash hx)
{
  TrueHash *h = (TrueHash*) hx ;
  free (h->keys) ;
  free (h->values) ;
  arrayDestroy (h->freeList) ;
  free (h) ;
  ++nDestroyed ;
}

void hashClear (Hash hx)
{ 
  TrueHash *h = (TrueHash*) hx ;
  h->n = 0 ;
  memset (h->keys, 0, sizeof(I64)*(1 << h->nbits)) ;
  h->freeList = arrayReCreate (h->freeList, 32, int) ;
}

/********************/

static void hashDouble (TrueHash *h)
{
  int oldsize, newsize ;
  I64 hash, delta = 0 ;
  I64 *oldKeys, *kp ;
  int *oldValues, i ;
  HashKey hk ;

  oldsize = 1 << h->nbits ;
  ++h->nbits ;
  newsize = 1 << h->nbits ;
  h->mask = (1 << h->nbits) - 1 ;
  h->guard = (1 << (h->nbits - 1)) ;
  
  oldKeys = h->keys ;
  h->keys  = (I64*) myalloc (sizeof(I64)*newsize) ;
  memset (h->keys, 0, sizeof(I64)*(1 << h->nbits)) ;
  oldValues = h->values ;
  h->values = (int*) myalloc (sizeof(int)*newsize) ;

  for (i = 0, kp = oldKeys ; i < oldsize ; ++i, ++kp)
    if (*kp && *kp != REMOVED)
      { hk.i = *kp ; HASH_FUNC(hk) ;
        while (true)
          if (!h->keys[hash])  /* don't need to test REMOVED */
	    { h->keys[hash] = *kp ;
	      h->values[hash] = oldValues[i] ;
	      --h->guard ;	/* NB don't need to change h->n */
	      ++nAdded ;
	      break ;
	    }
	  else
            { nBounced++ ;
	      if (!delta) DELTA(hk) ;
	      hash = (hash + delta) & h->mask ;
	    }
      }

  free (oldKeys) ;
  free (oldValues) ;
}

/************************ Searches  ************************************/

bool hashFind (Hash hx, HashKey k, int *index)
/* if found, returns index, else returns 0 */
{
  TrueHash *h = (TrueHash*) hx ;
  I64 hash, delta = 0 ;

  HASH_FUNC(k) ;
  while (true)
    if (h->keys[hash] == k.i)
      { nFound++ ;
	if (index) *index = h->values[hash] - 1 ;
	return true ;
      }
    else if (!h->keys[hash])
      { nNotFound++ ;
	return false ;
      }
    else 
      { nBounced++ ;
	if (!delta) DELTA(k) ;
	hash = (hash + delta) & h->mask ;
      }
}

/************************ insertions  ************************************/

     /* if already there returns false, else inserts and returns TRUE */
bool hashAdd (Hash hx, HashKey k, int *index)
{
  TrueHash *h = (TrueHash*) hx ;
  I64 hash, delta = 0 ;

  if (!h->guard)
    hashDouble (h) ;

  HASH_FUNC(k) ;
  while (true)
    if (!h->keys[hash] || h->keys[hash] == REMOVED)	/* free slot to fill */
      { if (!h->keys[hash]) --h->guard ;
	h->keys[hash] = k.i ;
	if (h->nFree)
	  h->values[hash] = arr(h->freeList, h->nFree--, int) ;
	else
	  h->values[hash] = ++h->n  ;
	nAdded++ ;
	if (index) *index = h->values[hash] - 1 ;
	return true ;
      }
    else if (h->keys[hash] == k.i)		/* already there */
      { ++nFound ;
	if (index) *index = h->values[hash] - 1 ;
	return false ;
      }
    else
      { nBounced++ ;
	if (!delta) DELTA (k) ;
	hash = (hash + delta) & h->mask ;
      }
}
 
/************************ Removals ************************************/

bool hashRemove (Hash hx, HashKey k)
{
  TrueHash *h = (TrueHash*) hx ;
  I64 hash, delta = 0 ;

  HASH_FUNC(k) ;
  while (true)
    if (h->keys[hash] == k.i)
      { h->keys[hash] = REMOVED ;
	array(h->freeList, ++h->nFree, int) = h->values[hash] ;
	++nFound ;
	return true ;
      }
    else if (!h->keys[hash])
      { nNotFound++ ;
	return false ;
      }
    else 
      { nBounced++ ;
	if (!delta) DELTA(k) ;
	hash = (hash + delta) & h->mask ;
      }
}

/********************* iterator through members of a Hash **********************/

void hashInitIterator (Hash hx)
{
  TrueHash *h = (TrueHash*) hx ;

  h->iter = -1 ;
}

bool hashNextKeyValue (Hash hx, HashKey *kp, int *ip)
{
  TrueHash *h = (TrueHash*) hx ;
  int size = 1 << h->nbits ;

  while (++h->iter < size)
    if (h->keys[h->iter] && h->keys[h->iter] != REMOVED)
      { kp->i = h->keys[h->iter] ;
	if (ip) *ip = h->values[h->iter] - 1 ;
	return true ;
      }

  return false ;	/* done if reached end of table */
}

/**********************************************************************/

int hashCount (Hash hx) { TrueHash *h = (TrueHash*) hx ; return h->n - h->nFree ; }

void hashStats (void)
{
  printf ("%d hashes (%d created, %d destroyed)\n", 
	  nCreated - nDestroyed, nCreated, nDestroyed) ;
  printf ("%lld added, %lld found, %lld bounced, %lld not found\n",
	  nAdded, nFound, nBounced, nNotFound) ;
}

/************************  end of file ********************************/
