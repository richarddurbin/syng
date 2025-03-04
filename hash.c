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
 * Last edited: Mar  4 01:41 2025 (rd109)
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
  I64 mask ;       		/* 2**nbits-1 */
  int n ;			/* number of items stored */
  int guard ;			/* number of slots to fill before doubling */
  HashKey *keys ;		/* array of keys stored */
  int *values ;			/* array of indices */
  Array freeList ;     		/* list of removed integers that can be reused */
  int nFree ;			/* number in free list */
  int iter ;			/* used in iterator */
} TrueHash ;

static const int IS5 = (sizeof(I64)*8)/5 ;
static const int IS7 = (sizeof(I64)*8)/7 ;

static inline I64 hashFunc (TrueHash *h, HashKey hk)
{ int z = IS5 ; I64 x = hk.i, hash ;
  for (hash = x, x >>= 5 ; z-- ; x >>= 5) hash ^= x ;
  return hash & h->mask ;
}

static inline I64 deltaFunc (TrueHash *h, HashKey hk)
{ int z = IS7 ; I64 x = hk.i, delta ;
  for (delta = x, x >>= 7 ; z-- ; x >>= 7) delta ^= x ;
  return (delta & h->mask) | 0x01 ;  /* delta odd is prime relative to  2^m */
}

#define HASH_FUNC(_key) { register int z = IS5, x = _key.i ; \
		          for (hash = x, x >>= 5 ; z-- ; x >>= 5) hash ^= x ; \
			  hash &= h->mask ;				\
                        }

#define DELTA(_key) { register int z = IS7, x = _key.i ; \
		      for (delta = x, x >>= 7 ; z-- ; x >>= 7) delta ^= x ; \
		      delta = (delta & h->mask) | 0x01 ; \
                    }

static HashKey REMOVED ;

static int nCreated = 0 ;
static int nDestroyed = 0 ;
static I64 nAdded = 0 ;
static I64 nBounced = 0 ;
static I64 nFound = 0 ;
static I64 nNotFound = 0 ;

/**************** create/destroy/clear ******************/

Hash hashCreate (int n)
{
  TrueHash *h = new0 (1, TrueHash) ;

  if (sizeof(I64) != sizeof(HashKey)) die ("type size mismatch in hashCreate") ;
  REMOVED = hashInt((I64MAX-1)^I64MAX) ;
  
  if (n < 64) n = 64 ;
  --n ;
  h->nbits = 1 ;	       /* make room, be twice as big as needed */
  while (n >>= 1) ++h->nbits ; /* number of left most bit + 1 */
  h->mask = (1 << h->nbits) - 1 ;
  h->guard = (1 << (h->nbits - 1)) ;
  h->keys = new0 (1 << h->nbits, HashKey) ;
  h->values = new (1 << h->nbits, int) ;
  h->n = 0 ;
  h->freeList = arrayCreate (32, int) ;
  h->nFree = 0 ;
  ++nCreated ;
  return (Hash) h ;
}

void hashDestroy (Hash hx)
{
  TrueHash *h = (TrueHash*) hx ;
  newFree (h->keys, 1 << h->nbits, HashKey) ;
  newFree (h->values, 1 << h->nbits, int) ;
  arrayDestroy (h->freeList) ;
  newFree (h, 1, TrueHash) ;
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
  int      oldsize, newsize ;
  I64      hash, delta = 0 ;
  HashKey *oldKeys, *kp ;
  int     *oldValues, i ;

  oldsize = 1 << h->nbits ;
  ++h->nbits ;
  newsize = 1 << h->nbits ;
  h->mask = (1 << h->nbits) - 1 ;
  h->guard = (1 << (h->nbits - 1)) ;
  
  oldKeys = h->keys ;
  h->keys  = new0 (newsize, HashKey) ;
  oldValues = h->values ;
  h->values = new (newsize, int) ;

  for (i = 0, kp = oldKeys ; i < oldsize ; ++i, ++kp)
    if (kp->i && kp->i != REMOVED.i)
      { hash = hashFunc(h, *kp) ;
        while (true)
          if (!h->keys[hash].i)  /* don't need to test REMOVED */
	    { h->keys[hash] = *kp ;
	      h->values[hash] = oldValues[i] ;
	      --h->guard ;	/* NB don't need to change h->n */
	      ++nAdded ;
	      delta = 0 ;
	      break ;
	    }
	  else
            { nBounced++ ;
	      if (!delta) delta = deltaFunc (h, *kp) ;
	      hash = (hash + delta) & h->mask ;
	    }
      }

  newFree (oldKeys, oldsize, HashKey) ;
  newFree (oldValues, oldsize, int) ;
}

/************************ Searches  ************************************/

bool hashFind (Hash hx, HashKey k, int *index)
/* if found, returns index, else returns 0 */
{
  TrueHash *h = (TrueHash*) hx ;
  I64 hash, delta = 0 ;

  hash = hashFunc (h,k) ;
  while (true)
    if (h->keys[hash].i == k.i)
      { nFound++ ;
	if (index) *index = h->values[hash] - 1 ;
	return true ;
      }
    else if (!h->keys[hash].i)
      { nNotFound++ ;
	return false ;
      }
    else 
      { nBounced++ ;
	if (!delta) delta = deltaFunc (h,k) ;
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

  hash = hashFunc (h,k) ;
  while (true)
    if (!h->keys[hash].i || h->keys[hash].i == REMOVED.i)	/* free slot to fill */
      { if (!h->keys[hash].i) --h->guard ;
	h->keys[hash] = k ;
	if (h->nFree)
	  h->values[hash] = arr(h->freeList, h->nFree--, int) ;
	else
	  h->values[hash] = ++h->n  ;
	nAdded++ ;
	if (index) *index = h->values[hash] - 1 ;
	return true ;
      }
    else if (h->keys[hash].i == k.i)		/* already there */
      { ++nFound ;
	if (index) *index = h->values[hash] - 1 ;
	return false ;
      }
    else
      { nBounced++ ;
	if (!delta) delta = deltaFunc (h,k) ;
	hash = (hash + delta) & h->mask ;
      }
}
 
/************************ Removals ************************************/

bool hashRemove (Hash hx, HashKey k)
{
  TrueHash *h = (TrueHash*) hx ;
  I64 hash, delta = 0 ;

  hash = hashFunc (h,k) ;
  while (true)
    if (h->keys[hash].i == k.i)
      { h->keys[hash].i = REMOVED.i ;
	array(h->freeList, ++h->nFree, int) = h->values[hash] ;
	++nFound ;
	return true ;
      }
    else if (!h->keys[hash].i)
      { nNotFound++ ;
	return false ;
      }
    else 
      { nBounced++ ;
	if (!delta) delta = deltaFunc (h,k) ;
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
    if (h->keys[h->iter].i && h->keys[h->iter].i != REMOVED.i)
      { kp->i = h->keys[h->iter].i ;
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
