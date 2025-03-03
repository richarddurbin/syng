/*  File: hash.h
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
 *-------------------------------------------------------------------
 * Description:
      - can hash integers, floats, pointers
      - returns an integer to be used as index into a local array
      - index of first item to be stored is 1, second is 2 etc.
      - support removal of objects, with free-list to reuse indices of removed objects
 * Exported functions:
 * HISTORY:
 * Last edited: Mar  1 23:03 2025 (rd109)
 * * Jan 24 07:01 2019 (rd109): restructured so hashAdd() and hashFind() return bool 
         and the index value is returned via a pointer, as for DICT.
 * Created: Thu Jan 13 10:57:20 2011 (rd)
 *-------------------------------------------------------------------
 */

#ifndef HASH_DEFINED
#define HASH_DEFINED

#include "utils.h"

typedef void* Hash ;
typedef union {I64 i ; double f ; struct { I32 ia, ib ; } ; const void* p ;} HashKey ;

/* define keys so as to allow ourselves to hash value 0 */
static inline HashKey hashInt (const I64 x)
	{ HashKey hk ; hk.i = x^I64MAX ; return hk ; }
static inline HashKey hashInt2 (const I32 x, const I32 y)
	{ HashKey hk ; hk.ia = x ; hk.ib = y ; hk.i ^= I64MAX ; return hk ; }
static inline HashKey hashFloat (const double x)
	{ HashKey hk ; hk.f = x ; hk.i ^= I64MAX ; return hk ; }
static inline HashKey hashPtr (const void *x)
	{ HashKey hk ; hk.p = x ; return hk ; }

Hash hashCreate (int n) ;
void hashDestroy (Hash h) ;
void hashClear (Hash h) ;
bool hashAdd  (Hash h, HashKey k, int *index) ; /* return true if added, fill index if non-zero */
bool hashFind (Hash h, HashKey k, int *index) ; /* return true if found, fill index if non-zero */
bool hashRemove (Hash h, HashKey k) ; /* if found, remove and return true, else false */
int  hashCount (Hash h) ;	      /* number of keys stored in hash */

/* iterator to get all key-value pairs, in arbitrary order */
/* note that if nothing is removed then values are incrementing from 1 to n_added */
/* note also that adding while iterating can invalidate the iterator */
void hashInitIterator (Hash h) ;
bool hashNextKeyValue (Hash h, HashKey *kp, int *ip) ;

void hashStats (void) ;		/* overall stats on package/performance */

#endif

/*********** end of file ************/
