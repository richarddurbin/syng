/*  File: arraysub.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1989-
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
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              Arbitrary length arrays
 *              These functions are declared in array.h
 * Exported functions:
 *              See Header file: array.h (includes lots of macros)
 * HISTORY:
 * Last edited: Sep 29 18:36 2024 (rd109)
 * * Aug 31 22:04 2024 (rd109): added ArrayPool
 * * Feb 14 11:21 2011 (rd): modified in 2009/10 by RD for stand-alone use
 * Created: Thu Dec 12 15:43:25 1989 (mieg)
 *-------------------------------------------------------------------
 */

#include "array.h"

/********** Array : class to implement variable length arrays **********/

static U64 totalAllocatedMemory = 0 ;
static U64 totalNumberCreated = 0 ;
static U64 totalNumberActive = 0 ;
static Array reportArray = 0 ;

#define arrayExists(a) ((a) && (a)->magic == ARRAY_MAGIC)

Array uArrayCreate (U64 n, U64 size)
{
  Array a = (Array) mycalloc (1, sizeof (struct ArrayStruct)) ;
  static bool isFirst = true ;

  if (isFirst)
    { isFirst = false ;
#ifdef ARRAY_REPORT
      if (ARRAY_REPORT_MAX) reportArray = arrayCreate (512, Array) ;
#endif
    }
  if (size <= 0) die ("negative size %d in uArrayCreate", size) ;
  if (n < 1)
    n = 1 ;
  totalAllocatedMemory += n * size ;

  a->magic = ARRAY_MAGIC ;
  a->base = mycalloc (n, size) ;
  a->dim = n ;
  a->max = 0 ;
  a->size = size ;
  ++totalNumberActive ;
#ifdef ARRAY_REPORT
  a->id = ++totalNumberCreated ;
  if (reportArray)
    { if (a->id < ARRAY_REPORT_MAX)
	array (reportArray, a->id, Array) = a ;
      else
	{ arrayDestroy (reportArray) ;
	  reportArray = 0 ;
	}
    }
#endif
  return a ;
}

/**************/

Array uArrayReCreate (Array a, U64 n, U64 size)
{
  if (!arrayExists(a))
    return  uArrayCreate (n, size) ;

  if (a->size != size)
    die ("type size mismatch in arrayReCreate: size %d != a->size %d", size, a->size) ;

  if (n < 1) n = 1 ;

  if (a->dim < n || (a->dim - n)*size > (1 << 20) ) /* free if save > 1 MB */
    { totalAllocatedMemory -= a->dim * size ; 
      myfree (a->base, a->dim*size) ;
      a->dim = n ;
      totalAllocatedMemory += a->dim * size ;
      a->base = mycalloc (n, size) ;
    }
  else
    memset (a->base, 0, n*size) ;

  a->max = 0 ;
  return a ;
}

/**************/

void arrayDestroy (Array a)
{
  if (!arrayExists (a))
    die ("arrayDestroy called on bad array %lx", (long unsigned int) a) ;

  totalAllocatedMemory -= a->dim * a->size ; 
  totalNumberActive-- ;
#ifdef ARRAY_REPORT
  if (reportArray)
    arr(reportArray, a->id, Array) = 0 ;
#endif
  a->magic = 0 ;
  myfree (a->base, a->dim*a->size) ;
  newFree (a, 1, struct ArrayStruct) ;
}

/**************/

Array arrayCopy (Array a) 
{
  Array new ;

  if (!arrayExists (a)) 
    die ("arrayCopy called on bad array %lx", (long unsigned int) a) ;
 
  new = uArrayCreate (a->dim, a->size) ;
  memcpy (new->base, a->base, a->dim * a->size) ;
  new->max = a->max ;
  return new;
}

/******************************/

void arrayExtend (Array a, U64 n) 
{
  if (!arrayExists (a))
    die ("arrayExtend called on bad array %lx", (long unsigned int) a) ;

  if (n < a->dim)
    return ;

  U64 oldDim = a->dim ;
  totalAllocatedMemory -= a->dim * a->size ;
  if (a->dim*a->size < 1 << 23)	/* 8MB */
    a->dim *= 2 ;
  else
    a->dim += 1024 + ((1 << 23) / a->size) ;
  if (n >= a->dim)
    a->dim = n + 1 ;

  totalAllocatedMemory += a->dim * a->size ;

  a->base = newResize (a->base, oldDim*a->size, a->dim*a->size, char) ;

  return;
}

/***************/

char *uArray (Array a, U64 i)
{
#ifdef ARRAY_CHECK
  if (!arrayExists (a)) die ("array() called on bad array %lx", (long unsigned int) a) ;
  if (i < 0) die("referencing array element %d < 0", i) ;
#endif
  if (i >= a->max)
    { if (i >= a->dim) arrayExtend (a,i) ;
      a->max = i+1 ;
    }
  return a->base + i*a->size ;
}

char *uArrCheck (Array a, U64 i)
{
  if (!arrayExists (a)) die ("array() called on bad array %lx", (long unsigned int) a) ;
  if (i < 0) die("referencing array element %d < 0", i) ;
  if (i >= a->max) die ("referencing array element %d >= max %d", i, a->max) ;

  return a->base + i*a->size ;
}

/***************/

char    *uArrayBlock (Array a, U64 i, U64 n)
{
#ifdef ARRAY_CHECK
  if (!arrayExists (a)) die ("arrayBlock() called on bad array %lx", (long unsigned int)a) ;
  if (i < 0) die ("arrayBlock() referencing array element %d < 0", i) ;
#endif
  if (i+n >= a->max)
    { if (i+n >= a->dim) arrayExtend (a,i+n) ;
      a->max = i+n+1 ;
    }
  return a->base + i*a->size ;
}

/*************** writing and reading to files ***************/

bool    arrayWrite (Array a, FILE *f)
{
  if (fwrite (a, sizeof(struct ArrayStruct), 1, f) != 1) return false ;
  if (fwrite (a->base, a->size, a->dim, f) != a->dim) return false ;
  return true ;
}

Array   arrayRead (FILE *f)
{
  Array a = new (1, struct ArrayStruct) ;
  if (fread (a, sizeof(struct ArrayStruct), 1, f) != 1) return 0 ;
  a->base = myalloc (a->size*a->dim) ;
  if (fread (a->base, a->size, a->dim, f) != a->dim)
    { myfree (a->base, a->size*a->dim) ;
      newFree (a, 1, struct ArrayStruct) ;
      return 0 ;
    }
#ifdef ARRAY_REPORT
  a->id = ++totalNumberCreated ;
  if (reportArray)
    { if (a->id < ARRAY_REPORT_MAX)
	array (reportArray, a->id, Array) = a ;
      else
	{ arrayDestroy (reportArray) ;
	  reportArray = 0 ;
	}
    }
#endif
  return a ;
}

/***********************************************/
       /* Finds Entry s from Array  a
        * sorted in ascending order of order()
        * If found, returns true and sets *ip
        * if not, returns false and sets *ip one step left
        */

/* bool arrayFind(Array a, void *s, U64 *ip, U64 (* order)(void*, void*)) */
bool arrayFind(Array a, void *s, U64 *ip, ArrayOrder *order)
{
  U64 ord ;
  U64 i = 0 , j, k ;

  if (!arrayExists (a)) 
    die ("arrayFind called on bad array %lx", (long unsigned int) a) ;

  j = arrayMax(a) ;
  if (!j || (ord = order(s,uArray(a,0))) < 0)
    { if (ip) *ip = -1 ;
      return false ;
    }   /* not found */

  if (ord == 0)
    { if (ip) *ip = 0 ;
      return true ;
    }

  if ((ord = order(s,uArray(a,--j))) > 0)
    { if (ip) *ip = j ;
      return false ;
    }
  
  if (ord == 0)
    { if (ip) *ip = j ;
      return true ;
    }

  while(true)
    { k = i + ((j-i) >> 1) ; /* midpoint */
      if ((ord = order(s, uArray(a,k))) == 0)
	{ if (ip) *ip = k ;
	  return true ;
	}
      if (ord > 0) i = k ;
      else j = k ;
      if (i == (j-1))
        break ;
    }
  if (ip)
    *ip = i ;
  return false ;
}

/**************************************************************/
       /* Removes Entry s from Array  a
        * sorted in ascending order of order()
        */

bool arrayRemove (Array a, void * s, int (* order)(const void*, const void*))
{
  U64 i;

  if (!arrayExists (a))
    die ("arrayRemove called on bad array %lx", (long unsigned int) a) ;

  if (arrayFind(a, s, &i,order))
    {
      /* memcpy would be faster but regions overlap
       * and memcpy is said to fail with some compilers
       */
      char *cp = uArray(a,i), *cq = cp + a->size ;
      U64 j = (arrayMax(a) - i)*(a->size) ;
      while (j--)
	*cp++ = *cq++ ;

      arrayMax(a)-- ;
      return true ;
    }
  else

    return false ;
}

/**************************************************************/
       /* Insert Segment s in Array  a
        * in ascending order of s.begin
        */

bool arrayInsert(Array a, void * s, int (*order)(const void*, const void*))
{
  U64 i, j, arraySize;

  if (!arrayExists (a))
    die ("arrayInsert called on bad array %x", (long unsigned int)a) ;

  if (arrayFind(a, s, &i,order))
    return false ;  /* no doubles */
  
  arraySize = arrayMax(a) ;
  j = arraySize + 1 ;

  uArray(a,j-1) ; /* to create space */

	/* avoid memcpy for same reasons as above */
  {
    char *cp, *cq ;
    U64 k ;
    
    if (arraySize > 0)
      { cp = uArray(a,j - 1) + a->size - 1 ;
	cq = cp - a->size ;
	k = (j - i - 1)*(a->size) ;
	while (k--)
	  *cp-- = *cq-- ;
      }
    
    cp = uArray(a,i+1) ; 
    cq = (char *) s ; 
    k = a->size ;
    while (k--)
      *cp++ = *cq++ ;
  }
  return true ;
}

/**************/

void arrayCompress(Array a)
{
  U64 i, j, k , as ;
  char *x, *y, *ab ;

  if (!arrayExists (a))
    die ("arrayCompress called on bad array %lx", (long unsigned int) a) ;

  if (arrayMax(a) < 2)
    return ;

  ab = a->base ; 
  as = a->size ;
  for (i = 1, j = 0 ; i < arrayMax(a) ; i++)
    { x = ab + i * as ; y = ab + j * as ;
      for (k = a->size ; k-- ;)		
	if (*x++ != *y++) 
	  goto different ;
      continue ;
      
    different:
      if (i != ++j)
	{ x = ab + i * as ; y = ab + j * as ;
	  for (k = a->size ; k-- ;)	 
	    *y++ = *x++ ;
	}
    }
  arrayMax(a) = j + 1 ;
}

/**************/

U64 arrayReportMark (void)
{
  if (reportArray)
    return arrayMax (reportArray) ;
  else
    return 0 ;
}

/**************/

void arrayReport (U64 j)
{
  U64 i ;
  Array a ;

  fprintf(stderr, "Array report: %llu created, %llu active, %llu MB allocated\n",   
	  totalNumberCreated, totalNumberActive, totalAllocatedMemory/(1024*1024)) ;

  if (reportArray)
    { i = arrayMax (reportArray) ;
      while (i-- && i > j)
	{ a = arr (reportArray, i, Array) ;
	  if (arrayExists(a))
	    fprintf (stderr, " array %llu size %llu max %llu\n", i, a->size, a->max) ;
	}
    }
}

/**************/

void arrayStatus (U64 *nmadep, U64 *nusedp, U64 *memAllocp, U64 *memUsedp)
{ 
  U64 i ;
  Array a ;

  *nmadep = totalNumberCreated ; 
  *nusedp = totalNumberActive ;
  *memAllocp = totalAllocatedMemory ;
  *memUsedp = 0 ;

  if (reportArray)
    for (i = 0 ; i < arrayMax(reportArray) ; ++i)
      if (arrayExists (a = arr(reportArray, i, Array)))
	*memUsedp += a->max * a->size ;
}

/************************  end of file ********************************/
 
