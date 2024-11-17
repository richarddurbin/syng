/*  File: array.h
 *  Authors: Richard Durbin (rd@sanger.ac.uk) and Jean Thierry-Mieg
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1989-
 * -------------------------------------------------------------------
 * This file is derived from the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * Description: header for array.c
 * Exported functions:
 *              the Array type and associated macros and functions
 * HISTORY:
 * Last edited: May 29 12:51 2023 (rd109)
 *-------------------------------------------------------------------
 */

#ifndef ARRAY_DEFINED
#define ARRAY_DEFINED

#include "utils.h"

/* #define ARRAY_CHECK either here or in a single file to
   check the bounds on arr() and arrp() calls
   if defined here can remove from specific C files by defining ARRAY_NO_CHECK
*/

/* #define ARRAY_CHECK */

typedef struct ArrayStruct
  { U64   magic ;
    char* base ;    /* char* since need to do pointer arithmetic in bytes */
    U64   dim ;     /* length of alloc'ed space in number of elements */
    U64   size ;    /* size of elements in bytes */
    U64   max ;     /* 1+largest element accessed via array() = number of active elements */
#ifdef ARRAY_REPORT
    U64   id ;      /* unique identifier */
#endif
  } *Array ;
 
    /* NB we need the full definition for arr() for macros to work
       do not use it in user programs - it is private.
    */

#define ARRAY_MAGIC 8918274

Array   uArrayCreate (U64 n, U64 size) ;
#define arrayCreate(n,type)	         uArrayCreate(n,sizeof(type))
Array   uArrayReCreate (Array a, U64 n, U64 size) ;
#define arrayReCreate(a,n,type)	         uArrayReCreate(a,n,sizeof(type))
void    arrayDestroy (Array a) ;
Array	arrayCopy (Array a) ;
void    arrayExtend (Array a, U64 n) ;

     /* array() and arrayp() will extend the array if necessary */

char    *uArray (Array a, U64 index) ;
#define array(ar,i,type)	(*(type*)uArray(ar,i))
#define arrayp(ar,i,type)	((type*)uArray(ar,i))
char    *uArrayBlock (Array a, U64 i, U64 n) ;
#define arrayBlock(ar,i,n,type) ((type*)uArrayBlock(ar,i,n)) /* use when memset(), fread() etc multiple items */

     /* only use arr() when there is no danger of needing expansion */

#if (defined(ARRAY_CHECK) && !defined(ARRAY_NO_CHECK))
char    *uArrCheck (Array a, U64 index) ;
#define arr(ar,i,type)	(*(type*)uArrCheck(ar,i))
#define arrp(ar,i,type)	((type*)uArrCheck(ar,i))
#else
#define arr(ar,i,type)	((*(type*)((ar)->base + (i)*(ar)->size)))
#define arrp(ar,i,type)	(((type*)((ar)->base + (i)*(ar)->size)))
#endif /* ARRAY_CHECK */

bool    arrayWrite (Array a, FILE *f) ;
Array   arrayRead (FILE *f) ;	/* returns 0 if fails */

#define arrayMax(ar)  ((ar)->max)

	/* JTM's package to hold sorted arrays of ANY TYPE */
typedef int ArrayOrder(const void*, const void*) ;             /* call back function prototype for sorting arrays */
#define arraySort(a,order)  qsort((a)->base, (a)->max, (a)->size, order)
bool    arrayInsert(Array a, void * s, ArrayOrder *order);
bool    arrayRemove(Array a, void * s, ArrayOrder *order);
void    arrayCompress(Array a) ;
bool    arrayFind(Array a, void *s, U64 *ip, ArrayOrder *order);

#ifdef ARRAY_REPORT
	/* status and memory monitoring */
#define ARRAY_REPORT_MAX 0	/* set to maximum number of arrays to keep track of */
void    arrayStatus (U64 *nmadep, U64* nusedp, U64 *memAllocp, U64 *memUsedp) ; /* memUsed only up to REPORT_MAX */
U64     arrayReportMark (void) ; /* returns current array number */
void    arrayReport (int j) ;	/* write stderr about all arrays since j */
#endif

#endif

/**************************** End of File ******************************/
