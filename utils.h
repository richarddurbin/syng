/*  File: utils.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
 * Description: includes standard system headers and own headers
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 19 05:03 2024 (rd109)
 * Created: Wed Jan  5 16:13:48 2011 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>		/* FILE etc. */
#include <stdlib.h>		/* malloc(), free(), ... notation */
#include <inttypes.h>		/* for standard size int types and their print macros */
#include <string.h>		/* memset() */
#include <limits.h>		/* INT_MAX etc. */
#include <stdbool.h>		/* bool, true, false */
#include <assert.h>

#ifndef UTILS_DEFINED
#define UTILS_DEFINED

typedef int8_t I8 ;
const static I8 I8MAX = 0x7f ;
typedef int16_t I16 ;
const static I16 I16MAX = 0x7fff ;
typedef int32_t I32 ;
const static I32 I32MAX = 0x7fffffff ;
typedef long long I64 ;
const static I64 I64MAX = 0x7fffffffffffffff ;

typedef uint8_t U8 ;
const static U8 U8MAX = 0xff ;
typedef uint16_t U16 ;
const static U16 U16MAX = 0xffff ;
typedef uint32_t U32 ;
const static U32 U32MAX = 0xffffffff ;
typedef unsigned long long U64 ;
const static U64 U64MAX = 0xffffffffffffffff ;
#endif

void die  (char *format, ...) ;
void warn (char *format, ...) ;

void *myalloc  (size_t size) ;
void *mycalloc (size_t number, size_t size) ;
void *myresize (void* x, size_t nOld, size_t nNew, size_t size) ;
void  myfree   (void* x, size_t size) ;
#define	new(n,type)                 (type*)myalloc((n)*sizeof(type))
#define	new0(n,type)	            (type*)mycalloc((n),sizeof(type))
#define newResize(x,nOld,nNew,type) (type*) myresize((x),(nOld),(nNew),sizeof(type))
#define newDouble(x,n,type)         myresize((x),(n),2*(n),sizeof(type)), (n) = 2*(n)
#define	newFree(x,n,type)           myfree((x),(n)*sizeof(type))

void  storeCommandLine (int argc, char *argv[]) ;
char *getCommandLine (void) ;

char *fgetword (FILE *f) ;	/* not threadsafe */
FILE *fzopen (const char* path, const char* mode) ; /* will open gzip files silently */
FILE *fopenTag (char* root, char* tag, char* mode) ;
char *fnameTag (char* root, char* tag) ; /* utility to return name as used by fopenTag() */

void timeUpdate (FILE *f) ;	/* print time usage since last call to file */
void timeTotal (FILE *f) ;	/* print full time usage since first call to timeUpdate */

/************************/
