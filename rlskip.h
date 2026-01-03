/*  File: rkskip.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: interface for run-length skip lists
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  1 18:03 2026 (rd109)
 * Created: Sun Nov 30 21:20:56 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h" // for type definitions

#ifndef RLSKIP_DEFINED
typedef void* RLSkip ;
#endif

int    rlsDestroy (RLSkip rls) ;
RLSkip rlsBuildFromI64 (int n, I64 *sym, I64 *runLen) ;    // n (sym, runLen) pairs
int    rlsCount   (RLSkip rls, U32 symbol) ;          // how many of symbol
int    rlsRank    (RLSkip rls, U32 i, U32 symbol) ;   // how many symbol up to and before i
int    rlsFind    (RLSkip rls, U32 i, U32 *symbol) ;  // gives symbol at i and returns rank of it
int    rlsAdd     (RLSkip *rls, U32 i, U32 symbol) ;  // add symbol at i and return rank
int    rlsLock    (RLSkip *rls) ;                     // pack memory and reduce
char  *rlsError   (RLSkip *rls, int err) ;            // returns error text for err

#define RLS_EMPTY 		-1 // all error codes must be negative
#define RLS_SYMBOL_OVERFLOW	-2
#define RLS_INDEX_OVERFLOW	-3
#define RLS_LOCKED		-4

// end of file
