/*  File: rskip.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: interface for run-length skip lists
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 27 23:37 2026 (rd109)
 * Created: Sun Nov 30 21:20:56 2025 (rd109)
 *-------------------------------------------------------------------
 */

// this implements run-length encoded skip lists
// we require the symbols to be integers from 0 to nSym-1, with no unused symbols
// there are two states: dynamic and locked
// dynamic rskips can be added to, whereas locked rskips can not be added to
// count and rank operations on locked rskips are log(r_s) where r_s is the count of symbol s

#include "utils.h" // for U32, I64 type definitions (and general Durbin utilities)

#ifndef RSKIP_DEFINED
typedef void* Rskip ;
#endif

Rskip rsCreateRaw (void) ;     // will create an empty raw rSkip
Rskip rsCreate    (int nSym, I32 *symbol) ; // optional preload nSym symbols, none if nSym == 0
int   rsDestroy   (Rskip rs) ; // destroy any Rskip, cleaning up memory

int   rsAdd    (Rskip *rs, U32 k, I32 symbol) ; // add symbol at k and return rank (see below)
// if rs was created with rsCreateRaw() then symbol must be in 0..nsym
// else symbol can be anything - linear search on symbol, adding at end if necessary
int   rsLength (Rskip rs) ;                     // total length = sum of runlengths
int   rsCount  (Rskip rs, I32 symbol) ;         // how many of symbol
int   rsRank   (Rskip rs, U32 k, I32 symbol) ;  // how many symbol up to (not including) k
int   rsFind   (Rskip rs, U32 k, I32 *symbol) ; // gives symbol at k and returns rank of it
int   rsNsym   (Rskip rs) ; // number of symbols currently used - constant time

// same for use with syng
Rskip rsCreateSyng (int nSym, I32 *symbol, U32 *offset) ; // no preloading if nSym == 0
Rskip rsBuildFixedSyng (int nSym, I32 *symbol, U32 *offset, int nRun, I64 *iSym, I64 *runLen) ;
Rskip rsBuildDynamicSyng (int nSym, I32 *symbol, U32 *offset, int nRun, I64 *iSym, I64 *runLen) ;
// use I64 for iSym and runLen here so can natively pass to/from ONEcode arrays
int   rsAddSyng    (Rskip *rs, U32 k, I32 symbol, U32 offset) ;
int   rsCountSyng  (Rskip rs, I32 symbol, U32 offset) ;
int   rsRankSyng   (Rskip rs, U32 k, I32 symbol, U32 offset) ;
int   rsFindSyng   (Rskip rs, U32 k, I32 *symbol, U32 *offset) ;
// rsLength() and rsNsym() work fine on Syng Rskips

// error codes - all must be negative
#define RS_ERROR_EMPTY 		-1
#define RS_ERROR_SYMBOL_OVERFLOW	-2
#define RS_ERROR_INDEX_OVERFLOW	-3
#define RS_ERROR_LOCKED		-4
char  *rsErrorText (void) ;                        // returns text for last error this thread

// end of file
