/*  File: rskip.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: code for run-length encoded skip lists
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 12 09:01 2026 (rd109)
 * Created: Sun Nov 30 21:42:51 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include <stdarg.h>

// Use a linear model with cumulative scores until MAX_LINEAR long (default 128)
// Beyond that we make Rskip lists, using Dynamic when building, or Fixed when fixed.
// Note that this limits us to maximum 2<<32 - 1 cells in the skip list.

// We will work with a list of nodes. The first one is special, holding other information
// which is referred to by union fields max, nSym, and in some cases start and/or free.
// If rs.linear->max then Linear, else if rs.dynamic->max than Dynamic, else Fixed.

// For the skiplists, we share nodes between an overall skiplist and embedded ones for each symbol.
// The structure is different between Dynamic and Fixed.  For Fixed, the skip list runs left to right
// with n->right pointing to the next right, and n->sRight to the next right of the same symbol.
// In this case the first column for each symbol must have maximum depth for that symbol.
// For Dynamic we maintain both left and right pointers, and the start for each symbol just points to
// the deepest node for the symbol, leftmost if there is more than one of the same depth.
// For Dynamic we store counts n->count and n->sCount, and for Fixed partial sums n->sum and n->sSum.
// For Fixed the symbol and down pointer share a slot, distinguished by n->sSum & FLAG32.

// In general we provide *Raw routines that access by index running from 0..nSym, and standard
// routines that access by I32 symbol.  The raw routines can start directly, whereas the standard
// ones need to find the index of the symbol in the symbol list, currently by linear search. So if
// possible to pre-load symbols in order of abundance; if lists are built up dynamically then usually
// high frequency symbols come early, so we get reasonable search times. I could maintain a hash
// table for this for large symbol lists...

// For use in syng (#define SYNG) we also keep for each symbol (a sync) a U16 offset, and we also need
// functionality to find the offset of a specific (sync, offset) in the list of these. In this case
// matching must be to both symbol and offset.

// Linear and Dynamic need directories to store ancillary information about the symbols etc.
// while Fixed only requires this #ifdef SYNG.  These directories follow the header. Actual data
// are in nodes running down from max-1.
// For Linear the directory contains the symbols, and #ifdef SYNG also the offsets and counts.
// For Dynamic, it uses full nodes and contains symbol, offset, count, start node (right)

#define MAX_LINEAR 128  // min 16, max 255 - define not constant to remove compiler warning!
static const int MAX_BIG    = (1<<16) - 1 ;       // largest count for Linear to store in bigCount
static const U32 FLAG32     = ((U32)1<<31) ;      // flag for bottom of column for Fixed
static const U32 MASK32     = ((U32)1<<31) - 1 ;

static int DEBUG = 0 ; // set to the entry you want to debug, -1 for general, or 0 for none

typedef union {
  struct {
    union { U8 iSym ; U8 max ; U8 free ; } ;
    union { U8 count ; U8 nSym ; U8 type ; } ; // if count == ff then take next node as U16 count
    // max, nSym in rs.linear[0], free in rs.linear[1]
  } ;
  union { U16 bigCount ; U16 offset ; } ; // offset is only used in syng directory
} Linear ; // we use two linear nodes for the header, and for each sym

typedef struct {
  U32 count ;
  union { U32 sCount ; U32 max ; U32 sum ; } ;
  union { U32 up ; U32 nSym ; } ;
  union { U32 down ; U32 free ; U32 depth ; } ; // down == 0 at bottom row
  union { U32 right ; U32 maxDepth ; } ;
  union { U32 sRight ; U32 start ; } ;
  U32 left, sLeft ;
  union { U32 before ; U32 offset ; } ; // before in nodes, offset in directory entries
  U32 sBefore ;
  union { U32 iSym ; I32 sym ; } ; // iSym in nodes, sym in directory
} Dynamic ;			   // version for building - max is number of nodes used
// nodes 1..nSym are directory entries; node[s+1].right points to actual top node for sIndex s
// node[s+1].depth = depth of top node for s, node[s+1].sum = total count of s
// node[s+1].count is used by syng code for the directory count for incoming not outgoing edges
// after building node[s+1].count == node[s+1].sum, but not during building
// start is the iSym whose column has maximal depth over all

typedef struct {
  U32 sum, sSum ; // sSum & FLAG32 indicates bottom row, in which case mask it with MASK32
  union { U32 right ; U32 max ; U32 offset ; } ;
  union { U32 sRight ; U32 start ; U32 count ; } ;  // for FIXED_SYNG directory count == sum
  union { U32 iSym ; U32 down ; U32 nSym ; I32 sym ; } ;
} Fixed ; 			// locked node

typedef union { Linear *linear ; Fixed *fixed ; Dynamic *dynamic ; } Rskip  ;

enum { LINEAR=1, LINEAR_RAW, LINEAR_SYNG, DYNAMIC, FIXED_RAW, FIXED, FIXED_SYNG } ;
static inline int rsType (Rskip rs) { return rs.linear[1].type ; }

typedef struct {
  I32 sym ;
  U16 count ;
  U16 offset ;
} LinearSyngDir ;
static inline LinearSyngDir *linearSyngDir (Rskip rs) { return (LinearSyngDir*) (rs.linear+2) ;  }
// a list of nSym of these following the header contains the directory information for syng

#define RSKIP_DEFINED

#include "rskip.h" // for interface declarations

/*************** first the error handling - a work in progress **************/

static _Thread_local char errorBuf[1024] ;

char *rsErrorText (void) { return errorBuf ; }

static int rsError (int errorCode, char *funcName, char *format, ...)
{
  va_list args ;
  int n = sprintf (errorBuf, "RS error in %s: ", funcName) ;
  va_start (args, format) ;
  switch (errorCode)
    {
    case RS_ERROR_EMPTY:
      n += sprintf (errorBuf+n, "null pointer ") ; vsprintf (errorBuf+n, format, args) ; break ;
    case RS_ERROR_SYMBOL_OVERFLOW:
      n += sprintf (errorBuf+n, "symbol overflow ") ; vsprintf (errorBuf+n, format, args) ; break ;
    case RS_ERROR_INDEX_OVERFLOW:
      n += sprintf (errorBuf+n, "index overflow ") ; vsprintf (errorBuf+n, format, args) ; break ;
    case RS_ERROR_LOCKED:
      n += sprintf (errorBuf+n, "object locked ") ; vsprintf (errorBuf+n, format, args) ; break ;
    default:
      die ("unknown RS error code %d in %s", errorCode, funcName) ;
    }
  va_end (args) ;
  return errorCode ;
}

/********************* next the basic functions **********************/

static Rskip rsNew (int type, U32 nSym, U32 max)
{
  Rskip rs ;
  switch (type)
    {
    case LINEAR: case LINEAR_RAW: case LINEAR_SYNG:
      { Linear *n = rs.linear = new0 (max, Linear) ;
        n->max = max ; n->nSym = nSym ; n[1].free = max - 1 ;
        break ;
      }
    case DYNAMIC:
      { Dynamic *n = rs.dynamic = new0 (max, Dynamic) ;
        n->max = max ; n->nSym = nSym ; n->free = max - 1 ;
        break ;
      }
    case FIXED: case FIXED_RAW: case FIXED_SYNG:
      { Fixed *n = rs.fixed = new0 (max, Fixed) ;
        n->max = max ; n->nSym = nSym ;
        break ;
      }
    }
  rs.linear[1].type = type ;
  return rs ;
}

Rskip rsCreateRaw (void) { return rsNew (LINEAR_RAW, 0, 8) ; }

Rskip rsCreate (int nSym, I32 *symbol)
{ Rskip rs ;
  int nNode = 2 * nSym ; // initial allowance for 2 per symbol
  int y = 2 + nSym + nNode, max = 16 ; while (max < y) max <<= 1 ;
  if (max <= MAX_LINEAR)
    { rs = rsNew (LINEAR, nSym, max) ;
      for (Linear *node = rs.linear+2 ; nSym-- ; node += 2) *(I32*)node = *symbol++ ;
    }
  else
    { nNode = 3*nSym ; // need to allow for some depth
      y = 1 + nSym + nNode ; max = 128 ; while (max < y) max <<= 1 ;
      rs = rsNew (DYNAMIC, nSym, max) ;
      for (Dynamic *node = rs.dynamic + 1 ; nSym-- ; ++node) node->sym = *symbol++ ;
    }
  return rs ;
}

Rskip rsCreateSyng (int nSym, I32 *symbol, U32 *offset)
{ Rskip rs ;
  int nNode = 2 * nSym ; // initial allowance for 2 per symbol
  int max = 16, y = 2 + 4*nSym + nNode ; while (max < y) max <<= 1 ;
  if (max <= MAX_LINEAR)
    { rs = rsNew (LINEAR_SYNG, nSym, max) ;
      for (LinearSyngDir *lsd = linearSyngDir (rs) ; nSym-- ; ++lsd)
	{ lsd->sym = *symbol++ ; lsd->offset = *offset++ ; }
    }
  else
    { nNode = 3*nSym ; // need to allow for some depth
      max = 128 ; y = 1 + nSym + nNode ; while (max < y) max <<= 1 ;
      rs = rsNew (DYNAMIC, nSym, max) ;
      for (Dynamic *node = rs.dynamic + 1 ; nSym-- ; ++node)
	{ node->sym = *symbol++ ; node->offset = *offset++ ; }
    }
  return rs ;
}
  
int rsDestroy (Rskip rs)
{
  if (!rs.linear) return rsError (RS_ERROR_EMPTY, "rsDestroy", "") ;
  if (rs.linear->max) newFree (rs.linear, rs.linear->max, Linear) ;
  else if (rs.dynamic->max) newFree (rs.dynamic, rs.dynamic->max, Dynamic) ;
  else newFree (rs.fixed, rs.fixed->max, Fixed) ;
  return 0 ;
}

int rsNsym (Rskip rs)
{
  if (!rs.linear) return 0 ;
  if (rs.linear->max) return rs.linear->nSym ;
  else if (rs.dynamic->max) return rs.dynamic->nSym ;
  else return rs.fixed->nSym ;
}

int rsSize (Rskip rs, int *linearSize, int *skipSize) // returns max, fills others with byte size
{
  if (!rs.linear) return 0 ;
  if (rs.linear->max)
    { if (linearSize) *linearSize = rs.linear->max * sizeof(LINEAR) ;
      return rs.linear->max ;
    }
  else if (rs.dynamic->max)
    { if (skipSize) *skipSize = rs.dynamic->max * sizeof(DYNAMIC) ;
      return rs.dynamic->max ;
    }
  else
    { if (skipSize) *skipSize = rs.dynamic->max * sizeof(FIXED) ;
      return rs.fixed->max ;
    }
}

/********************* some debugging support  **********************/

static void inline printDynamic (Dynamic *n, int i)
{
  if (n->down) // not at bottom
    { Dynamic *n1 = n ;
      int d = 1, i1 = i ; while (n1->down) { int i2 = n1->down ; n1 += (i2-i1) ; i1 = i2 ; ++d ; }
      printf ("    %3d: iSym %2d depth %2u down %3u up %3u left %3u before %3u right %3u count %3u "
	      "sLeft %3u sBefore %3u sRight %3u sCount %3u\n",
	      i, n->iSym, d, n->down, n->up, n->left, n->before, n->right, n->count, 
	      n->sLeft, n->sBefore, n->sRight, n->sCount) ;
    }
  else  // at bottom
    printf ("    %3d: iSym %2u depth  1 down   0 up %3u left %3u before %3u right %3u count %3u "
	    "sLeft %3u sBefore %3u sRight %3u sCount %3u\n",
	    i, n->iSym, n->up, n->left, n->before, n->right, n->count, 
	    n->sLeft, n->sBefore, n->sRight, n->sCount) ;
}

static void inline printFixed (Fixed *n, int i)
{
  if (n->sSum & FLAG32) // at bottom
    printf ("    %u: iSym %u sum %u right %u sSum %u sRight %u\n",
           i, n->iSym, n->sum, n->right, n->sSum & MASK32, n->sRight) ;
  else
    printf ("    %u:        sum %u right %u sSum %u sRight %u down %u\n",
           i, n->sum, n->right, n->sSum, n->sRight, n->down) ;
}

void rsPrint (Rskip rs)
{
  if (rs.linear->max)
    { U32 i ;
      printf ("Linear max %u nSym %u free %u\n", rs.linear->max, rs.linear->nSym, rs.linear[1].free) ;
      int nSym = rs.linear->nSym ;
      Linear *node = rs.linear + 2 ;
      if (rsType(rs) == LINEAR_SYNG)
	{ LinearSyngDir *lsd = linearSyngDir (rs) ;
	  for (i = 0 ; i < nSym ; ++i, ++lsd)
	    printf ("    DIR %2u sym %d offset %u count %u\n", i, lsd->sym, lsd->offset, lsd->count) ;
	}
      else if (rsType(rs) == LINEAR)
	for (i = 0 ; i < nSym ; ++i, node += 2)
	  printf ("    DIR %2u sym %d\n", i, *(I32*)node) ;
      int sum = 0, sSum[nSym] ; memset (sSum, 0, nSym*sizeof(int)) ;
      int j = 0 ;
      node = rs.linear + rs.linear->max - 1 ;
      for (i = rs.linear->max ; --i > rs.linear[1].free ; --node)
        { int count = (node->count == 255) ? node[-1].bigCount : node->count ;
          sum += count ; sSum[node->iSym] += count ;
          printf ("    %d: i %u iSym %u count %u sum %u sSum %u\n",
		  ++j, i, node->iSym, count, sum, sSum[node->iSym]) ;
          if (node->count == 255) { --i ; --node ; }
        }
    }
  else if (rs.dynamic->max)
    { U32 i ;
      Dynamic *node = rs.dynamic ;
      printf ("Dynamic max %u nSym %u start %u maxDepth %u free %u\n",
	      node->max, node->nSym, node->start, node->maxDepth, node->free) ;
      if (rs.dynamic->nSym < 3000) // 30
	for (i = 1 ; i <= rs.dynamic->nSym ; ++i)
	  printf ("    DIR %3d: sym %8d offset %3u count %4u right %3u depth %2u sum %d\n", i,
		  node[i].sym, node[i].offset, node[i].count, node[i].right, node[i].down, node[i].sum) ;
      if (rs.dynamic->max < 50000) // 500
	for (i = rs.dynamic->free + 1 ; i < rs.dynamic->max ; ++i)
	  printDynamic (&node[i], i) ;
    }
  else
    { U32 i ;
      Fixed *node = rs.fixed ;
      printf ("Fixed max %u nSym %u start %u\n", node->max, node->nSym, node->start) ;
      for (++node, i = 1 ; i <= rs.fixed->nSym ; ++i, ++node)
	printf ("    DIR %3d: sym %8d offset %3u count %u sum %d\n",
		i, node->sym, node->offset, node->count, node->sum) ;
      for ( ; i < rs.fixed->max ; ++i, ++node) printFixed (node, i) ;
    }
}

static bool rsCheck (Rskip rs)
{
  if (!rs.linear) die ("rsCheck called on null pointer") ;
  bool isBad = false ;

  if (rs.linear->max)
    { Linear *node = rs.linear ;
      int nSym = rs.linear->nSym ;
      int kCount[nSym] ; memset (kCount, 0, nSym*sizeof(int)) ;
      for (int i = rs.linear->max ; --i > rs.linear[1].free ; )
	{ int is = node[i].iSym, ic = node[i].count ;
	  if (ic == 255) ic = node[--i].bigCount ;
	  if (is >= nSym)
	    { printf ("--linear node %d has iSym %d not in [0,%d)\n", i, is, nSym) ;
	      isBad = true ;
	    }
	  else
	    kCount[is] += ic ;
	}
    }
  else if (rs.dynamic->max)
    { Dynamic *node = rs.dynamic ;
      int nSym = rs.dynamic->nSym ;

      // check directory entries
      for (int i = 1 ; i <= nSym ; ++i)
	if (node[i].sum)
	  { if (!node[i].right)
	      { printf ("-- dir node %d has null right\n", i) ; isBad = true ; }
	    // verify depth by walking down from the top node
	    int d = 1, j = node[i].right ; while (node[j].down) { j = node[j].down ; ++d ; }
	    if (d != node[i].depth)
	      { printf ("-- dir node %d depth %d != actual %d\n", i, node[i].depth, d) ;
		isBad = true ;
	      }
	  }

      // find the leftmost node
      int i = rs.dynamic->start ;
      while (true)
	if (node[i].left) i = node[i].left ;
	else if (node[i].down) i = node[i].down ;
	else break ;

      // check all bottom nodes
      for ( ; i ; i = node[i].right)
	{ if (node[i].down)
	    { printf ("-- bottom row node %d has down=%u\n", i, node[i].down) ; isBad = true ; }
	  if (node[i].iSym >= nSym)
	    { printf ("-- node %d has sym=%u >= nSym=%d\n", i, node[i].sym, nSym) ; isBad = true ; }
	  // now check column upwards
	  for (int j = i ; j ; j = node[j].up)
	    { if (node[j].right)
		{ if (node[node[j].right].left != j)
		    { printf ("-- node %d right %d left mismatch\n", j, node[j].right) ; isBad = true ; }
		  if (node[node[j].right].before != node[j].count)
		    { printf ("-- node %d count %d != right %d before %d\n", j, node[j].count,
			      node[j].right, node[node[j].right].before) ; isBad = true ; }
		}
	      if (node[j].sRight)
		{ if (node[node[j].sRight].sLeft != j)
		    { printf ("-- node %d sRight %d sLeft mismatch\n", j, node[j].sRight) ;
		      isBad = true ; }
		  if (node[node[j].sRight].sBefore != node[j].sCount)
		    { printf ("-- node %d sCount %d != sRight %d sBefore %d\n", j, node[j].sCount,
			      node[j].sRight, node[node[j].sRight].sBefore) ; isBad = true ; }
		  if (node[node[j].sRight].iSym != node[j].iSym)
		    { printf ("-- node %d sym %d != sRight %d sym %d\n", j, node[j].iSym,
			      node[j].sRight, node[node[j].sRight].iSym) ; isBad = true ; }
		}
	      //if (!node[j].count)
	      //  { printf ("--node %d has no counts\n", j) ; isBad = true ; }
	    }
	}
    }

  if (isBad)
    { printf (" !!!! rsCheck() failed:\n") ;
      rsPrint (rs) ;
      return false ;
    }
  else
    return true ;
}

/************* now the core builders - use I64 inputs so native for ONEcode *************/

static void checkSymRuns (int nSym, int nRun, I64 *iSym, I64 *runLen) // die on fail
{
  for (int i = 0 ; i < nRun ; ++i) // check iSym are in bounds and runLen are positive
    { if (iSym[i] < 0 || iSym[i] >= nSym)
	die ("checkSymRuns: iSym[%d] %lld not in 0..%d", i, iSym[i], nSym) ;
      if (runLen[i] < 0)
	die ("checkSymRuns: runLen[%d] %lld <= 0", i, runLen[i]) ;
    }
}

static int linearSize (int type, int nSym, int nRun, I64 *runLen)
{
  I64 sum = 0, size = 2 + nRun ; // 2 for the header, nRun provisionally for runs
  for (int i = 0 ; i < nRun ; ++i)
    { if (runLen[i] >= 255) { ++size ; if (runLen[i] > MAX_BIG) { size = MAX_LINEAR+1 ; break ; } }
      sum += runLen[i] ;
    }
  if (sum > MAX_BIG) size = MAX_LINEAR+1 ;
  if (type == LINEAR) size += 2*nSym ; // for symbol directory
  else if (type == LINEAR_SYNG) size += 4*nSym ; // for symbol, offset, count directory
  return size ;
}

static bool fillLinearNodes (Rskip rs, int nRun, I64 *iSym, I64 *runLen)
// used both when constructing from external data and when rebuilding/extending during rsAdd()
{
  checkSymRuns (rs.linear->nSym, nRun, iSym, runLen) ;
  Linear *node = rs.linear + rs.linear->max - 1 ;
  Linear *min = rs.linear + 2 ;
  U64 total = 0 ;

  if (rsType(rs) == LINEAR_SYNG) min += 4*rs.linear->nSym ;
  else if (rsType(rs) == LINEAR) min += 2*rs.linear->nSym ;

  for (int i = 0 ; i < nRun ; ++i, --node)
    { if (!runLen[i]) continue ;
      if (node < min) return false ;
      node->iSym = iSym[i] ;
      total += runLen[i] ;
      if (runLen[i] < 255)
	node->count = runLen[i] ;
      else if (runLen[i] <= MAX_BIG)
	{ node->count = 255 ; if (--node < min) return false ;
	  node->bigCount = runLen[i] ;
	}
      else return false ;
    }
  if (total > MAX_BIG)
    return false ;
  
  rs.linear[1].free = node - rs.linear ; // next available free node
      
  return true ;
}

// some utilities for the Fixed and Dynamic graph representations

static inline int randomDepth (void) // generate randomDepth with very simple pseudorandom generator
{ int d = 1 ;
  static U64 rng = 18982392197 ; // a prime - each thread's copy starts with this
  static const int MAX_DEPTH = 32 ; // maximum allowed depth
  rng *= 479 ;  // 479 is prime
  while (((rng >> 24) & 0x7) < 3) { ++d ; rng *= 479 ; }
  // deterministic pseudorandom with probability 3/8 ~= 1/e
  if (d > MAX_DEPTH) d = MAX_DEPTH ;
  return d ;
}

static Rskip buildFixed (U8 type, int nSym, int nRun, I64 *iSym, I64 *runLen)
// makes the rs and fills the nodes - doesn't fill the directory for types FIXED and FIXED_SYNG
{
  checkSymRuns (nSym, nRun, iSym, runLen) ;
  // static int callCount = 0 ; ++callCount ; // for debugging
  // first count some things we will need
  int i, s ;         // generic node counter and symbol index
  int nNode = (type == FIXED_RAW) ? 0 : nSym ; // number of nodes needed - start with the directory
  int iDepth[nRun] ;    // number of nodes/layers for run i
  int maxDepth = 0 ; // maximum of iDepth[]
  int sFirst[nSym] ;    // the lowest run index i at which symbol s occurs
  int sMaxDepth[nSym] ; // max depth for symbol s
  for (s = 0 ; s < nSym ; ++s) { sMaxDepth[s] = 0 ; sFirst[s] = -1 ; }
  for (i = nRun ; i-- ; )  // must be backwards so that sFirst[] works
    { sFirst[iSym[i]] = i ;
      iDepth[i] = randomDepth () ;
      nNode += iDepth[i] ;
      if (iDepth[i] > maxDepth) maxDepth = iDepth[i] ;
      if (iDepth[i] > sMaxDepth[iSym[i]]) sMaxDepth[iSym[i]] = iDepth[i] ;
    }
  for (s = 0 ; s < nSym ; ++s)
    { if (sFirst[s] == -1) die ("unused symbol %d < %d", s, nSym) ; // confirm no unused symbols
      // now ensure that the first column for s has maximal depth of all columns for s
      nNode += sMaxDepth[s] - iDepth[sFirst[s]] ; iDepth[sFirst[s]] = sMaxDepth[s] ;
    }
  // and similarly that the first column of all has maximal depth of all columns
  nNode += maxDepth - iDepth[0] ; iDepth[0] = maxDepth ; sMaxDepth[iSym[0]] = maxDepth ;

  // next make the rs object with the right amount of memory
  Rskip rs = rsNew (type, nSym, 1+nNode) ;
  Fixed *node = rs.fixed ;
  node[0].sum = node[0].sSum = 0 ; // rsNew's type assignment corrupts byte 3 of sum; clear both

  // we will need stack and sStack to set up right and sRight pointers
  int stack[maxDepth] ; memset (stack, 0, maxDepth*sizeof(int)) ; 
  int sStack[nSym][maxDepth] ; memset (sStack, 0, nSym*maxDepth*sizeof(int)) ;
  // initialise stacks to point to node[0], which is the header node
  // we will use this as a generic empty node with sum = sSum = 0
  // this allows the update for the first node for each symbol skiplist to be 
  //   the same as for all others in the lines below starting nd->sum = and nd->sSum =

  int d, r ;             // generic variables for a depth and a run length
  Fixed *nd ;            // generic node pointer
  int freeNode = nNode ; // allocate nodes downwards from nNode
      
  // now loop over the runs - overall nlog_n complexity
  int sOff = (type == FIXED_RAW) ? 1 : 1+nSym ; // start of first columns, leaving 1..nSym for directory
  for (i = 0 ; i < nRun ; ++i)
    { s = iSym[i] ;
      r = runLen[i] ;
      U32 next = (i == sFirst[s]) ? s+sOff : freeNode-- ;
      for (d = iDepth[i] ; d-- ; )
	{ nd = &node[next] ;
	  nd->sum = node[stack[d]].sum ; node[stack[d]].right = next ; stack[d] = next ;
	  nd->sSum = node[sStack[s][d]].sSum ; node[sStack[s][d]].sRight = next ; sStack[s][d] = next ;
          if (d) { next = freeNode-- ; nd->down = next ; }
	}
      nd->sSum |= FLAG32 ; // marker for bottom of the column
      nd->iSym = s ;
      // now go back up the stacks adding on r
      for (d = 0 ; d < maxDepth ; ++d) node[stack[d]].sum += r ;
      for (d = 0 ; d < sMaxDepth[s] ; ++d) node[sStack[s][d]].sSum += r ;
    }

  rs.fixed->sum = rs.fixed->sSum = 0 ;       // need to clear (NB corrupts type byte which overlaps sum)
  rs.linear[1].type = type ;                // restore type after clearing sum
  rs.fixed->max = 1+nNode ; // have to set these now since we used them in stack initialisation
  rs.fixed->start = iSym[0] + sOff ; // actual start
  rs.fixed->nSym = nSym ;
  if (type == FIXED_SYNG || type == FIXED) // record the counts per symbol
    for (s = 0 ; s < nSym ; ++s)
      node[1+s].sum = node[1+s].count = node[sStack[s][sMaxDepth[s]-1]].sSum & MASK32 ;
  return rs ;
}

static Rskip buildDynamic (U8 type, int nSym, int nRun, I64 *iSym, I64 *runLen) // mirrors buildFixed()
{
  checkSymRuns (nSym, nRun, iSym, runLen) ;
  // preliminaries are the same as buildFixed() except we scale max to next power of 2
  int i, s ;         // generic node counter and symbol index
  int nNode = nSym ; // number of nodes needed - allow for directory
  int iDepth[nRun] ; // number of nodes/layers for run i
  int maxDepth = 0 ; // maximum of iDepth[]
  int symMax = -1 ;  // symbol value for max depth
  int sMaxDepth[nSym] ; // max depth for symbol s
  int sFirst[nSym] ; // the lowest run index i at which symbol s has depth sMaxDepth[nsym]
  for (s = 0 ; s < nSym ; ++s) { sMaxDepth[s] = 0 ; sFirst[s] = -1 ; }

  for (i = nRun ; i-- ; ) // plan the columns - must be backwards so that sFirst[] works
    { iDepth[i] = randomDepth () ;
      nNode += iDepth[i] ;
      if (iDepth[i] >= maxDepth) { maxDepth = iDepth[i] ; symMax = iSym[i] ; }
      if (iDepth[i] >= sMaxDepth[iSym[i]]) { sMaxDepth[iSym[i]] = iDepth[i] ; sFirst[iSym[i]] = i ; }
    }

  int size = 128 ; while (size <= nNode) size <<= 1 ; // <= to allow for header
  nNode = size ;

  // next make the rs object with the right amount of memory
  Rskip rs = rsNew (DYNAMIC, nSym, 1+nNode) ;
  Dynamic *node = rs.dynamic ; node->count = 0 ; // clear type == rs.linear[1].type
  
  // preliminaries for building, stack and sStack as in buildFixed()
  int stack[maxDepth] ; memset (stack, 0, maxDepth*sizeof(int)) ;
  int sStack[nSym][maxDepth] ; memset (sStack, 0, nSym*maxDepth*sizeof(int)) ;
  int sTop[nSym] ; memset (sTop, 0, nSym*sizeof(int)) ; // top node index per symbol
  int sum = 0, sSum[nSym] ; memset (sSum, 0, nSym*sizeof(int)) ;

  int d, r ;             // generic variables for a depth and a run length
  int freeNode = nNode ; // allocate nodes downwards from nNode

  // now loop over the runs - overall nlog_n complexity
  for (i = 0 ; i < nRun ; ++i)
    { s = iSym[i] ;
      r = runLen[i] ;
      U32 next = freeNode-- ; // all column nodes from free pool
      if (i == sFirst[s]) sTop[s] = next ; // record top of leftmost max-depth column for s
      for (d = iDepth[i] ; d-- ; ) // NB rely on nodes being cleared to 0 for default values
	{ Dynamic *nd = &node[next] ;
	  nd->left = stack[d] ; node[stack[d]].right = next ; stack[d] = next ;
	  nd->sLeft = sStack[s][d] ; node[sStack[s][d]].sRight = next ; sStack[s][d] = next ;
	  nd->before = nd->left ? node[nd->left].count : sum ;
	  nd->sBefore = nd->sLeft ? node[nd->sLeft].sCount : sSum[s] ;
	  nd->iSym = s ;
          if (d) { node[freeNode].up = next ; next = freeNode-- ; nd->down = next ; }
	}
      sum += r ; sSum[s] += r ;
      // now go back up the stacks adding r to the counts
      for (d = 0 ; d < maxDepth ; ++d) node[stack[d]].count += r ;
      for (d = 0 ; d < sMaxDepth[s] ; ++d) node[sStack[s][d]].sCount += r ;
    }

  rs.dynamic->count = rs.dynamic->sCount = 0 ;      // need to clear
  rs.dynamic->max = 1+nNode ;      // have to set these after building, since space used in building
  rs.dynamic->maxDepth = maxDepth ;
  rs.dynamic->start = sTop[symMax] ;
  rs.dynamic->nSym = nSym ;
  rs.dynamic->free = freeNode ;
  rs.linear[1].type = type ;

  // set up directory entries: node[s+1] for each symbol s
  for (s = 0 ; s < nSym ; ++s)
    { node[s+1].right = sTop[s] ;      // index of top node for this symbol
      node[s+1].depth = sMaxDepth[s] ; // depth of that column
      node[s+1].sum = sSum[s] ;        // total number of this symbol
      node[s+1].sym = s ;              // default symbol id
    }
  return rs ;
}

/************* now the external interfaces to these *************/

Rskip rsBuildFixedSyng (int nSym, I32 *symbol, U32 *offset, int nRun, I64 *iSym, I64 *runLen)
// this is used when building from data stored in a .1gbwt OneFile
{
  Rskip rs ; rs.linear = 0 ;
  int size = linearSize (LINEAR_SYNG, nSym, nRun, runLen) ;
  if (size <= MAX_LINEAR) // build in Linear nodes
    { rs = rsNew (LINEAR_SYNG, nSym, size) ; 
      if (fillLinearNodes (rs, nRun, iSym, runLen))
	{ for (LinearSyngDir *lsd = linearSyngDir (rs) ; nSym-- ; ++lsd)
	    { lsd->sym = *symbol++ ; lsd->offset = *offset++ ; }
	}
      else rs.linear = 0 ;
    }

  if (!rs.linear)
    { rs = buildFixed (FIXED_SYNG, nSym, nRun, iSym, runLen) ;
      for (Fixed *node = rs.fixed + 1 ; nSym-- ; ++node)
	{ node->sym = *symbol++ ; node->offset = *offset++ ; }
    }

  return rs ;
}

Rskip rsBuildDynamicSyng (int nSym, I32 *symbol, U32 *offset, int nRun, I64 *iSym, I64 *runLen)
// very similar to above, except expand to a power of 2 if Linear, and Dynamic not Fixed if not Linear
{
  Rskip rs ;
  int size = linearSize (LINEAR_SYNG, nSym, nRun, runLen) ;
  if (size <= MAX_LINEAR) // build in Linear nodes
    { if (size & (size-1)) // not a power of 2
	size = 1 + (size | (size >> 1) | (size >> 2) | (size >> 4)) ; // make a power of 2
      rs = rsNew (LINEAR_SYNG, nSym, size) ; 
      if (!fillLinearNodes (rs, nRun, iSym, runLen)) die ("screwup filling linear") ;
      for (LinearSyngDir *lsd = linearSyngDir (rs) ; nSym-- ; ++lsd)
	{ lsd->sym = *symbol++ ; lsd->offset = *offset++ ; }
    }
  else
    { rs = buildDynamic (DYNAMIC, nSym, nRun, iSym, runLen) ;
      for (Dynamic *node = rs.dynamic + 1 ; nSym-- ; ++node)
	{ node->sym = *symbol++ ; node->offset = *offset++ ; }
    }

  return rs ;
}

/****************************** code to add ****************************/

static inline int linearSpace (Rskip rs) // returns number of Linear structs available
{
  switch (rsType(rs))
    {
    case LINEAR_RAW:  return rs.linear[1].free - 1 ;
    case LINEAR:      return rs.linear[1].free - (1 + 2*rs.linear->nSym) ;
    case LINEAR_SYNG: return rs.linear[1].free - (1 + 4*rs.linear->nSym) ;
    }
  return 0 ;
}

// first we have a version that converts linear to I64 list then rebuilds

static Rskip rebuildAddLinear (Rskip rs, int k, U32 kSym)
// NB this increments the dynamic directory sum because buildDynamic() does that
{
  static int callCount = 0 ; ++callCount ;
  // if (callCount == DEBUG) printf ("rebuildAddLinear callCount %d\n", callCount) ;
  int     max = rs.linear->max ; // current size
  Linear *node = rs.linear + max - 1 ; // start at end
  I64     iSym[max+2] ; memset (iSym, 0, (max+2)*sizeof(I64)) ;
  I64     runLen[max+2] ; memset (runLen, 0, (max+2)*sizeof(I64)) ;
  int     nRun = 0, sum = 0, nAdded = 0 ; // nRun will be the eventual number of runs
  if (callCount == DEBUG)
    { printf ("rebuildAddLinear callCount %d k %d kSym %d\n", callCount, k, kSym) ;
      rsPrint (rs) ;
    }
  for (int i = max ; --i > rs.linear[1].free ; --node)
    { if (sum == k && node->iSym != kSym) // insert here
	{ iSym[nRun] = kSym ; runLen[nRun++] = 1 ;
	  nAdded = 1 ;
	  sum = -(1<<30) ; // prevent adding again
	}
      iSym[nRun] = node->iSym ;
      if (node->count < 255)
	{ sum += node->count ;
	  if (sum < k) // copy
	    runLen[nRun++] = node->count ;
	  else if (iSym[nRun] == kSym) // increment
	    { runLen[nRun++] = node->count + 1 ;
	      if (node->count == 254) nAdded = 1 ; // will need to extend to bigCount
	      sum = -(1<<30) ; // prevent adding again
	    }
	  else if (sum > k) // insert
	    { runLen[nRun++] = node->count - (sum - k) ;
	      iSym[nRun] = kSym ; runLen[nRun++] = 1 ;
	      iSym[nRun] = node->iSym ; runLen[nRun++] = sum - k ;
	      nAdded = 2 ;
	      sum = -(1<<30) ; // prevent adding again
	    }
	  else // (sum == k && iSym[n] != kSym) - copy
	    runLen[nRun++] = node->count ;
	}
      else
	{ int iS = node->iSym, count = (--node)->bigCount ; --i ;
	  sum += count ;
	  if (sum < k)
	    runLen[nRun++] = count ;
	  else if (iS == kSym) // increment
	    { runLen[nRun++] = count + 1 ;
	      if (count == MAX_BIG) nAdded = MAX_LINEAR+1 ; // force Dynamic
	      sum = -(1<<30) ; // prevent re-adding
	    }
	  else if (sum > k) // insert
	    { runLen[nRun] = k - (sum - count) ;
	      if (runLen[nRun++] < 255) --nAdded ; // this one is smaller than it was
	      iSym[nRun] = kSym ; runLen[nRun++] = 1 ; ++nAdded ;
	      iSym[nRun] = iS ; runLen[nRun] = sum - k ;
	      nAdded += (runLen[nRun++] < 255) ? 1 : 2 ;
	      sum = -(1<<30) ; // prevent re-adding
	    }
	  else // (sum == k && iS != kSym) - copy
	    runLen[nRun++] = count ;
	}
    }
  if (sum == k) // add to end
    { iSym[nRun] = kSym ; runLen[nRun++] = 1 ; ++nAdded ; }

  // now we have built iSym, nRun we can extend the directory which might overwrite up to 4 nodes
  int nSym = rs.linear->nSym ;
  if (kSym == nSym)
    { if (rsType(rs) == LINEAR)
	{ ++rs.linear->nSym ; nAdded += 2 ; } // NB don't change local Sym here, so can use below
      else if (rsType(rs) == LINEAR_SYNG)
	{ ++rs.linear->nSym ; nAdded += 4 ; }
    }

  // now we refill if it still fits, or try to rebuild as linear if possible
  bool isDone = false ;
  if (nAdded <= linearSpace(rs))
    { if (fillLinearNodes (rs, nRun, iSym, runLen))
	isDone = true ;
    }
  else if (rs.linear->max <= MAX_LINEAR/2)
    { Rskip rsOut = rsNew (rsType(rs), rs.linear->nSym, 2*rs.linear->max) ;
      if (rsType(rs) == LINEAR)
	memcpy (rsOut.linear+2, rs.linear+2, 2*nSym*sizeof(Linear)) ;
      else if (rsType(rs) == LINEAR_SYNG)
	memcpy (rsOut.linear+2, rs.linear+2, 4*nSym*sizeof(Linear)) ;
      if (fillLinearNodes (rsOut, nRun, iSym, runLen))
	{ rsDestroy (rs) ;
	  rs = rsOut ;
	  isDone = true ;
	}
      else rsDestroy (rsOut) ;
    }

  if (!isDone) // we have to convert to dynamic
    { Rskip rsOut = buildDynamic (DYNAMIC, rs.linear->nSym, nRun, iSym, runLen) ;
      if (rsType(rs) == LINEAR)
	for (int i = 0 ; i < nSym ; ++i) rsOut.dynamic[1+i].sym = *(I32*)(rs.linear + 2 + 2*i) ;
      else if (rsType(rs) == LINEAR_SYNG)
	{ LinearSyngDir *lsd = linearSyngDir (rs) ;
	  Dynamic *node = rsOut.dynamic + 1 ;
	  for (int i = 0 ; i < nSym ; ++i, ++lsd, ++node)
	    { node->sym = lsd->sym ; node->offset = lsd->offset ; node->count = lsd->count ; }
	} 
      rsDestroy (rs) ;
      rs = rsOut ;
    }

  if (DEBUG && !rsCheck(rs))
    die ("rsCheck failure debugAddLinear: callCount %d k %d kSym %d", callCount, k, kSym) ;
  
  return rs ;
}

static bool doubleLinear (Rskip *rsp) // extends space
{
  // Strategy is to double the size of the Dynamic array, copy the header and directory
  // and move everything else up to the new space at the end, also increasing free accordingly.
  Rskip rs = *rsp ;
  int oldMax = rs.linear->max ;
  int newMax = oldMax * 2 ;
  int nSym = rs.linear->nSym ;

  if (newMax > MAX_LINEAR) return false ;
  
  // allocate new array with doubled size
  Rskip rs2 = rsNew (rsType(rs), nSym, newMax) ;

  Linear *nOld = rs.linear ; // origin of old data - indices start at 1
  Linear *nNew = rs2.linear ; // origin of new data

  // copy header and directory entries
  int nHeader = 2 ;
  if (rsType(rs) == LINEAR_SYNG) nHeader += 4*nSym ;
  else if (rsType(rs) == LINEAR) nHeader += 2*nSym ;
  memcpy (nNew, nOld, nHeader*sizeof(Linear)) ; 
  rs2.linear->max = newMax ; // need to change this

  // calculate offset for moving nodes: they move them from the node above free upwards
  int offset = newMax - oldMax ;
  rs2.linear[1].free = rs.linear[1].free + offset ;
  memcpy (nNew + rs2.linear[1].free + 1, nOld + rs.linear[1].free + 1,
	  (oldMax - rs.linear[1].free - 1)*sizeof(Linear)) ;

  rsDestroy (rs) ;
  *rsp = rs2 ;
  return true ;
}

static Rskip doubleDynamic (Rskip *rsp) // extends space
{
  // Strategy is to double the size of the Dynamic array, leave 0..nSym-1 alone,
  // and move everything else up to the new space at the end, also increasing free
  // accordingly. All pointers need to be updated accordingly.
  Rskip rs = *rsp ;
  int oldMax = rs.dynamic->max ;
  int newMax = oldMax * 2 ;
  int nSym = rs.dynamic->nSym ;
  
  // allocate new array with doubled size
  Rskip rs2 = rsNew (DYNAMIC, nSym, newMax) ;

  // copy header node (node 0)
  Dynamic *nOld = rs.dynamic ; // origin of old data - indices start at 1
  Dynamic *nNew = rs2.dynamic ; // origin of new data

  // copy header and directory entries 1..nSym (leave them in place)
  memcpy (nNew, nOld, (1+nSym)*sizeof(Dynamic)) ; 
  rs2.dynamic->max = newMax ; // need to change this

  // calculate offset for moving nodes: they move them from the node above free upwards
  int offset = newMax - oldMax ;
  rs2.dynamic->free = rs.dynamic->free + offset ;
  memcpy (nNew + rs2.dynamic->free + 1, nOld + rs.dynamic->free + 1,
	  (oldMax - rs.dynamic->free - 1)*sizeof(Dynamic)) ;

  // update directory entries (1..nSym): only .right is a node pointer (.down is depth)
  for (int i = 1 ; i <= nSym ; ++i)
    if (nNew[i].right) nNew[i].right += offset ;
  // update all pointers in the used nodes (free+offset..newMax)
  for (int i = rs2.dynamic->free + 1 ; i < newMax ; ++i)
    { Dynamic *node = &nNew[i] ;
      if (node->right) node->right += offset ;
      if (node->left) node->left += offset ;
      if (node->sRight) node->sRight += offset ;
      if (node->sLeft) node->sLeft += offset ;
      if (node->up) node->up += offset ;
      if (node->down) node->down += offset ;
    }
  if (rs2.dynamic->start) rs2.dynamic->start += offset ; // also need this

  rsDestroy (rs) ;
  *rsp = rs2 ;
  return rs2 ;
}

// some functions we need for the dynamic case

static inline int dynamicSpace (Rskip rs) { return rs.dynamic->free - rs.dynamic->nSym ; }

static inline void addColumn (Rskip rs, U32 kSym, int depth, int k,
			      int iL, int iR, int xL, int xR, int isL, int isR, int xsL, int xsR)
{
  Dynamic *node = rs.dynamic ;

  // determine whether this column should become the new start for kSym
  bool newSymStart = false ;
  if (kSym == rs.dynamic->nSym)
    { newSymStart = true ;
      rs.dynamic->nSym = kSym+1 ;
      node[kSym+1].sym = kSym ; // default
      node[kSym+1].sCount = 0 ;
    }
  else
    { int oldDepth = node[kSym+1].down ; // O(1) lookup from directory
      if (depth == oldDepth) // check if the new column is to the left of the current start
	{ int curTop = node[kSym+1].right ;
	  for (int i = isR ; i ;)
	    if (i == curTop) { oldDepth = 0 ; break ; } // force newSymStart
	    else if (node[i].up) i = node[i].up ;
	    else i = node[i].sRight ;
	}
      if (depth > oldDepth) newSymStart = true ;
    }

  // check if this should also become the overall start (maxDepth)
  if (newSymStart)
    { int maxDepth = rs.dynamic->maxDepth ;
      if (depth == maxDepth)
	{ int curStart = rs.dynamic->start ;
	  for (int i = iR ; i ;)
	    if (i == curStart) { maxDepth = 0 ; break ; }
	    else if (node[i].up) i = node[i].up ;
	    else i = node[i].right ;
	}
      if (depth > maxDepth)
	{ rs.dynamic->maxDepth = depth ;
	  rs.dynamic->start = rs.dynamic->free - (depth-1) ; // this will be the top node
	}
    }

  // now create the column and link it into the existing structure
  int iLast ;
  for (int d = 0 ; d < depth ; ++d)
    { int iNew = rs.dynamic->free-- ; // all nodes from free pool
      Dynamic *nNew = &node[iNew] ;
      nNew->sym = kSym ;
      if (d) { nNew->down = iLast ; node[iLast].up = iNew ; }
      nNew->left = iL ; nNew->before = xL ;
      if (iL) { node[iL].right = iNew ; node[iL].count = xL ; }
      nNew->right = iR ; nNew->count = xR ;
      if (iR) { node[iR].left = iNew ; node[iR].before = xR ; }
      nNew->sLeft = isL ; nNew->sBefore = xsL ;
      if (isL) { node[isL].sRight = iNew ; node[isL].sCount = xsL ; }
      nNew->sRight = isR ; nNew->sCount = xsR ;
      if (isR) { node[isR].sLeft = iNew ; node[isR].sBefore = xsR ; }
      while (iL && !node[iL].up) { xL += node[iL].before ; iL = node[iL].left ; }
      if (iL) iL = node[iL].up ;
      while (iR && !node[iR].up) { xR += node[iR].count ; iR = node[iR].right ; }
      if (iR) iR = node[iR].up ;
      while (isL && !node[isL].up) { xsL += node[isL].sBefore ; isL = node[isL].sLeft ; }
      if (isL) isL = node[isL].up ;
      while (isR && !node[isR].up) { xsR += node[isR].sCount ; isR = node[isR].sRight ; }
      if (isR) isR = node[isR].up ;
      iLast = iNew ;
    }

  // update directory if this is the new start for kSym
  if (newSymStart)
    { node[kSym+1].right = iLast ;
      node[kSym+1].down = depth ;
    }
}

// add kSym at k and return rank

static int addDirect (Rskip *rsp, U32 k, U32 kSym)
{
  Rskip rs = *rsp ;
  int nSym = rsNsym (rs) ;
  static int callCount = 0 ; ++callCount ;

  // linear case first - do increments here and expansions in rebuildAddLinear
  if (rs.linear->max)
    { int     i, sum = 0, sSum = 0 ;
      Linear *node = rs.linear + rs.linear->max - 1 ;
      int     sym, count ;
      for (i = rs.linear->max ; --i > rs.linear[1].free ; --node)
	{ sym = node->iSym ;
	  count = node->count ; if (count == 255) { count = (--node)->bigCount ; --i ; }
	  if (sum == k) break ; // last sym did not match kSym (or at start)
	  sum += count ;
	  if (sym == kSym)
	    { sSum += count ;
	      if (sum >= k) break ; // we can increment this node
	    }
	  else
	    if (sum > k) break ; // we will need to insert (option 3 below)
	}
      if (i == rs.linear[1].free) // at or past end of list
	{ if (sum == k) // doesn't match the last node - add onto the end
	    *rsp = rebuildAddLinear (rs, k, kSym) ; // need to expand/convert
	  else
	    return rsError (RS_ERROR_INDEX_OVERFLOW, "rsAdd", "index %u > sum %d", k, sum) ;
	}
      else if (sym == kSym) // just increment the counter, except for edge cases
	{ if (count < 254) ++node->count ;
	  else if (count > 254 && count < MAX_BIG) ++node->bigCount ;
	  else *rsp = rebuildAddLinear (rs, k, kSym) ; // need to expand/convert
	  sSum -= (sum - k) ;
	  if (sSum < 0) die ("sSum %d < 0 in addDirect callCount %d", sSum, callCount) ;
	}
      else
	*rsp = rebuildAddLinear (rs, k, kSym) ; // need to expand/convert
      if (DEBUG && !rsCheck (*rsp))
	die ("rsCheck failed after linear addDirect call %d k %d kSym %d", callCount, k, kSym) ;
      return sSum ;
    }

  // next the dynamic case
  else if (!rs.linear->max && rs.dynamic->max)
    { Dynamic *node = rs.dynamic ;
      if (callCount == DEBUG)
	{ printf ("addDirect dynamic callCount %d k %u kSym %u\n", callCount, k, kSym) ;
	  rsPrint (rs) ;
	}

#define I_NEED_THIS // though I think I should be able to remove it given changes elsewhere...
#ifdef I_NEED_THIS
      // special case the initial node
      if (!rs.dynamic->start)
	{ if (k) die ("rsAddDirect problem - k %d with empty Dynamic", k) ;
	  int newDepth = randomDepth() ;
	  if (newDepth > dynamicSpace (rs))
	    { doubleDynamic (rsp) ;
	      return addDirect (rsp, k, kSym) ; // recurse
	    }
	  if (kSym == nSym) ++rs.dynamic->nSym ;
	  U32 next = rs.dynamic->start = rs.dynamic[1+kSym].right = rs.dynamic->free-- ;
	  for (int d = newDepth ; d-- ;)
	    { Dynamic *nd = &node[next] ;
	      nd->iSym = kSym ; nd->count = nd->sCount = 1 ;
	      if (d)
		{ node[rs.dynamic->free].up = next ; next = rs.dynamic->free-- ; nd->down = next ; }
	    }
	  ++rs.dynamic[1+kSym].sum ;
	  rs.dynamic->maxDepth = rs.dynamic[1+kSym].depth = newDepth ;
	  return 0 ; // rank has to be 0
	}
#endif

      int iLeft = 0, iRight = rs.dynamic->start, sum = node[iRight].before ;
      // iLeft, iRight will be to the left and right of the new column
      // typically node[iLeft].right = iRight, but at the edges one or the other will be 0
      // we want sum = count up to and not including iRight - we need sum >= k
      if (sum >= k) // go left
	{ while (true) // LOGARITHMIC in r
	    if (sum - node[iRight].before > k ||
		sum - node[iRight].before == k && k) // need to protect against 0 here
	      { sum -= node[iRight].before ; iRight = node[iRight].left ; }
	    else if (node[iRight].down) iRight = node[iRight].down ;
	    else break ;
	  iLeft = node[iRight].left ; // can be 0
	}
      else // go right
	{ iLeft = iRight ;
	  while (true) // LOGARITHMIC in r
	    if (sum + node[iLeft].count < k)
	      { sum = sum + node[iLeft].count ; iLeft = node[iLeft].right ;
		if (iLeft == 0)
		  { rsPrint (rs) ;
		    die ("rsAddDirect callCount %d k %d off end sum = %d length %d",
			 callCount, k, sum, rsLength(rs)) ;
		  }
	      }
	    else if (node[iLeft].down) iLeft = node[iLeft].down ;
	    else break ; // at bottom in correct node
	  sum += node[iLeft].count ;
	  iRight = node[iLeft].right ; // can be 0
	}
      // printf ("rsAddDirect DYNAMIC callCount %d k %d sum %d iLeft %d iRight %d\n", callCount, k, sum, iLeft, iRight) ;
      if (k > sum) die ("k %u > sum %d in rsAddDirect dynamic", k, sum) ;
      if (!iLeft && !k) sum = 0 ; // deal with 0 case protected above

      // now build new columns if necessary, not yet incrementing for the new count
      int isLeft, isRight ; // we will also need these below
      int sSum = 0 ;
      if (sum != k) // within iLeft
	{ if (node[iLeft].iSym == kSym)
	    { isLeft = iLeft ; isRight = node[isLeft].sRight ; sSum = k - sum + node[iLeft].count ;
	      if (sSum < 0)
		die ("-ve sSum %d k %d sum %d iLeft %d node[iLeft].count %u\n",
		     sSum, k, sum, iLeft, node[iLeft].count) ;
	    }
	  else
	    { int newDepth = randomDepth() ;
	      if (newDepth > dynamicSpace (rs))
		{ doubleDynamic (rsp) ; return addDirect (rsp, k, kSym) ; }
	      addColumn (rs, node[iLeft].sym, newDepth, k,
			 iLeft, iRight, k - sum + node[iLeft].count, sum - k,
			 iLeft, node[iLeft].sRight, k - sum + node[iLeft].count,
			 node[iLeft].sCount - (k - sum + node[iLeft].count)) ;
	      if (DEBUG && !rsCheck (rs))
		die ("failed rsCheck after split addColumn: callCount %d", callCount) ;
	      iRight = node[iLeft].right ; // this is the bottom of the new column
	      sum = k ; // by construction - NB will trigger next block via "if (sum == k)"
	      // we don't need to set isLeft and isRight here - we will do that in the next block
	    }
	}

      if (sum == k) // on the boundary between iLeft and iRight - note fall through from above
	{ if (iLeft && node[iLeft].sym == kSym)
	    { isLeft = iLeft ; isRight = node[isLeft].sRight ; sSum = node[iLeft].count ; }
	  else if (iRight && node[iRight].sym == kSym) // can increment iRight
	    { isLeft = iLeft = iRight ; iRight = node[iLeft].right ; isRight = node[isLeft].sRight ; }
	  else // insert new column between iLeft and iRight
	    { int newDepth = randomDepth() ;
	      if (newDepth + (kSym == nSym) > dynamicSpace (rs))
		{ doubleDynamic (rsp) ; return addDirect (rsp, k, kSym) ; }
	      // before we add the column we must find isLeft and isRight
	      if (kSym == rs.dynamic->nSym || !rs.dynamic[1+kSym].sum)
		isLeft = isRight = 0 ;
	      else
		{ for (isLeft = iLeft ; isLeft ; isLeft = node[isLeft].left) // LINEAR in ALPHABET
		    if (node[isLeft].sym == kSym) break ;
		  if (isLeft)
		    isRight = node[isLeft].sRight ;
		  else
		    { for (isRight = iRight ; isRight ; isRight = node[isRight].right)
			if (node[isRight].sym == kSym) break ;
		      if (!isRight)
			{ printf ("SCREWUP: addDirect callCount %d, k %u, kSym %u, nSym %u\n",
				  callCount, k, kSym, rs.dynamic->nSym) ;
			  rsPrint (rs) ;
			  die ("screwup finding isLeft %d, isRight %d - iLeft %d iRight %d",
			       isLeft, isRight, iLeft, iRight) ;
			}
		      // isLeft is 0 in this case
		    }
		}
	      
	      if (callCount == DEBUG)
		printf ("before insert addColumn iLeft %d iRight %d isLeft %d isRight %d\n",
			iLeft, iRight, isLeft, isRight) ;

	      addColumn (rs, kSym, newDepth, k,
			 iLeft, iRight, iLeft ? node[iLeft].count : 0, 0,
			 isLeft, isRight, isLeft ? node[isLeft].sCount : 0, 0) ;
	      	      if (DEBUG && !rsCheck (rs))
	      		die ("failed rsCheck after insert addColumn: callCount %d", callCount) ;
	      // now reset iLeft and isLeft to the bottom of the new column
	      isLeft = iLeft = iLeft ? node[iLeft].right : node[iRight].left ;
	      // iRight = node[iLeft].right ; isRight = node[isLeft].sRight ; // should not change
	    }
	}

      if (callCount == DEBUG)
	printf ("  before increment iLeft %d iRight %d isLeft %d isRight %d\n",
		iLeft, iRight, isLeft, isRight) ;

      // now increment upwards - first iLeft, iRight
      while (iLeft) // increment, update and increment iRight, climb
	{ ++node[iLeft].count ;
	  iRight = node[iLeft].right ; if (iRight) ++node[iRight].before ;
	  while (iLeft && !node[iLeft].up) iLeft = node[iLeft].left ;
	  if (iLeft) iLeft = node[iLeft].up ;
	}
      while (iRight) // climb then increment
	{ while (iRight && !node[iRight].up) iRight = node[iRight].right ;
	  if (iRight) { iRight = node[iRight].up ; if (iRight) ++node[iRight].before ; }
	}
      // then isLeft, isRight exactly the same, but incrementing sSum as moving left
      while (isLeft)
	{ ++node[isLeft].sCount ;
	  isRight = node[isLeft].sRight ; if (isRight) ++node[isRight].sBefore ;
	  while (isLeft && !node[isLeft].up)
	    { sSum += node[isLeft].sBefore ; isLeft = node[isLeft].sLeft ; }
	  if (isLeft) isLeft = node[isLeft].up ;
	}
      while (isRight)
	{ while (isRight && !node[isRight].up) isRight = node[isRight].sRight ;
	  if (isRight) { isRight = node[isRight].up ; if (isRight) ++node[isRight].sBefore ; }
	}
      ++node[kSym+1].sum ; // update directory total
      if (DEBUG && !rsCheck (rs))
	die ("failed rsCheck after increments: callCount %d", callCount) ;
      return sSum ; // this is the only ultimate real return location for Dynamic
    }
  else return rsError (RS_ERROR_LOCKED, "rsAdd", "") ;
}

/****************** directory actions for syng ******************************/

bool  rsDirSyng (Rskip rs, int iSym, I64 *sym, I64 *offset, I64 *count)
{
  if (!rs.linear || iSym < 0 || iSym >= rsNsym(rs)) return false ;
  if (rs.linear->max)
    { if (rsType(rs) != LINEAR_SYNG) return false ;
      LinearSyngDir *lsd = linearSyngDir(rs) + iSym ;
      if (sym) *sym = lsd->sym ;
      if (offset) *offset = lsd->offset ;
      if (count) *count = lsd->count ;
    }
  else if (rs.dynamic->max)
    { Dynamic *node = rs.dynamic + 1 + iSym ;
      if (sym) *sym = node->sym ;
      if (offset) *offset = node->offset ;
      if (count) *count = node->count ;
    }
  else
    { Fixed *node = rs.fixed + 1 + iSym ;
      if (sym) *sym = node->sym ;
      if (offset) *offset = node->offset ;
      if (count) *count = node->count ;
    }
  return true ;
}

int rsDirRankSyng (Rskip rs, I32 symbol, U32 offset)
{
  int i, sum = 0, nSym = rsNsym(rs) ;
  if (rs.linear->max)
    { if (rsType(rs) != LINEAR_SYNG) die ("non-syng type %d in rsDirRankSyng", rsType(rs)) ;
      LinearSyngDir *lsd = linearSyngDir(rs) ;
      for (i = 0 ; i < nSym ; ++i)
	if (lsd[i].sym == symbol && lsd[i].offset == offset) // found it
	  return sum ;
	else
	  sum += lsd[i].count ;
      die ("failed to fine symbol %d offset %u in linear rskip", symbol, offset) ;
    }
  else if (rs.dynamic->max) // dynamic
    { Dynamic *node = rs.dynamic + 1 ;
      for (i = 0 ; i < nSym ; ++i)
	if (node[i].sym == symbol && node[i].offset == offset) // found it
	  return sum ;
	else
	  sum += node[i].count ;
      die ("failed to fine symbol %d offset %u in linear rskip", symbol, offset) ;
    }
  else 
    die ("non-syng type %d in rsDirRankSyng", rsType(rs)) ;

  return -1 ;
}

int rsDirAddSyng (Rskip *rsp, I32 symbol, U32 offset)  // increment count
{
  Rskip rs = *rsp ;
  int i, sum = 0, nSym = rsNsym(rs) ;
  if (rs.linear->max)
    { if (rsType(rs) != LINEAR_SYNG) die ("non-syng type %d in rsDirAddSyng", rsType(rs)) ;
      LinearSyngDir *lsd = linearSyngDir(rs) ;
      for (i = 0 ; i < nSym ; ++i)
	if (lsd[i].sym == symbol && lsd[i].offset == offset) // found it
	  { if (lsd[i].count == MAX_BIG) // convert to dynamic and recurse
	      { int size = rs.linear->max - rs.linear[1].free ;
		I64 nRun, *iSym = new(size, I64), *runLen = new(size, I64) ;
		rsLinearise (rs, &nRun, iSym, runLen) ;
		*rsp = buildDynamic (DYNAMIC, nSym, nRun, iSym, runLen) ; // forces DYNAMIC
		Dynamic *node = rsp->dynamic + 1 ;
		for (i = 0 ; i < nSym ; ++i)
		  { node[i].sym = lsd[i].sym ; node[i].offset = lsd[i].offset ;
		    node[i].count = lsd[i].count ;
		  }
		rsDestroy (rs) ;
		newFree (iSym, size, I64) ; newFree (runLen, size, I64) ;
		return rsDirAddSyng (rsp, symbol, offset) ;
	      }
	    ++lsd[i].count ;
	    return sum ;
	  }
	else
	  sum += lsd[i].count ;
      
      // only get here if we didn't find it - need to increment nSym and create it
      if (linearSpace(rs)*sizeof(Linear) >= sizeof(LinearSyngDir)) // just add
	{ lsd[nSym].sym = symbol ; lsd[nSym].offset = offset ; lsd[nSym].count = 1 ;
	  ++rs.linear->nSym ;
	}
      else if (doubleLinear (rsp)) // was able to double and stay linear - add to new directory
	{ lsd = linearSyngDir(*rsp) ;
	  lsd[nSym].sym = symbol ; lsd[nSym].offset = offset ; lsd[nSym].count = 1 ;
	  ++(*rsp).linear->nSym ;
	}
      else // must convert to dynamic
	{ int size = rs.linear->max - rs.linear[1].free ;
	  I64 nRun, *iSym = new(size, I64), *runLen = new(size, I64) ;
	  rsLinearise (rs, &nRun, iSym, runLen) ;
	  I32 *aSym = new(nSym+1, I32) ; U32 *aOff = new(nSym+1, U32), *aCnt = new(nSym+1, U32) ;
	  for (i = 0 ; i < nSym ; ++i)
	    { aSym[i] = lsd[i].sym ; aOff[i] = lsd[i].offset ; aCnt[i] = lsd[i].count ; }
	  aSym[nSym] = symbol ; aOff[nSym] = offset ; aCnt[nSym] = 1 ;
	  rsDestroy (rs) ;
	  rs = *rsp = rsBuildDynamicSyng (nSym+1, aSym, aOff, nRun, iSym, runLen) ;
	  for (i = 0 ; i <= nSym ; ++i) rs.dynamic[1+i].count = aCnt[i] ;
	  newFree (iSym, size, I64) ; newFree (runLen, size, I64) ;
	  newFree (aSym, nSym+1, I32) ; newFree (aOff, nSym+1, U32) ; newFree (aCnt, nSym+1, U32) ;
	}
    }
  else if (rs.dynamic->max) // dynamic
    { Dynamic *node = rs.dynamic + 1 ;
      for (i = 0 ; i < nSym ; ++i)
	if (node[i].sym == symbol && node[i].offset == offset) // found it
	  { ++node[i].count ; return sum ; }
	else
	  sum += node[i].count ;

      // again, only get here if we didn't find it - much simpler in this case
      if (dynamicSpace (rs) == 0) rs = doubleDynamic (rsp) ;
      node = rs.dynamic + 1 ;
      node[nSym].sym = symbol ; node[nSym].offset = offset ; node[nSym].count = 1 ;
      ++rs.dynamic->nSym ;
    }
  else 
    die ("non-syng type %d in rsDirAddSyng", rsType(rs)) ;

  return sum ;
}

U32 rsDirSum (Rskip rs)
{
  U32 sum = 0 ;
  if (rsType(rs) == LINEAR_SYNG)
    { LinearSyngDir *lsd = linearSyngDir(rs) ;
      for (int i = 0 ; i < rs.linear->nSym ; ++i) sum += lsd[i].count ;
    }
  else if (rsType(rs) == DYNAMIC)
    { Dynamic *node = rs.dynamic + 1 ;
      for (int i = 0 ; i < rs.dynamic->nSym ; ++i) sum += node[i].count ;
    }
  else die ("unsupported type %d in rsDirSetCount0", rsType(rs)) ;
  return sum ;
}

void rsDirSetCount (Rskip rs, U32 iSym, U32 count)
{
  if (rsType(rs) == LINEAR_SYNG)  linearSyngDir(rs)[iSym].count = count ;
  else if (rsType(rs) == DYNAMIC) rs.dynamic[1+iSym].count = count ;
  else die ("unsupported type %d in rsDirSetCount0", rsType(rs)) ;
}

/****************** public interface size, length, count, rank, find, add *********************/

// first a couple of routines to look up the symbol, and in the case of syng also the offset

static inline int kSymFind (Rskip rs, I32 symbol)
{
  int i, nSym = rsNsym(rs) ;
  switch (rsType(rs))
    {
    case LINEAR_RAW: case FIXED_RAW:
      if (symbol >= 0 && symbol <= rsNsym (rs)) return symbol ;
      else die ("bad raw symbol %d not in 1..%d", symbol, rsNsym(rs)) ;
    case LINEAR:
      for (i = 0 ; i < nSym ; ++i) if (rs.linear[2+i].bigCount == symbol) return i ;
      return nSym ;
    case LINEAR_SYNG:
      { LinearSyngDir *lsd = linearSyngDir (rs) ;
        for (i = 0 ; i < nSym ; ++i, ++lsd) if (lsd->sym == symbol) return i ;
        return nSym ;
      }
    case DYNAMIC:
      for (i = 0 ; i < nSym ; ++i) if (rs.dynamic[1+i].sym == symbol) return i ;
      return nSym ;
    case FIXED: case FIXED_SYNG:
      for (i = 0 ; i < nSym ; ++i) if (rs.fixed[1+i].sym == symbol) return i ;
      return nSym ;
    default: die ("kSymFind called on bad rs type %d", rsType(rs)) ;
    }
  return 0 ;
}

static inline int kSymSyng (Rskip rs, I32 symbol, U32 offset)
{
  int i, nSym = rsNsym(rs) ;
  switch (rsType(rs))
    {
    case LINEAR_SYNG:
      { LinearSyngDir *lsd = linearSyngDir (rs) ;
        for (i = 0 ; i < nSym ; ++i, ++lsd)
	  if (lsd->sym == symbol && lsd->offset == offset) return i ;
        return nSym ;
      }
    case DYNAMIC:
      { Dynamic *node = rs.dynamic + 1 ;
        for (i = 0 ; i < nSym ; ++i, ++node)
	  if (node->sym == symbol && node->offset == offset) return i ;
        return nSym ;
      }
    case FIXED_SYNG:
      { Fixed *node = rs.fixed + 1 ;
        for (i = 0 ; i < nSym ; ++i, ++node)
	  if (node->sym == symbol && node->offset == offset) return i ;
        return nSym ;
      }
    default: die ("kSymSyng called on bad rs type %d", rsType(rs)) ;
    }
  return 0 ;
}

int rsLength (Rskip rs)
{
  if (!rs.linear) return rsError (RS_ERROR_EMPTY, "rsLength", "") ;
  if (rs.linear->max) // linear
    { int sum = 0 ;
      Linear *node = rs.linear + rs.linear->max - 1 ;
      for (int i = rs.linear->max ; --i > rs.linear[1].free ; --node)
	if (node->count < 255) sum += node->count ;
	else { sum += (--node)->bigCount ; --i ; }
      return sum ;
    }
  else if (rs.dynamic->max) // dynamic
    { Dynamic *node = rs.dynamic ;
      int i = rs.dynamic->start, sum = node[i].before ;
      for ( ; i ; i = node[i].right) sum += node[i].count ;
      return sum ;
    }
  else // fixed
    { Fixed *node = rs.fixed ;
      U32 i ;
      for (i = rs.fixed->start ; node[i].right ; i = node[i].right) { ; }
      return node[i].sum ;
    }
}

static int count (Rskip rs, int kSym)	// how many of symbol
{
  // static int callCount = 0 ; ++callCount ;
  if (kSym == rsNsym(rs)) return 0 ; // new symbol
  if (rs.linear->max) // linear
    { int i, sSum = 0 ;
      Linear *node = rs.linear + rs.linear->max - 1 ;
      for (i = rs.linear->max ; --i > rs.linear[1].free ; --node)
	if (node->count < 255) { if (node->iSym == kSym) sSum += node->count ; }
	else { --i ; if ((node--)->iSym == kSym) sSum += node->bigCount ; }
      return sSum ;
    }
  else if (rs.dynamic->max) // dynamic - O(1) lookup from directory
    return rs.dynamic[kSym+1].sum ;
  else // fixed
    { if (rsType(rs) == FIXED_SYNG || rsType(rs) == FIXED) return rs.fixed[kSym+1].sum ;
      U32 i ;
      Fixed *node = rs.fixed ;
      for (i = kSym+1 ; node[i].sRight ; i = node[i].sRight) { ; }
      return node[i].sSum & MASK32 ;
    }
}

int rsCount (Rskip rs, I32 symbol)
{ int kSym = kSymFind (rs, symbol) ; return count (rs, kSym) ; }

int rsCountSyng (Rskip rs, I32 symbol, U32 offset)
{ int kSym = kSymSyng (rs, symbol, offset) ; return count (rs, kSym) ; }

static int rank (Rskip rs, U32 k, int kSym)  // how many of symbol up to (not including) k
{
  // static int callCount = 0 ; ++callCount ;
  if (kSym == rsNsym(rs)) return 0 ; // new symbol
  if (rs.linear->max) // linear
    { int     i, sum = 0, sSum = 0 ;
      Linear *node = rs.linear + rs.linear->max - 1 ;
      for (i = rs.linear->max ; --i > rs.linear[1].free ; --node)
	{ int sym = node->iSym, count = node->count ;
	  if (count == 255) { count = (--node)->bigCount ; --i ; }
	  sum += count ;
	  if (sym == kSym) sSum += count ;
	  if (sum > k) return (sym == kSym) ? sSum - (sum - k) : sSum ;
	}
      if (sum == k) return sSum ; // for k == length condition, which is legal
    }
  else if (rs.dynamic->max) // dynamic - must work in the main list
    { if (!k) return 0 ; // by definition - protects against case noted below
      Dynamic *node = rs.dynamic ;
      U32 i = rs.dynamic->start, sum = node[i].before ;
      if (sum >= k) // go left
      // 260308 - maybe I could change >= to > in the line above and 2 lines down, and be cleaner?
	{ while (true) // LOGARITHMIC in r
	    if (sum - node[i].before >= k) { sum -= node[i].before ; i = node[i].left ; }
	    else if (node[i].down) i = node[i].down ;
	    else break ;
	  sum -= node[i].before ; i = node[i].left ; // NB node[i] defined because of check above
	}
      else // go right
	while (true) // LOGARITHMIC in r
	  if (sum + node[i].count < k)
	    { sum = sum + node[i].count ; i = node[i].right ;
	      if (i == 0) return count(rs, kSym) ;
	    }
	  else if (node[i].down) i = node[i].down ;
	  else break ; // at bottom in correct node

      U32 sSum = 0 ;
      if (node[i].iSym == kSym) sSum = k - sum ;
      else
	{ for (i = node[i].left ; i && node[i].iSym != kSym ; i = node[i].left) { ; }
	  if (i) sSum = node[i].count ;
	}
      while (i)
	if (node[i].up) i = node[i].up ;
	else { sSum += node[i].sBefore ; i = node[i].sLeft ; }
      return sSum ;
    }
  else // fixed - can follow the symbol sub-list, with care
    { Fixed *node = rs.fixed ;
      U32 i = 1+kSym ; if (rsType(rs) == FIXED_SYNG || rsType(rs) == FIXED) i += rsNsym(rs) ;
      int sSumLast = 0 ;
      while (true) // LOGARITHMIC
        if (node[i].sum < k)
	  { U32 j = node[i].sRight ;
            if (!j) 
              { if (!(node[i].sSum & FLAG32)) i = node[i].down ; 
                else return node[i].sSum & MASK32 ; 
              }
            else if (node[j].sum < k || (node[i].sSum & FLAG32)) 
              { sSumLast = node[i].sSum & MASK32 ; 
                i = j ; 
	      }
            else 
              i = node[i].down ;
          }
	else if (!(node[i].sSum & FLAG32))
	  i = node[i].down ;
	else
	  break ; // at bottom in correct node
      int val = (node[i].sSum & MASK32) - (node[i].sum - k) ; 
      if (val < sSumLast) // i.e. k is before this cell
	return sSumLast ;
      else
	return val ;
    }
  die ("shouldn't get here in rsRank") ; return -1 ;
}

int rsRank (Rskip rs, U32 k, I32 symbol)
{ int kSym = kSymFind (rs, symbol) ; return rank (rs, k, kSym) ; }

int rsRankSyng (Rskip rs, U32 k, I32 symbol, U32 offset)
{ int kSym = kSymSyng (rs, symbol, offset) ; return rank (rs, k, kSym) ; }

static int find (Rskip rs, U32 k, int *kSym) // quite like rank()
{
  // static int callCount = 0 ; ++callCount ;
  if (rs.linear->max) // linear
    { U32     i, sum = 0, sSum = 0 ;
      Linear *node = rs.linear + rs.linear->max - 1 ;
      int     sSums[rs.linear->nSym] ; memset (sSums, 0, rs.linear->nSym * sizeof(int)) ;
      for (i = rs.linear->max ; --i > rs.linear[1].free ; --node)
	{ int sym = node->iSym, count = node->count ;
	  if (count == 255) { count = (--node)->bigCount ; --i ; }
          sum += count ;
	  sSums[sym] += count ;
	  if (sum > k)
	    { if (kSym) *kSym = sym ;
	      return sSums[sym] - (sum - k) ;
	    }
	}
    }
  else if (rs.dynamic->max) // dynamic
    { Dynamic *node = rs.dynamic ;
      U32 i = rs.dynamic->start, sum = node[i].before ;
      if (sum > k) // go left
	{ while (true) // LOGARITHMIC in r
	    if (sum - node[i].before > k) { sum -= node[i].before ; i = node[i].left ; }
	    else if (node[i].down) i = node[i].down ;
	    else break ;
	  sum -= node[i].before ; i = node[i].left ; // NB i can't be 0
	}
      else // go right
	while (true) // LOGARITHMIC in r
	  if (sum + node[i].count <= k)
	    { sum = sum + node[i].count ; i = node[i].right ;
	      if (!i) die ("input k = %u > max in rsFind() dynamic", k) ;
	    }
	  else if (node[i].down) i = node[i].down ;
	  else break ; // at bottom in correct node

      if (kSym) *kSym = node[i].iSym ;
      U32 sSum = k - sum ;
      while (i)
	if (node[i].up) i = node[i].up ;
	else { sSum += node[i].sBefore ; i = node[i].sLeft ; }
      return sSum ;
    }
  else // fixed
    { Fixed *node = rs.fixed ;
      U32 i = rs.fixed->start ;
      while (true) // LOGARITHMIC
	if (node[i].sum <= k)
          { i = node[i].right ;
	    if (!i) die ("input k = %u > max in rsFind() fixed, i %u", k, i) ;
          }
	else if (!(node[i].sSum & FLAG32)) i = node[i].down ;
	else break ; // at bottom in correct node
      if (kSym) *kSym = node[i].iSym ;
      return (node[i].sSum & MASK32) - (node[i].sum - k) ;
    }
  die ("shouldn't get here in rsFind") ; return -1 ;
}

int rsFind (Rskip rs, U32 k, I32 *symbol)
{ int kSym, rank = find (rs, k, &kSym) ;
  switch (rsType(rs))
    {
    case LINEAR:  *symbol = *(I32*)(&rs.linear[2+2*kSym]) ; break ;
    case FIXED:   *symbol = rs.fixed[1+kSym].sym ; break ;
    case DYNAMIC: *symbol = rs.dynamic[1+kSym].sym ; break ;
    case LINEAR_RAW: case FIXED_RAW: *symbol = kSym ; break ;
    }
  return rank ;
}

int rsFindSyng (Rskip rs, U32 k, I32 *symbol, U32 *offset)
{ int kSym, rank = find (rs, k, &kSym) ;
  switch (rsType(rs))
    {
    case LINEAR_SYNG:
      { LinearSyngDir *lsd = linearSyngDir(rs) ; *symbol = lsd[kSym].sym ; *offset = lsd[kSym].offset ; }
      break ;
    case FIXED_SYNG:
      *symbol = rs.fixed[1+kSym].sym ; *offset = rs.fixed[1+kSym].offset ;
      break ;
    case DYNAMIC:
      *symbol = rs.dynamic[1+kSym].sym ; *offset = rs.dynamic[1+kSym].offset ;
      break ;
    }
  return rank ;
}

int rsAdd (Rskip *rsp, U32 k, I32 symbol)
{
  int kSym = kSymFind (*rsp, symbol) ;
  if (kSym < rsNsym(*rsp)) return addDirect (rsp, k, kSym) ; // NB must call addDirect AFTER test
  // else this was a new symbol - must add it to the directory
  addDirect (rsp, k, kSym) ;
  if (rsType(*rsp) == LINEAR) *(I32*)(&rsp->linear[2+2*kSym]) = symbol ;
  else if (rsType(*rsp) == DYNAMIC) rsp->dynamic[1+kSym].sym = symbol ;
  return 0 ; // there can not have been anything before, so rank was 0
}

int rsAddSyng (Rskip *rsp, U32 k, I32 symbol, U32 offset)
{
  int kSym = kSymSyng (*rsp, symbol, offset) ;
  if (kSym < rsNsym(*rsp)) return addDirect (rsp, k, kSym) ; // NB must call addDirect AFTER test
  // else this was a new symbol - must add it to the directory
  addDirect (rsp, k, kSym) ;
  if (rsType(*rsp) == LINEAR_SYNG)
    { LinearSyngDir *lsd = linearSyngDir(*rsp) ;
      lsd[kSym].sym = symbol ; lsd[kSym].offset = offset ; lsd[kSym].count = 0 ;
    }
  else if (rsType(*rsp) == DYNAMIC)
    { Dynamic *node = rsp->dynamic + 1 ;
      node[kSym].sym = symbol ; node[kSym].offset = offset ; node[kSym].count = 0 ;
    }
  return 0 ; // there can not have been anything before, so rank was 0
}

/************ linearisation - reversal of rsBuild*Syng()  ****************/

bool rsLinearise  (Rskip rs, I64 *nRun, I64 *iSym, I64 *runLen)
{
  if (!rs.linear) return false ;
  I64 n = 0 ;
  if (rs.linear->max)
    { Linear *node = rs.linear + rs.linear->max - 1 ;
      for (int i = rs.linear->max ; --i > rs.linear[1].free ; --node)
	{ if (iSym) *iSym++ = node->iSym ;
	  if (node->count < 255) { if (runLen) *runLen++ = node->count ; }
	  else { --node ; --i ; if (runLen) *runLen++ = node->bigCount ; }
	  ++n ;
	}
    }
  else if (rs.dynamic->max)
    { Dynamic *node = rs.dynamic ;
      
      // first find the leftmost node
      int i = rs.dynamic->start ;
      while (true)
	if (node[i].left) i = node[i].left ;
	else if (node[i].down) i = node[i].down ;
	else break ;
      // then fill
      for ( ; i ; i = node[i].right)
	{ if (iSym) *iSym++ = node[i].iSym ;
	  if (runLen) *runLen++ = node[i].count ;
	  ++n ;
	}
    }
  else
    { Fixed *node = rs.fixed ;
      
      // first find go to bottom of the start column
      int i = rs.dynamic->start ; while (node[i].down) i = node[i].down ;
      // then fill
      for ( ; i ; i = node[i].right)
	{ if (iSym) *iSym++ = node[i].iSym ;
	  if (runLen) *runLen++ = node[i].count ;
	  ++n ;
	}
    }
  if (nRun) *nRun = n ;
  return true ;
}

/***************************** test package ***************************/

#ifdef TEST

// compile with: cc -g -DTEST -o rtest rskip.c utils.c ONElib.c -lz  

#include "ONElib.h"

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  --argc ; ++argv ;
  bool isDynamic = false, isBuild = false, isTime = false, isBlank = false ;
  bool isLength = false, isCount = false, isRank = false, isFind = false, isLinear = false ;
  while (argc && *argv[0] == '-')
    if (!strcmp (argv[0], "-dynamic")) { isDynamic = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-build")) { isBuild = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-length")) { isLength = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-count")) { isCount = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-rank")) { isRank = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-find")) { isFind = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-time")) { isTime = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-blank")) { isBlank = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-linear")) { isLinear = true ; --argc ; ++argv ; }
    else die ("unknown option %s", argv[0]) ;
  if (argc != 1) die ("%s XX.1gbwt", argv[-1]) ;

  if (isBlank) isLength = isCount = isRank = isFind = isTime = false ; // meaningless in this case

  OneFile *of = oneFileOpenRead (argv[0], 0, "gbwt", 1) ;
  if (!of || strcmp (of->subType, "gbwt")) die ("can't open .1gbwt file %s", argv[0]) ;
  I64 nc, nC, maxc, maxC ;
  if (oneStats (of, 'c', &nc, &maxc, 0) && oneStats (of, 'C', &nC, &maxC, 0))
    { nC += nc ;
      if (maxc > maxC) maxC = maxc ;
    }
  else die ("failed to count the number of GBWT entries") ;
  if (!oneGoto (of, 'V', 1)) die ("can't locate to a graph node 'V' line in %s", of->fileName) ;

  if (RAND_MAX != ((long long)1<<31)-1) die ("RAND_MAX too small %'lld\n", (long long) RAND_MAX) ;

  U64 *runHist = new0(32, U64), *lenHist = new0(32, U64) ;
  I64 *B, *C ;
  U64 nLinear = 0, totSymLinear = 0, totMaxLinear = 0, totRunsLinear = 0, totLenLinear = 0, maxLenLinear = 0, nBigLinear = 0 ;
  U64 nDynamic = 0, totSymDynamic = 0, totMaxDynamic = 0, totRunsDynamic = 0, totLenDynamic = 0, maxLenDynamic = 0 ;
  U64 nFixed = 0, totSymFixed = 0, totMaxFixed = 0, totRunsFixed = 0, totLenFixed = 0, maxLenFixed = 0 ;
  U64 totCount0 = 0, totRank0 = 0, totFind = 0 ;
  int i, n, nVertex = 0, callCount = 0 ;
  int *lens = new (nC, int) ;
  int *D = new (maxC, int) ;
  Rskip *aRS = new (nC, Rskip) ;
  U32 buildListSize = 1<<20, *buildList = new (buildListSize, U32) ;
  U32 *buildCount = new (buildListSize, U32) ;
  static const int MAX_SYM = 1<<16 ;
  I32 *symbols = new (MAX_SYM, I32) ; for (i = 0 ; i < MAX_SYM ; ++i) symbols[i] = i ;
  U32 *zeros = new0 (MAX_SYM, U32) ;
  while (oneReadLine (of))
    switch (of->lineType)
      {
      case 'V': ++nVertex ; break ;
      case 'B': case 'b': B = oneIntList(of) ; break ;
      case 'C': case 'c':
        C = oneIntList(of) ;

	// debug
	++callCount ;
	// printf ("V %d count %d runs %d\n", nVertex, callCount, (int)oneLen(of)) ;
	if (callCount == DEBUG)
	  { printf ("  B") ; for (i = 0 ; i < oneLen(of) ; ++i) printf (" %d",(int)B[i]) ;
	    printf ("\n  C") ; for (i = 0 ; i < oneLen(of) ; ++i) printf (" %d",(int)C[i]) ;
	    putchar ('\n') ;
	  }
	
	int nRuns = oneLen(of) ;
	int kRuns = 0 ; n = nRuns ; while (n >>= 1) ++kRuns ;

	int len = 0 ;
	for (i = 0 ; i < nRuns ; ++i) len += C[i] ;
	lens[callCount-1] = len ;
	int kLen = 0 ; n = len ; while (n >>= 1) ++kLen ;

	int nSym = 0 ; for (i = 0 ; i < nRuns ; ++i) if (B[i] > nSym) nSym = B[i] ; ++nSym ;
	if (nSym > MAX_SYM) die ("nSym %d > MAX_SYM %d", nSym, MAX_SYM) ;
	
	// now make the rs
        Rskip rs ; rs.linear = 0 ;
	if (isBuild)
	  { // strategy is to make an array of symbols of length len
	    // then we need to keep track of which symbols we have used, and map them to
	    // incremental values
	    //if (nRuns > 10000)
	    //  printf ("callCount %d nRuns %d len %d nSym %d\n", callCount, nRuns, len, nSym) ;
	    if (len*2 > buildListSize) // need the 2* for buildCount
	      { newFree (buildList, buildListSize, U32) ; newFree (buildCount, buildListSize, U32) ;
		buildListSize = len*2 ; buildList = new(buildListSize, U32) ;
		buildCount = new (buildListSize, U32) ;
	      }
	    for (U32 *u = buildList, i = 0 ; i < nRuns ; ++i)
	      for (int j = 0 ; j < C[i] ; ++j) *u++ = B[i]+1 ; // NB set > 0
	    memset (buildCount, 0, (2*len)*sizeof(U32)) ; // how many in [x & (x-1), x)
	    int nMap = 0 ;
	    rs = rsCreateSyng (nSym, symbols, zeros) ;
	    for (i = 0 ; i < len ; ++i) // sample 
	      {	int j = rand() % len ; while (j < len && !buildList[j]) ++j ;
		if (j == len) for (j = 0 ; j < len && !buildList[j] ; ) ++j ;
		if (j == len) die ("j search failed: i %d len %d", i, len) ;
		// the total of this j search above will be O(NlogN) - OK
		int b = buildList[j] ;
		buildList[j] = 0 ; // flag as no longer available
		U32 k = 0, x = 1 ;
		for (U32 j0 = j, x = 1 ; x < len ; x <<= 1) // logarithmic - very cool
		  if (j & x) { k += buildCount[j0] ; j0 &= ~x ; }
		  else ++buildCount[j0+x] ;
                if (!isBlank && rsAddSyng (&rs, k, b-1, 0) < 0) // an error
		  die ("rsAddDirect failed") ;
		if (isLinear && rsType(rs) == DYNAMIC) break ;
	      }
	  }
	else if (!isBlank)
	  { if (isDynamic) rs = rsBuildDynamicSyng (nSym, symbols, zeros, nRuns, B, C) ;
	    else rs = rsBuildFixedSyng (nSym, symbols, zeros, nRuns, B, C) ;
	  }
	
	if (!isBlank)
          { if (callCount == DEBUG) rsPrint (rs) ;
	    if (rs.linear->max)
	      { ++nLinear ;
		totSymLinear += rs.linear->nSym ;
		totMaxLinear += rs.linear->max ;
		totRunsLinear += nRuns ;
		totLenLinear += len ;
		if (len > maxLenLinear) maxLenLinear = len ;
		for (i = 0 ; i < nRuns ; ++i) if (C[i] > 256) ++nBigLinear ;
	      }
	    else if (rs.dynamic->max)
	      { ++nDynamic ;
		totSymDynamic += rs.dynamic->nSym ;
		totMaxDynamic += rs.dynamic->max ;
		totRunsDynamic += nRuns ;
		totLenDynamic += len ;
		if (len > maxLenDynamic) maxLenDynamic = len ;
	      }
	    else
	      { ++nFixed ;
		totSymFixed += rs.fixed->nSym ;
		totMaxFixed += rs.fixed->max ;
		totRunsFixed += nRuns ;
		totLenFixed += len ;
		if (len > maxLenDynamic) maxLenDynamic = len ;
	      }
	  }
	
	// Verify correctness by computing expected values from B and C arrays
	if (isLength)
	  { int actual = rsLength (rs) ;
	    if (actual != len)
	      { rsPrint (rs) ;
		die ("V %d callCount %d length() mismatch: expected %d, got %d",
		     nVertex, callCount, len, actual) ;
	      }
	  }
	if (isCount)
	  { int expected = 0 ;
	    for (i = 0 ; i < nRuns ; ++i) if (B[i] == 0) expected += C[i] ;
	    int actual = rsCountSyng (rs, 0, 0) ;
	    if (actual != expected)
	      { rsPrint (rs) ;
		die ("V %d callCount %d count(0) mismatch: expected %d, got %d",
		     nVertex, callCount, expected, actual) ;
	      }
	    totCount0 += actual ;
	  }
	if (isRank)
	  { int expected = 0, sum = 0 ;
	    int target = len/2 ;
	    for (i = 0 ; i < nRuns ; ++i)
	      { if (B[i] == 0) expected += C[i] ;
		sum += C[i] ;	
		if (sum >= target)
		  { if (B[i] == 0) expected -= (sum - target) ;
		    break ;
		  }
	      }
	    int actual = rsRankSyng (rs, target, 0, 0) ;
	    if (actual != expected)
	      { rsPrint (rs) ;
		die ("V %d callCount %d rank(%d,0) mismatch: expected %d, got %d",
		     nVertex, callCount, target, expected, actual) ;
	      }
	    totRank0 += actual ;
	  }
	if (isFind)
	  { I32 expectedSym = 0 ;
	    int expectedRank = 0, sum = 0 ;
	    int target = len/2 ;
	    for (i = 0 ; i < nRuns ; ++i)
	      { sum += C[i] ;
		if (sum > target)
		  { expectedSym = B[i] ;
		    for (int j = 0 ; j <= i ; ++j)
		      if (B[j] == expectedSym) expectedRank += C[j] ;
		    expectedRank -= (sum - target) ;
		    break ;
		  }
	      }
	    I32 actualSym ; U32 actualOffset ;
	    int actualRank = rsFindSyng (rs, target, &actualSym, &actualOffset) ;
	    if (actualSym != expectedSym || actualRank != expectedRank)
	      { rsPrint (rs) ;
		die ("V %d callCount %d rsFind(%d) mismatch: expected sym=%d rank=%d, got sym=%d rank=%d",
		     nVertex, callCount, target, expectedSym,  expectedRank, actualSym, actualRank) ;
	      }
	    totFind += actualRank ;
	  }
	// rsDestroy (rs) ;
	aRS[callCount-1] = rs ;
        break ;
      }
  oneFileClose (of) ;

  if (callCount != nC) die ("nC %d != callCount %d", (int)nC, (int)callCount) ;

  if (nLinear)
    printf ("linear <= %d: %'lld totRuns %'lld avRuns %.1f totLen %'lld avLen %.1f avSym %.1f nBigCount %'lld\n",
	    MAX_LINEAR, (long long)nLinear, (long long)totRunsLinear, totRunsLinear/(1.*nLinear),
	    (long long)totLenLinear, totLenLinear/(1.*nLinear), totSymLinear/(1.*nLinear),
	    (long long)nBigLinear) ;
  if (nFixed)
    printf ("fixed: %'lld totRuns %'lld avRuns %.1f totLen %'lld avLen %.1f avSym %.1f\n",
	    (long long)nFixed, (long long)totRunsFixed, totRunsFixed/(1.*nFixed),
	    (long long)totLenFixed, totLenFixed/(1.*nFixed), totSymFixed/(1.*nFixed)) ;
  if (nDynamic)
    printf ("dynamic: %'lld totRuns %'lld avRuns %.1f totLen %'lld avLen %.1f avSym %.1f\n",
	    (long long)nDynamic, (long long)totRunsDynamic, totRunsDynamic/(1.*nDynamic),
	    (long long)totLenDynamic, totLenDynamic/(1.*nDynamic), totSymDynamic/(1.*nDynamic)) ;
  double fac = 1./(1<<30) ;
  printf ("GBWT SPACE ~= %.3f linear + %.3f fixed + %.3f dynamic = %.3f GB\n",
	  fac * totMaxLinear * sizeof(Linear),
	  fac * totMaxFixed * sizeof(Fixed),
	  fac * totMaxDynamic * sizeof(Dynamic),
	  fac*(totMaxLinear*sizeof(Linear)+totMaxFixed*sizeof(Fixed)+totMaxDynamic*sizeof(Dynamic))) ;

  double total = nLinear + nFixed + nDynamic ;
  if (isCount) printf ("Total rsCount(0) = %'lld average %.1f\n", 
                       (long long)totCount0, totCount0/total) ;
  if (isRank)  printf ("Total rsRank(len/2,0) = %'lld average %.1f\n", 
                       (long long)totRank0, totRank0/total) ;
  if (isFind)  printf ("Total rsFind(len/2,&symbol) = %'lld average %.1f\n", 
                       (long long)totFind, totFind/total) ;
  
  if (isTime)
    { timeUpdate (stdout) ;
      I64 n = 0, tot = 0 ;
      int j ;
      for (i = 0 ; i < nC ; ++i)
	if (aRS[i].linear && aRS[i].linear->max)
	  for (j = 0 ; j < aRS[i].linear->nSym ; ++j)
	    { tot += rsCount (aRS[i], j) ; ++n ; }
      printf ("linear count %'lld %'lld : ", (long long) n, (long long)tot) ; timeUpdate (stdout) ;
      n = tot = 0 ;
      for (i = 0 ; i < nC ; ++i)
	if (aRS[i].linear && aRS[i].linear->max)
	  for (j = 0 ; j < lens[i] ; ++j)
	    { tot += rsFind (aRS[i], j, 0) ; ++n ; }
      printf ("linear find %'lld %'lld : ", (long long) n, (long long)tot) ; timeUpdate (stdout) ;
      n = tot = 0 ;
      for (i = 0 ; i < nC ; ++i)
	if (aRS[i].linear && !aRS[i].linear->max)
	  for (j = 0 ; j < aRS[i].fixed->nSym ; ++j)
	    { tot += rsCount (aRS[i], j) ; ++n ; }					       
      printf ("skip count %'lld %'lld : ", (long long) n, (long long)tot) ; timeUpdate (stdout) ;
      n = tot = 0 ;
      for (i = 0 ; i < nC ; ++i)
	if (aRS[i].linear && !aRS[i].linear->max)
	  for (j = 0 ; j < lens[i] ; ++j)
	    { tot += rsFind (aRS[i], j, 0) ; ++n ; }
      printf ("skip find %'lld %'lld : ", (long long) n, (long long)tot) ; timeUpdate (stdout) ;
    }

  // free stored items
  if (!isBlank) for (i = 0 ; i < nC ; ++i) rsDestroy (aRS[i]) ;
  newFree (aRS, nC, Rskip) ;
  
  //
  //  i = 0 ; n = 1 ;
  //  printf ("           n\tcountMax\t     run\t     len\n") ;
  //  while (i < 25)
  //    { printf ("%12d\t%8lld\t%8lld\n", n, (long long)runHist[i], (long long)lenHist[i]) ;
  //    ++i ; n <<= 1 ;
  //    }

  timeTotal (stderr) ;
  return 0 ;
}

#endif

// end of file
