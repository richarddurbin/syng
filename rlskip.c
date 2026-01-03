/*  File: rlskip.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: code for run-length encoded skip lists
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  2 10:07 2026 (rd109)
 * Created: Sun Nov 30 21:42:51 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

// Use a linear model with cumulative scores until MAX_LINEAR long (default 16).
// Beyond that we make RLSkip lists, using Build32 when building, or Node32 when fixed.
// Note that this limits us to maximum 2<<32 - 1 cells in the skip list.

// We will work with a list of nodes. The first one is special, holding other information
// which is referred to by union fields max, start, nSym and free.
// If rls.linear->max then Linear, else if rls.build32->max than Build32, else Node32.

// For the skiplists, we share nodes between an overall skiplist and embedded ones for each symbol.
// n->right points to the next right, and n->sRight to the next right of the same symbol.
// When building we store counts n->count and n->sCount
// else we store partial sums n->sum and n->sSum.
// The symbol and the down pointer share a slot, since the symbol is only needed at level 0.
// This is distinguished by n->sCount == 0 when building, else by n->sSum & FLAG.
// When building we also need an up pointer, for reasons that will become clear below.

// The first sMax nodes of the skiplists are the start nodes for symbols 0..sMax-1
// (this saves on an extra set of pointers at the cost of a little gymnastics).

static const int MAX_LINEAR = 16 ;
static const U32 FLAG       = (1<<31) ;
static const U32 MASK       = ((unsigned)1<<31)-1 ;
static const U32 TERM       = ((long long)1<<30)-1 ; // end sentinel for ->right pointers

typedef union {
  struct {
    union { U8 sym ; U8 max ; } ;
    union { U8 count ; U8 maxSym ; } ;   // if count == ff then take next node as U16 count
  } ;
  U16 bigCount ;
} LinearNode ;

typedef struct {
  U32 count ;
  union { U32 sCount ; U32 max ; } ; // sCount==0 indicates bottom row, since there sCount==count
  union { U32 right ; U32 maxDepth ; } ;
  union { U32 sRight ; U32 start ; } ;
  union { U32 up ; U32 nSym ; } ;
  union { U32 sym ; U32 down ; U32 free ; } ; // rls.build32->free is the first free node
} Build32 ;			// version for building

typedef struct {
  U32 sum, sSum ; // sSum & FLAG indicates bottom row, in which case mask it with MASK
  union { U32 right ; U32 max ; } ;
  union { U32 sRight ; U32 start ; } ; 
  union { U32 sym ; U32 down ; } ;
} Node32 ; 			// locked node

typedef union { LinearNode *linear ; Node32 *node32 ; Build32 *build32 ; } RLSkip  ;

// start for each symbol i is i+1'th node

#define RLSKIP_DEFINED

#include "rlskip.h" // for interface declarations

  
int rlsDestroy (RLSkip rls)
{
  if (!rls.linear) return RLS_EMPTY ;
  if (rls.linear->max) newFree (rls.linear, rls.linear->max+1, LinearNode) ;
  else if (rls.build32->max) newFree (rls.build32, rls.build32->max+1, Build32) ;
  else newFree (rls.node32, rls.node32->max+1, Node32) ;
  return 0 ;
}

int rlsSize (RLSkip rls)
{
  if (!rls.linear) return 0 ;
  if (rls.linear->max) return (rls.linear->max+1) * sizeof (LinearNode) ;
  else if (rls.build32->max) return (rls.build32->max+1) * sizeof (Build32) ;
  else return (rls.node32->max+1) * sizeof (Node32) ;
}

// next a very simple pseudo-random number generator - good enough for us here

#include <stdlib.h> // for random number seed
static U64 rng = 18982392197 ; // a prime - each thread's copy starts with this

static inline bool eFlip (void)
{ rng *= 479 ; // a prime
  return (((rng >> 20) & 0x7) < 3) ; // deterministic pseudorandom with probability 3/8 ~= 1/e 
}

// some debugging functions

static void inline printNode32 (Node32 *n, U32 i)
{
  if (n->sSum & FLAG) // at bottom
    printf (" %u: sym %u sum %u right %u sSum %u sRight %u\n",
           i, n->sym, n->sum, n->right, n->sSum & MASK, n->sRight) ;
  else
    printf (" %u:        sum %u right %u sSum %u sRight %u down %u\n",
           i, n->sum, n->right, n->sSum, n->sRight, n->down) ;
}

static void rlsPrint (RLSkip rls)
{
  if (rls.linear->max)
    { U32 i ;
      LinearNode *nd = rls.linear ; nd++ ;
      printf ("LinearNode max %u maxSym %u\n", rls.linear->max, rls.linear->maxSym) ;
      int j = 0, sum = 0, sSum[rls.linear->maxSym] ; memset (sSum, 0, rls.linear->maxSym*sizeof(int)) ;
      for (i = 0 ; i < rls.linear->max ; ++i, ++nd)
        { int count = (nd->count == 255) ? nd[1].bigCount : nd->count ;
          sum += count ;
          sSum[nd->sym] += count ;
          printf (" %d: i %u sym %u count %u sum %u sSum %u\n", ++j, i, nd->sym, count, sum, sSum[nd->sym]) ;
          if (nd->count == 255) { ++i ; ++nd ; }
        }
    }
  else if (rls.build32->max) 
    { ; // not implemented yet
    }
  else
    { U32 i ;
      Node32 *nd = rls.node32 ; nd++ ;
      printf ("Node32 start %u max %u\n", rls.node32->start, rls.node32->max) ;
      for (i = 0 ; i < rls.node32->max ; ++i, ++nd)
        printNode32 (nd, i) ;
    }
}

static U32 DEBUG_COUNT = -1 ;

RLSkip rlsBuildFromI64 (int n, I64 *sym, I64 *runLen) // n is count of (sym, runLen) pairs
// this is used when building from data stored in a .1gbwt OneFile
{ static int callCount = 0 ; ++callCount ;
  int i, max = n ;
  if (max < MAX_LINEAR) // build in Linear nodes
    { for (i = 0 ; i < n ; ++i) if (runLen[i] >= 255) ++max ;
      if (max < MAX_LINEAR)
	{ LinearNode *node = new0 (max+1, LinearNode) ;
	  RLSkip rls ; rls.linear = node++ ;
	  rls.linear->max = max ;
	  for (i = 0 ; i < n ; ++i, ++node)
	    { node->sym = sym[i] ;
	      if (sym[i] >= rls.linear->maxSym) rls.linear->maxSym = sym[i]+1 ;
	      if (runLen[i] < 255)
		node->count = runLen[i] ;
	      else
		{ node->count = 255 ; node++ ;
		  node->bigCount = runLen[i] ;
		}
	    }
	  return rls ;
	}
    }
  // if we get here, then max >= MAX_LINEAR - make a skipList

  // first count number of nodes needed, and some other things we will need
  max = 0 ;       // start again to count nodes needed
  int maxDepth = 0, maxSym = -1 ;
  static int maxMaxSym = 0 ;
  for (i = 0 ; i < n ; ++i) if (sym[i] > maxSym) maxSym = sym[i] ;
  ++maxSym ; // now maxSym is number of symbols
  if (maxSym > maxMaxSym) maxMaxSym = maxSym ;
  int iDepth[n] ; // number of nodes/layers for run i
  int s ;         // a symbol
  int sFirst[maxSym] ; // the lowest run index i at which symbol s occurs
  int sMaxDepth[maxSym] ; // max depth for symbol s
  for (s = 0 ; s < maxSym ; ++s) { sMaxDepth[s] = 0 ; sFirst[s] = -1 ; }
  for (i = n ; i-- ; ) // must be backwards so that sFirst[] works
    { sFirst[sym[i]] = i ;
      iDepth[i] = 1 ; while (eFlip ()) ++iDepth[i] ; // skiplist randomness
      max += iDepth[i] ;
      if (iDepth[i] > maxDepth) maxDepth = iDepth[i] ;
      if (iDepth[i] > sMaxDepth[sym[i]]) sMaxDepth[sym[i]] = iDepth[i] ;
    }
  for (s = 0 ; s < maxSym ; ++s)
    { if (sFirst[s] == -1) die ("unused symbol %d < %d", s, maxSym) ; // confirm no unused symbols
      // now ensure that the first column for s has maximal depth of all columns for s
      max += sMaxDepth[s] - iDepth[sFirst[s]] ; iDepth[sFirst[s]] = sMaxDepth[s] ;
    }
  // and similarly that the first column of all has maximal depth of all columns
  max += maxDepth - iDepth[0] ; iDepth[0] = maxDepth ; sMaxDepth[sym[0]] = maxDepth ;

  // next make the rls object with the right amount of memory
  Node32 *node = new0 (max+1, Node32) ;
  RLSkip rls ; rls.node32 = node++ ;
  
  // preliminaries for building it
  int stack[maxDepth] ; // we will need these to set up right pointers
  int sStack[maxSym][maxDepth] ;
  int d, r ;   // generic variables for a depth and a run length
  Node32 *nd ; // generic node pointer

  // initialise stacks to point to node[-1], which is the header node
  // we will use this as a generic empty node with sum = sSum = 0
  // this allows the update for the first node for each symbol skiplist to be 
  //   the same as for all others in the lines below starting nd->sum = and nd->sSum =
  for (d = 0 ; d < maxDepth ; ++d)
    { stack[d] = -1 ;
      for (s = 0 ; s < maxSym ; ++s) sStack[s][d] = -1 ;
    }

  // allocate nodes upwards, starting after the symbols, so last element is allocated last
  int freeNode = maxSym ; 
  
  // now loop over the runs - overall nlog_n complexity
  for (i = 0 ; i < n ; ++i)
    { s = sym[i] ;
      r = runLen[i] ;
      U32 next = (i == sFirst[s]) ? s : freeNode++ ;
      for (d = iDepth[i] ; d-- ; )
	{ nd = &node[next] ;
	  nd->sum = node[stack[d]].sum ; node[stack[d]].right = next ; stack[d] = next ;
	  nd->sSum = node[sStack[s][d]].sSum ; node[sStack[s][d]].sRight = next ; sStack[s][d] = next ;
          nd->right = TERM ; nd->sRight = TERM ;
          if (d) { next = freeNode++ ; nd->down = next ; }
	}
      nd->sSum |= FLAG ; // marker for bottom of the column
      nd->sym = s ;
      // now go back up the stacks adding on r
      for (d = 0 ; d < maxDepth ; ++d) node[stack[d]].sum += r ;
      for (d = 0 ; d < sMaxDepth[s] ; ++d) node[sStack[s][d]].sSum += r ;
    }

  rls.node32->max = max ;      // have to set these after building
  rls.node32->start = sym[0] ; // since we used them in stack initialisation
  return rls ;
}

int rlsCount (RLSkip rls, U32 symbol)	// how many of symbol
{
  static int callCount = 0 ; ++callCount ;
  // fprintf (stderr, "DEBUG rlsCount call %d for symbol %u\n", callCount, symbol) ;
  if (rls.linear->max) // linear
    { int i, sSum = 0 ;
      LinearNode *node = rls.linear ; node++ ;
      for (i = 0 ; i < rls.linear->max ; ++i, ++node)
	if (node->count < 255) { if (node->sym == symbol) sSum += node->count ; }
	else { ++i ; if ((node++)->sym == symbol) sSum += node->bigCount ; }
      return sSum ;
   
}
  else if (rls.build32->max) // build32
    { U32 i, sSum = 0 ;
      Build32 *node = rls.build32 ; node++ ;
      for (i = symbol ; node[i].sRight != TERM ; i = node[i].sRight) 
        sSum += node[i].sCount ;
      return sSum + node[i].sCount ;
    }
  else // node32
    { U32 i ;
      Node32 *node = rls.node32 ; node++ ;
      for (i = symbol ; node[i].sRight != TERM ; i = node[i].sRight) { ; }
      return node[i].sSum & MASK ;
    }
}

int rlsRank (RLSkip rls, U32 k, U32 symbol)  // how many of symbol up to (not including) k
{
  static int callCount = 0 ; ++callCount ;
  // fprintf (stderr, "DEBUG rlsRank call %d for k=%u symbol=%u\n", callCount, k, symbol) ;
  if (callCount == DEBUG_COUNT) rlsPrint (rls) ;
  if (rls.linear->max) // linear
    { int i, sum = 0, sSum = 0 ;
      LinearNode *node = rls.linear ; node++ ;
      for (i = 0 ; i < rls.linear->max ; ++i, ++node)
	{ int sym = node->sym, count = node->count ;
	  if (count == 255) { count = (++node)->bigCount ; ++i ; }
	  sum += count ;
	  if (sym == symbol) sSum += count ;
	  if (sum >= k) return (sym == symbol) ? sSum - (sum - k) : sSum ;
	}
    }
  else if (rls.build32->max) // build32 - must work in the main list
    { U32 sum = 0, sSum = rlsCount (rls, symbol) ; // LOGARITHMIC
      Build32 *node = rls.build32 ; node++ ;
      U32 i = rls.build32->start ;
      while (true) // LOGARITHMIC
	if (sum + node[i].count < k)
	  { sum = sum + node[i].count ;
	    i = node[i].right ;
	    if (i == TERM) return sSum ;
	  }
	else if (node[i].sCount)
	  i = node[i].down ;
	else
	  break ; // at bottom in correct node
      if (node[i].sym == symbol)
	sSum -= (sum + node[i].count - k) ;
      else // go right until we find symbol
	for (i = node[i].right ; i ; i = node[i].right) // ALPHABET_LINEAR
	  if (node[i].sym == symbol) break ;
      while (i != TERM) // now subtract any remaining counts - LOGARITHMIC
	if (node[i].up) i = node[i].up ;
	else
	  { sSum -= node[i].sCount ;
	    i = node[i].right ;
	  }
      return sSum ;
    }
  else // node32 - can follow the symbol sub-list, with care
    { Node32 *node = rls.node32 ; node++ ;
      U32 i = symbol ;
      int sSumLast = 0 ;
      while (true) // LOGARITHMIC
        if (node[i].sum < k)
	  { U32 j = node[i].sRight ;
            if (j == TERM) 
              { if (!(node[i].sSum & FLAG)) i = node[i].down ; 
                else return node[i].sSum & MASK ; 
              }
            else if (node[j].sum < k || (node[i].sSum & FLAG)) 
              { sSumLast = node[i].sSum & MASK ; 
                i = j ; 
        }
            else 
              i = node[i].down ;
          }
	else if (!(node[i].sSum & FLAG))
	  i = node[i].down ;
	else
	  break ; // at bottom in correct node
      int val = (node[i].sSum & MASK) - (node[i].sum - k) ; 
      if (val < sSumLast) // i.e. k is before this cell
	return sSumLast ;
      else
	return val ;
    }
  die ("shouldn't get here in rlsRank") ; return -1 ;
}

int rlsFind (RLSkip rls, U32 k, U32 *symbol) // quite like rank()
{
  if (rls.linear->max) // linear
    { U32 i, sum = 0, sSum = 0 ;
      LinearNode *node = rls.linear ; node++ ;
      int sSums[rls.linear->maxSym] ; memset (sSums, 0, rls.linear->maxSym * sizeof(int)) ;
      for (i = 0 ; i < rls.linear->max ; ++i, ++node)
	{ int sym = node->sym, count = node->count ;
	  if (count == 255) { count = (++node)->bigCount ; ++i ; }
          sum += count ;
	  sSums[sym] += count ;
	  if (sum >= k)
	    { if (symbol) *symbol = sym ;
	      return sSums[sym] - (sum - k) ;
	    }
	}
    }
  else if (rls.build32->max) // build32
    { U32 sum = 0 ; // LOGARITHMIC
      Build32 *node = rls.build32 ; node++ ;
      U32 i = rls.build32->start ;
      while (true) // LOGARITHMIC
	if (sum + node[i].count < k)
	  { sum = sum + node[i].count ;
	    i = node[i].right ;
	    if (i == TERM) die ("input k = %u > max in rlsFind() node32", k) ;
	  }
	else if (node[i].sCount) i = node[i].down ;
	else break ; // at bottom in correct node
      if (symbol) *symbol = node[i].sym ;
      U32 sSum = rlsCount (rls, node[i].sym) ; // LOGARITHMIC
      sSum -= (sum + node[i].count - k) ;
      while (i != TERM) // now subtract any remaining counts - LOGARITHMIC
	if (node[i].up) 
          i = node[i].up ;
	else
	  { sSum -= node[i].sCount ;
	    i = node[i].right ;
	  }
      return sSum ;
    }
  else // node32
    { Node32 *node = rls.node32 ; node++ ;
      U32 i = rls.node32->start ;
      while (true) // LOGARITHMIC
	if (node[i].sum < k)
          { i = node[i].right ;
	    if (i == TERM) die ("input k = %u > max in rlsFind() node32, i %u TERM %u", k) ;
          }
	else if (!(node[i].sSum & FLAG))
	  i = node[i].down ;
	else break ; // at bottom in correct node
      if (symbol) *symbol = node[i].sym ;
      return (node[i].sSum & MASK) - (node[i].sum - k) ;
    }
  die ("shouldn't get here in rlsFind") ; return -1 ;
}

#ifdef OLD

static int map32 (RLSkip *rls, U32 i, U32 *symbol, bool isFind)
{ // used for both rank() and find() when building
  // need to first find i in the full skipList,
  // then find the next node matching symbol at each level (stored in post[])
  // then run the symbol's skipList using sRight until post[] then down, summing sCount
  // in total three log time steps, plus for rank() a linear search for symbol if necessary
  Build32 *x = arrp(rls->a, rls->start, Build32) ; // search full skipList
  while (true)
    if (!x)
      return ERROR_INDEX_OVERFLOW ;
    else if (i >= x->count)
      { i -= x->count ;
	x = arrp(rls->a, x->right, Build32) ;
      }
    else if (x->sCount) // not yet at row 0
      x = arrp(rls->a, x->down, Build32) ;
  // now at the right node

  bool isSymbolMatch = true ; // default
  if (isFind) *symbol = x->symbol ;
  else if (x->symbol != *symbol) // linear search to find a matching symbol - E(t) ~ A
    { isSymbolMatch = false ;
      while (x && x->symbol != *symbol) x = arrp(rls->a, x->right, Build32) ;
      if (!x) return rlsCount (rls, *symbol) ; // no matches - we are done
    }
  
  // next build the post list
  int maxDepth = 0 ; // first find what depth I need
  { Build32 *y = arrp(rls->a, arr(rls->first, *symbol, U32), Build32) ;
    while (y->sCount) { ++maxDepth ; y = arrp(rls->a, y->down, Build32) ; }
  }
  // introduce a new block so as to declare post[] on the stack
  { U32 post[maxDepth+1] ;
    int depth = 0 ;
    while (x)
      { post[depth++] = x - arrp(rls->a, 0, Build32) ;
	while (x && !x->up) x = arrp(rls->a, x->sRight, Build32) ;
	if (x) x = arrp(rls->a, x->up, Build32) ;
      }
    while (depth <= maxDepth) post[depth++] = 0 ;
    // OK - we have our stack, with known depth.  Run back down.
    
    // finally search 
    Build32 *z = arrp(rls->a, arr(rls->first, *symbol, U32), Build32) ;
    depth = maxDepth ; // this is the depth of z, by construction
    int sum = 0 ;
    while (true)
      if (z->sRight != post[depth])
	{ sum += z->sCount ? z->sCount : z->count ;
	  z = arrp(rls->a, z->sRight, Build32) ;
	}
      else if (z->sCount) // this means that we are not at the bottom - also depth > 1
	{ z = arrp(rls->a, z->down, Build32) ;
	  --depth ;
	}
    // now we are done
    if (isSymbolMatch)
      { z = arrp(rls->a, z->sRight, Build32) ;
	return sum + i ;
      }
    else
      return sum ;
  }
}

int rlskipRank (RLSkip *rls, U32 i, U32 symbol)  // return how many symbol up to and before i
{
  if (!rls) return ERROR_NO_RLS ;
  if (symbol >= arrayMax(rls->first)) return ERROR_SYMBOL_OVERFLOW ;
  if (rls->isLocked) // single scan in symbol's skipList - super efficient
    { Node32 *x = arrp(rls->a, arr(rls->first, symbol, U32), Node32) ;
      Node32 *last = x ;
      while (true)
	if (!x)
	  return ERROR_INDEX_OVERFLOW ;
	else if (i >= x->sum)
	  { last = x ;
	    x = arrp(rls->a, x->sRight, Node32) ;
	  }
	else if (x->sSum & 0x80000000) // not yet at row 0
	  x = arrp(rls->a, x->down, Node32) ;
      if (x->sSum + i - x->sum > last->sSum) return (int)(x->sSum + i - x->sum) ; // i is within x
      else return (int)last->sSum ; // i is between last and x
    }
  else // three scans needed, plus potentially a linear search of expected size A
    return map32 (rls, i, &symbol, false) ;
}

int  rlskipFind (RLSkip *rls, U32 i, U32 *symbol) // fills symbol at i and returns rank of it
{
  if (!rls) return ERROR_NO_RLS ;
  if (*symbol >= arrayMax(rls->first)) return ERROR_SYMBOL_OVERFLOW ;
  if (rls->isLocked) // single scan in full skilList - logarithmically efficient
    { Node32 *x = arrp(rls->a, rls->start, Node32) ;
      while (true)
	if (!x)
	  return ERROR_INDEX_OVERFLOW ;
	else if (i >= x->sum)
	  x = arrp(rls->a, x->sRight, Node32) ;
	else if (x->sSum & 0x80000000) // not yet at row 0
	  x = arrp(rls->a, x->down, Node32) ;
      *symbol = x->symbol ;
      return (int)(x->sSum + i - x->sum) ;
    }
  else // three scans needed, but no linear search
    return map32 (rls, i, symbol, true) ;
}

int  rlskipAdd (RLSkip *rls, U32 i, U32 symbol)  // add symbol at i and return rank
{
  if (!rls) return ERROR_NO_RLS ;
  if (rls->isLocked) return ERROR_LOCKED ;
  return 0 ;
}

int rlskipLock (RLSkip *rls) // pack memory and reduce
{ die ("not implemented yet") ;
  if (!rls) return ERROR_NO_RLS ;
  if (rls->isLocked) return 0 ; // we are done
  
  rls->isLocked = true ;
  return 0 ;
}

char *rlskipError (RLSkip *rls, int err) // returns error text for last error
{ thread_local char errBuf[256] ;
  switch (err)
    {
    RLS_EMPTY:
      return "rls pointer is null" ;
    RLS_SYMBOL_OVERFLOW:
      sprintf (errBuf, "symbol >= max %d", (int)arrayMax(rls->first)) ; break ;
    RLS_INDEX_OVERFLOW:
      sprintf (errBuf, "index >= max %d", (int)arrayMax(rls->a)) ; break ;
    RLS_LOCKED:
      return "attempt to modify a locked rls object" ;
    }
  return errBuf ;
}

#endif // OLD

/***************************** test package ***************************/

#ifdef TEST

// compile with: cc -g -DTEST -o rltest rlskip.c array.c utils.c ONElib.c -lz  

#include "ONElib.h"
#include "array.h"

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  --argc ; ++argv ;
  bool isBuild = false, isCount = false, isRank = false, isFind = false ;
  while (argc && *argv[0] == '-')
    if (!strcmp (argv[0], "-build")) { isBuild = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-count")) { isCount = true ; isBuild = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-rank")) { isRank = true ; isBuild = true ; --argc ; ++argv ; }
    else if (!strcmp (argv[0], "-find")) { isFind = true ; isBuild = true ; --argc ; ++argv ; }
    else die ("unknown option %s", argv[0]) ;
  if (argc != 1) die ("%s XX.1gbwt", argv[-1]) ;

  OneFile *of = oneFileOpenRead (argv[0], 0, "gbwt", 1) ;
  if (!of || strcmp (of->subType, "gbwt")) die ("can't open .1gbwt file %s", argv[0]) ;
  if (!oneGoto (of, 'V', 1)) die ("can't locate to a graph node 'V' line in %s", of->fileName) ;

  U64 *runHist = new0(32, U64), *lenHist = new0(32, U64) ;
  I64 *B, *C ;
  U64 nLinear = 0, totRunsLinear = 0, totLenLinear = 0, maxLenLinear = 0, nBigLinear = 0 ;
  U64 nSmall = 0, totRunsSmall = 0, totLenSmall = 0, maxRunsSmall = 0 ;
  U64 nMedium = 0, totRunsMedium = 0, totLenMedium = 0, maxRunsMedium = 0 ;
  U64 nLarge = 0, totRunsLarge = 0, totLenLarge = 0, maxRunsLarge = 0 ;
  U64 totSize = 0, totCount0 = 0, totRank0 = 0, totFind = 0 ;
  int i, n, count = 0 ;
  while (oneReadLine (of))
    switch (of->lineType)
      {
      case 'B': case 'b': B = oneIntList(of) ; break ;
      case 'C': case 'c':
        ++count ;
	C = oneIntList(of) ;
	int nRuns = oneLen(of) ;
	int kRuns = 0 ; n = nRuns ; while (n >>= 1) ++kRuns ;
	int len = 0 ;
	for (i = 0 ; i < nRuns ; ++i) len += C[i] ;
	int kLen = 0 ; n = len ; while (n >>= 1) ++kLen ;
	if (nRuns < MAX_LINEAR) // linear option: 2 bytes per node (symbol + count)
	  { ++nLinear ;
	    totRunsLinear += nRuns ;
	    totLenLinear += len ;
	    if (len > maxLenLinear) maxLenLinear = len ;
	    for (i = 0 ; i < nRuns ; ++i) if (C[i] > 256) ++nBigLinear ;
	  }
/* // comment out this section for now
	else if (len < 256) // byte option: 5/6 bytes per node
	  { ++nSmall ;
	    totRunsSmall += nRuns ;
	    totLenSmall += len ;
	    if (nRuns > maxRunsSmall) maxRunsSmall = nRuns ;
	  }
	else if (len < (1 << 16)) // 2-byte option: 10/12 bytes per node
	  { ++nMedium ;
	    totRunsMedium += nRuns ;
	    totLenMedium += len ;
	    if (nRuns > maxRunsMedium) maxRunsMedium = nRuns ;
	  }
*/
	else // 4-byte option
	  { ++nLarge ;
	    totRunsLarge += nRuns ;
	    totLenLarge += len ;
	    if (nRuns > maxRunsLarge) maxRunsLarge = nRuns ;
	  }
        RLSkip rls ;
        if (isBuild)
          { rls = rlsBuildFromI64 (nRuns, B, C) ;
            totSize += rlsSize (rls) ;

            // Verify correctness by computing expected values from B and C arrays
            if (isCount)
              { int expected = 0 ;
                for (i = 0 ; i < nRuns ; ++i) if (B[i] == 0) expected += C[i] ;
                int actual = rlsCount (rls, 0) ;
                if (actual != expected)
                  die ("rlsCount(0) mismatch: expected %d, got %d", expected, actual) ;
                totCount0 += actual ;
              }
            if (isRank)
              { int expected = 0, sum = 0 ;
                for (i = 0 ; i < nRuns ; ++i)
                  { if (B[i] == 0) expected += C[i] ;
                    sum += C[i] ;
                    if (sum >= len/2)
                      { if (B[i] == 0) expected -= (sum - len/2) ;
                        break ;
                      }
                  }
                int actual = rlsRank (rls, len/2, 0) ;
                if (actual != expected)
                  { fprintf(stderr, "DEBUG: count=%d nRuns=%d len=%d k=%d\n", count, nRuns, len, len/2) ;
                    fprintf(stderr, "B: ") ; for (i = 0 ; i < nRuns ; ++i) fprintf(stderr, "%2lld ", B[i]) ; fprintf(stderr, "\n") ;
                    fprintf(stderr, "C: ") ; for (i = 0 ; i < nRuns ; ++i) fprintf(stderr, "%2lld ", C[i]) ; fprintf(stderr, "\n") ;
                    die ("rlsRank(%d,0) mismatch: expected %d, got %d", len/2, expected, actual) ;
                  }
                totRank0 += actual ;
              }
            if (isFind)
              { U32 expectedSym = 0 ;
                int expectedRank = 0, sum = 0 ;
                for (i = 0 ; i < nRuns ; ++i)
                  { sum += C[i] ;
                    if (sum >= len/2)
                      { expectedSym = B[i] ;
                        for (int j = 0 ; j <= i ; ++j)
                          if (B[j] == expectedSym) expectedRank += C[j] ;
                        expectedRank -= (sum - len/2) ;
                        break ;
                      }
                  }
                U32 actualSym ;
                int actualRank = rlsFind (rls, len/2, &actualSym) ;
                if (actualSym != expectedSym || actualRank != expectedRank)
                  die ("rlsFind(%d) mismatch: expected sym=%u rank=%d, got sym=%u rank=%d",
                       len/2, expectedSym, expectedRank, actualSym, actualRank) ;
                totFind += actualRank ;
              }
            rlsDestroy (rls) ;
          }
        break ;
      }
  oneFileClose (of) ;

  printf ("linear < %d: %lld totRuns %lld av %.1f totLen %lld avLen %.1f nBigCount %lld\n", MAX_LINEAR, 
	  (long long)nLinear, (long long)totRunsLinear, totRunsLinear/(1.*nLinear),
	  (long long)totLenLinear, totLenLinear/(1.*nLinear), (long long)nBigLinear) ;
  printf ("small: %lld totRuns %lld av %.1f totLen %lld avLen %.1f avCount %.1f\n",
	  (long long)nSmall, (long long)totRunsSmall, totRunsSmall/(1.*nSmall),
	  (long long)totLenSmall, totLenSmall/(1.*nSmall), totLenSmall/(1.*totRunsSmall)) ;
  printf ("medium: %lld totRuns %lld av %.1f totLen %lld avLen %.1f avCount %.1f\n",
	  (long long)nMedium, (long long)totRunsMedium, totRunsMedium/(1.*nMedium),
	  (long long)totLenMedium, totLenMedium/(1.*nMedium), totLenMedium/(1.*totRunsMedium)) ;
  printf ("large: %lld totRuns %lld av %.1f totLen %lld avLen %.1f avCount %.1f\n",
	  (long long)nLarge, (long long)totRunsLarge, totRunsLarge/(1.*nLarge),
	  (long long)totLenLarge, totLenLarge/(1.*nLarge), totLenLarge/(1.*totRunsLarge)) ;
  double fac = 1./(1<<30) ;
  printf ("GBWT SPACE ~= %.2f linear + %.2f small + %.2f medium + %.2f large = %.2f GB\n",
	  fac*(2*totRunsLinear + nBigLinear),
	  8*fac*totRunsSmall, 16*fac*totRunsMedium, 32*fac*totRunsLarge,
	  fac*(2*totRunsLinear + nBigLinear + 8*totRunsSmall + 16*totRunsMedium + 32*totRunsLarge)) ;
  double total = nLinear + nSmall + nMedium + nLarge ;
  double totRuns = totRunsLinear + totRunsSmall + totRunsMedium + totRunsLarge ;
  if (isBuild) printf ("Total RLSkip size = %.2f GB\n", totSize * fac) ;
  if (isCount) printf ("Total rlsCount(0) = %lld average %.1f\n", 
                       (long long)totCount0, totCount0/total) ;
  if (isRank)  printf ("Total rlsRank(len/2,0) = %lld average %.1f\n", 
                       (long long)totRank0, totRank0/total) ;
  if (isFind)  printf ("Total rlsFind(len/2,&symbol) = %lld average %.1f\n", 
                       (long long)totFind, totFind/total) ;
  
  //
  //  i = 0 ; n = 1 ;
  //  printf ("           n\tcountMax\t     run\t     len\n") ;
  //  while (i < 25)
  //    { printf ("%12d\t%8lld\t%8lld\t%8lld\n", n, (long long)eHist[i], (long long)rHist[i], (long long)lenHist[i]) ;
  //    ++i ; n <<= 1 ;
  //    }

  timeTotal (stderr) ;
  return 0 ;
}

#endif

// end of file
