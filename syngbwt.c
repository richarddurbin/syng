/*  File: syngbwt.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 23 08:54 2025 (rd109)
 * * Apr  5 11:38 2025 (rd109): added FM, threaded read and write, sort on write
 * Created: Mon Sep  9 11:34:51 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "syng.h"

typedef struct {
  I32 in, out ;           // syncs
  I32 inCount, outCount ; // number of paths entering from in (+ve dirn), outCount from out (-ve dirn)
  I32 inOff, outOff ;     // sync offsets
} SimpleNode ;

typedef struct {
  I32 sync ;
  I32 count ;
} SyncCount ;  // also used for run-length encoded GBWT

typedef struct {
  bool         is2 : 1 ;
  unsigned int cum : 31 ;
  union {
    struct { unsigned int index1 : 6  ; unsigned int occ1 : 26 ; } ;
    struct { unsigned int index2 : 16 ; unsigned int occ2 : 16 ; } ; 
  } fm ;
} FM ;         // alternate used for FM-encoded GBWT if syngBWT->isFM

static inline I32 fmIndex (FM *f) { return (I32) (f->is2 ? f->fm.index2 : f->fm.index1) ; }
static inline I32 fmOcc (FM *f)   { return (I32) (f->is2 ? f->fm.occ2 : f->fm.occ1) ; }

typedef struct {
  I32        inN, outN ;   // how many distinct syncs coming in, going out
  I32        inG, outG ;   // how big are the run-length encoded GBWTs
  SyncCount *sc ;          // inList, then outList, inGBWT, then outGBWT, then inOff, then outOff
} ComplexNode ;

typedef struct {
  I64        size ;
  U8        *data ;
} PackedNode ;

typedef union {
  SimpleNode   simple ;
  ComplexNode  complex ;
  PackedNode   packed ;
} Node ;

#define NODE_SIMPLE     0x01
#define NODE_COMPLEX    0x02
#define NODE_PACKED     0x04 // not used yet, but packing/unpacking code is below
#define NODE_SKIP8	0x08
#define NODE_SKIP16	0x10
#define NODE_SKIP32	0x20

static Node nodeCreate (I32 k, I32 in, I32 out, I32 inOff,I32 outOff) ;
static I32  nodePack (Node *n, U8 *status) ;   // returns new size in bytes
static I32  nodeUnpack (Node *n, U8 *status) ; // returns new size in bytes

static inline int intPut (U8 *u, I32 val) ;   // at end of file for integer byte packing
static inline int intGet (U8 *u, I32 *pval) ; // at end of file for integer byte unpacking

/*********************** create and destroy ****************************/

SyngBWT *syngBWTcreate (int fixedLen, I64 max)
{
  if (!fixedLen) die ("syngBWT does not yet support variable length sequences") ;
  SyngBWT *sb = new0 (1, SyngBWT) ;
  sb->fixedLen = fixedLen ;
  if (!max) max = 1 << 20 ;
  sb->node = arrayCreate (max, Node) ;
  sb->status = arrayCreate (max, U8) ;
  return sb ;
}

void syngBWTdestroy (SyngBWT *sb)
{
  int i ;
  for (i = 1 ; i < arrayMax (sb->status) ; ++i)
    if (arr(sb->status,i,U8) & NODE_COMPLEX) free (arr(sb->node,i,Node).complex.sc) ;
  arrayDestroy (sb->node) ;
  arrayDestroy (sb->status) ;
  if (sb->length) arrayDestroy (sb->length) ;
  if (sb->startHash) hashDestroy (sb->startHash) ;
  if (sb->startHashCount) arrayDestroy (sb->startHashCount) ;
  newFree (sb, 1, SyngBWT) ;
}

/*********************** add - the main function ***********************/

// first a utility

static inline void scSplit (ComplexNode *cn, SyncCount **inList, SyncCount **outList,
			    SyncCount **inGBWT, SyncCount **outGBWT,
			    I32 **inOffList, I32 **outOffList)
{
  if (!cn->sc) return ; // don't do this more than once
  SyncCount *t ;
  // don't need to reallocated *inList, because it takes over cn->sc
  t = new(cn->outN+1,SyncCount) ; memcpy(t,*outList,cn->outN*sizeof(SyncCount)) ; *outList = t ;
  t = new(cn->inG+2,SyncCount) ;  memcpy(t,*inGBWT,cn->inG*sizeof(SyncCount)) ;   *inGBWT = t ;
  t = new(cn->outG+2,SyncCount) ; memcpy(t,*outGBWT,cn->outG*sizeof(SyncCount)) ; *outGBWT = t ;
  I32 *tt = new (cn->inN+1,I32) ; memcpy(tt, *inOffList,cn->inN*sizeof(I32)) ; *inOffList = tt ;
  tt = new (cn->outN+1,I32) ; memcpy(tt, *outOffList, cn->outN*sizeof(I32)) ; *outOffList = tt ;
  cn->sc = 0 ;
  // do this now, with 0 so it doesn't actually free the memory, since we know size
  newFree (0, cn->inN + 2*(cn->outN + cn->inG + cn->outG) + 5, SyncCount) ;
  newFree (0, (cn->inN + cn->outN + 1)/2, SyncCount) ; // offList old space
  newFree (0, (cn->inN + cn->outN + 3)/2, SyncCount) ; // offList split space
}

static U32 startCount (SyngBWT *sb, I32 k, bool isAdd)
{
  if (!sb->startHash)
    { sb->startHash = hashCreate (8192) ;
      sb->startHashCount = arrayCreate (8192, U32) ;
    }
  int index ;
  
  if (isAdd)
    { hashAdd (sb->startHash, hashInt(k), &index) ;
      return array(sb->startHashCount, index, U32)++ ;
    }
  else
    { if (!hashFind (sb->startHash, hashInt(k), &index)) return 0 ;
      return arr(sb->startHashCount, index, U32) ;
    }
}

static void startCountAdd (SyngBWT *sb, I32 k, I32 count)
{
  if (!sb->startHash)
    { sb->startHash = hashCreate (8192) ;
      sb->startHashCount = arrayCreate (8192, U32) ;
    }
  int index ;
  hashAdd (sb->startHash, hashInt(k), &index) ;
  array(sb->startHashCount, index, U32) += count ;
}

// this is the core operation to build the GBWT, by inserting path (in, k, out) at k
// j is the offset in the list of paths from in to k; returns next j

// #define DEBUG_ADD // use for detailed debugging - very verbose!

static void nodePrint (Node *n, U8 *s, bool isFM)
{
  if (!*s)
    printf ("\n  empty node\n") ;
  else if (*s & NODE_SIMPLE)
    { SimpleNode *sn = &(n->simple) ;
      printf ("\n  simple node: in %d [%d] (%d) out %d [%d] (%d)\n",
	      sn->in, sn->inCount, sn->inOff, sn->out, sn->outCount, sn->outOff) ;
    }
  else if (*s & NODE_COMPLEX)
    { int i ;
      ComplexNode *cn = &(n->complex) ;
      printf ("\n  complex node: in %d out %d inG %d outG %d", cn->inN,cn->outN,cn->inG,cn->outG) ;
      SyncCount *inList = cn->sc, *outList = inList + cn->inN ;
      SyncCount *inGBWT = outList + cn->outN, *outGBWT = inGBWT + cn->inG ;
      I32 *inOffList = (I32*) (outGBWT + cn->outG), *outOffList = inOffList + cn->inN ;
      printf ("\n    in:") ;
      for (i=0; i<cn->inN; ++i) printf(" %d(%d) [%d]", inList[i].sync,inOffList[i],inList[i].count) ;
      printf ("\n    out:") ;
      for (i=0; i<cn->outN; ++i) printf(" %d(%d) [%d]",outList[i].sync,outOffList[i],outList[i].count) ;
      if (isFM)
	{ FM *inFM = (FM*)inGBWT, *outFM = (FM*)outGBWT ;
	  printf ("\n    inFM:") ;
	  for (i=0; i<cn->inG; ++i)  printf(" %d [%d,%d]",fmIndex(inFM+i),inFM[i].cum,fmOcc(inFM+i)) ;
	  printf ("\n    outFM:") ;
	  for (i=0; i<cn->outG; ++i) printf(" %d [%d,%d]",fmIndex(outFM+i),outFM[i].cum,fmOcc(outFM+i)) ;
	  putchar ('\n') ;
	}
      else
	{ printf ("\n    inG:") ;
	  for (i=0; i<cn->inG; ++i)  printf(" %d [%d]", inGBWT[i].sync, inGBWT[i].count) ;
	  printf ("\n    outG:") ;
	  for (i=0; i<cn->outG; ++i) printf(" %d [%d]", outGBWT[i].sync, outGBWT[i].count) ;
	  putchar ('\n') ;
	}
    }
  else die ("\n  nodePrint: unrecognisable node type") ;
  fflush (stdout) ;
}

static void syncPrint (SyngBWT *sb, I32 k)
{
  printf ("node %d", k) ;
  if (k < 0) k = -k ;
  Node *n = arrayp(sb->node, k, Node) ;
  U8   *s = arrayp(sb->status, k, U8) ;
  nodePrint (n, s, sb->isFM) ;
}

// #define TRACE_NODE 352

static I32 syngBWTadd (SyngBWT *sb, I32 k, I32 in, I32 inOff, I32 j, I32 out, I32 outOff) 
{
  if (sb->isFM) die ("syngBWTadd() can't operate in FM mode") ;
    
#ifdef DEBUG_ADD  
  printf ("adding %d %d in %d %d out %d %d - ", k, j, in, inOff, out, outOff) ;
#endif
  bool isPositive = true ;
  if (k < 0)
    { isPositive = false ;
      k = -k ;
      I32 t = in ; in = -out ; out = -t ; // NB in this case we will come from -out to -k to -in
      t = inOff ; inOff = outOff ; outOff = t ;
    }

  Node *n = arrayp(sb->node, k, Node) ;
  U8   *s = arrayp(sb->status, k, U8) ;
#ifdef TRACE_NODE  
  if (k == TRACE_NODE)
    { if (!*s) putchar ('\n') ;
      printf ("adding %d %c %d in %d %d out %d %d", k, isPositive?'+':'-', j, in, inOff, out, outOff) ;
    }
#endif
  if (!*s)               // make a new simple node
    { assert (j == 0) ;
      *n = nodeCreate (isPositive ? k : -k, in, inOff, out, outOff) ;
      *s = NODE_SIMPLE ;
#ifdef DEBUG_ADD  
      printf ("new simple node - jNext 0\n") ;
#endif
#ifdef TRACE_NODE  
      if (k == TRACE_NODE) { printf (" - new simple node") ; nodePrint (n, s, sb->isFM) ; }
#endif
      return 0 ;
    }
  if (*s & NODE_SIMPLE)    // try to update an old simple node
    { SimpleNode *sn = &(n->simple) ;
      if (in == sn->in && inOff == sn->inOff && out == sn->out && outOff == sn->outOff)
	{ if (isPositive)
	    { ++sn->inCount ;
#ifdef DEBUG_ADD  
	      printf (" adding to simple node inCount %d - ", sn->inCount) ;
#endif
	    }
	  else
	    { ++sn->outCount ;
#ifdef DEBUG_ADD  
	      printf (" adding to simple node outCount %d - ", sn->outCount) ;
#endif
	    }
#ifdef DEBUG_ADD
	  printf (" jNext %d\n", j) ;
#endif
#ifdef TRACE_NODE
	  if (k == TRACE_NODE) { printf (" - update simple node") ; nodePrint (n, s, sb->isFM) ; }
#endif
	  return j ;
	}
      else                  // expand out to Complex then fall through to code below
	{ SimpleNode sn = n->simple ; // make a copy this time
	  ComplexNode *cn = &n->complex ;
	  cn->inN = cn->outN = 1 ;
	  cn->inG = cn->outG = 1 ;
	  SyncCount *sc = cn->sc = new (5, SyncCount) ;
	  sc->sync = sn.in ; sc->count = sn.inCount ; ++sc ;    // inList
	  sc->sync = sn.out ; sc->count = sn.outCount ; ++sc ;  // outList
	  sc->sync = 0 ; sc->count = sn.outCount ; ++sc ;        // inGBWT
	  sc->sync = 0 ; sc->count = sn.inCount ; ++sc ;         // outGBWT
	  I32* soff = (I32*)sc ; *soff = sn.inOff ; ++soff ;     // inOffList
	  *soff = sn.outOff ;                                    // outOffList
	  *s = NODE_COMPLEX ;
#ifdef DEBUG_ADD  
	  printf ("unpacking simple node - ") ;
#endif
#ifdef TRACE_NODE  
	  if (k == TRACE_NODE) { printf (" - unpack simple node") ; nodePrint (n,s, sb->isFM) ; }
#endif
	}
    }
  else if (*s & NODE_PACKED)
    { nodeUnpack (n, s) ;
#ifdef DEBUG_ADD  
      printf ("unpacking packed node - ") ;
#endif
    }
  
  // now node must be complex
  ComplexNode *cn = &(n->complex) ;

  // strategy is to conceptually separate out the SyncCount lists as below
  // update them in place if we can
  // if we need to extend one or more of them scSplit() will reassign them and set cn->sc to 0
  SyncCount *inList = cn->sc, *outList = inList + cn->inN ;
  SyncCount *inGBWT = outList + cn->outN, *outGBWT = inGBWT + cn->inG ;
  I32 *inOffList = (I32*) (outGBWT + cn->outG), *outOffList = inOffList + cn->inN ;

  // next find in within inList, or add to the end of the list if necessary
  int inK = 0, inPre = 0 ; // inK is index of in within inList, and inPre is sum(count) before inK
  while (inK < cn->inN && (inList[inK].sync != in || inOffList[inK] != inOff))
    inPre += inList[inK++].count ;
  if (inK == cn->inN)      // add in to inList
    { scSplit (cn, &inList, &outList, &inGBWT, &outGBWT, &inOffList, &outOffList) ;
      inList[inK].sync = in ; inList[inK].count = 0 ; inOffList[inK] = inOff ;
      ++cn->inN ;
#ifdef DEBUG_ADD  
      printf ("adding %d to inList - ", in) ;
#endif
    }
  
  // now do the same for out
  int outK = 0, outPre = 0 ;
  while (outK < cn->outN && (outList[outK].sync != out || outOffList[outK] != outOff))
    outPre += outList[outK++].count ;
  if (outK == cn->outN) // add a sync to outList
    { scSplit (cn, &inList, &outList, &inGBWT, &outGBWT, &inOffList, &outOffList) ;
      outList[outK].sync = out ; outList[outK].count = 0 ; outOffList[outK] = outOff ;
      ++cn->outN ;
#ifdef DEBUG_ADD  
      printf ("adding %d to outList - ", out) ;
#endif
    }

  // update the list count, only for the side we are actually coming from
  if (isPositive) ++inList[inK].count ; else ++outList[outK].count ;

  // now we need to update the GBWT
  SyncCount *g0 = isPositive ? outGBWT : inGBWT, *g = g0 ;  // start of the GBWT to insert into
  int target = isPositive ? outK : inK ;                    // the GBWT sync value to match
  int pre = j + (isPositive ? inPre : outPre) ;             // number of elts before we insert
  I32 jNext = 0 ;                                           // next value of j to return
  int sumCount = 0 ;
  while (sumCount + g->count < pre)                         // find GBWT block for j
    { sumCount += g->count ;
      if (g->sync == target) jNext += g->count ;
      ++g ;
    }
  if (g-g0 >= (isPositive ? cn->outG : cn->inG))
    { nodePrint (n, s, sb->isFM) ;
      die ("addNode %c outG %d inG %d g-g0 %d target %d pre %d j %d jNext %d sumCount %d g[-1].count %d",
	 isPositive?'+':'-',cn->outG,cn->inG,(int)(g-g0), target, pre, j, jNext, sumCount, g[-1].count) ;
    }
  if (g->sync == target)                                    // within matching block - easy
    { jNext -= sumCount - pre ;
      g->count++ ;
#ifdef DEBUG_ADD  
      printf ("increment block %d - ", (int)(g-g0)) ;
#endif
    }
  else if (sumCount + g->count == pre)                      // at the end of this GBWT block
    { if (++g - g0 < (isPositive?cn->outG:cn->inG) && g->sync == target) // can increment next block
	{ g->count++ ;                                      // jNext is already correct
#ifdef DEBUG_ADD  
	  printf ("increment %s-block %d - ", isPositive?"out":"in", (int)(g-g0)) ;
#endif
	}
      else                                                  // add one row to the GBWT
	{ scSplit (cn, &inList, &outList, &inGBWT, &outGBWT, &inOffList, &outOffList) ;
	  if (isPositive)                                   // remember we incremented g!
	    { int z ; for (z = cn->outG ; z > g - g0 ; --z) outGBWT[z] = outGBWT[z-1] ;
	      outGBWT[z].sync = target ; outGBWT[z].count = 1 ; // again jNext is correct
#ifdef DEBUG_ADD  
	      printf ("adding out-block %d - ", z) ;
#endif
	      ++cn->outG ;
	    }
	  else
	    { int z ; for (z = cn->inG ; z > g - g0 ; --z) inGBWT[z] = inGBWT[z-1] ;
	      inGBWT[z].sync = target ; inGBWT[z].count = 1 ;
#ifdef DEBUG_ADD  
	      printf ("adding in-block %d - ", z) ;
#endif
	      ++cn->inG ;
	    }
	}
    }
  else if (pre == 0) // at very start - like above case but don't increment g
    { scSplit (cn, &inList, &outList, &inGBWT, &outGBWT, &inOffList, &outOffList) ;
      if (isPositive)
	{ int z ; for (z = cn->outG ; z > g - g0 ; --z) outGBWT[z] = outGBWT[z-1] ;
	  outGBWT[z].sync = target ; outGBWT[z].count = 1 ; // again jNext is correct
#ifdef DEBUG_ADD  
	  printf ("adding out-block %d - ", z) ;
#endif
	  ++cn->outG ;
	}
      else
	{ int z ; for (z = cn->inG ; z > g - g0 ; --z) inGBWT[z] = inGBWT[z-1] ;
	  inGBWT[z].sync = target ; inGBWT[z].count = 1 ;
#ifdef DEBUG_ADD  
	  printf ("adding in-block %d - ", z) ;
#endif
	  ++cn->inG ;
	}
    }
  else // we are in the middle of a GBWT block - we need to add two rows: new one and reversal
    { scSplit (cn, &inList, &outList, &inGBWT, &outGBWT, &inOffList, &outOffList) ;
      if (isPositive)
	{ int z ; for (z = cn->outG+1 ; z > g - g0 + 2 ; --z) outGBWT[z] = outGBWT[z-2] ;
	  outGBWT[z].sync = outGBWT[z-2].sync ; outGBWT[z].count = sumCount + g->count - pre ; --z ;
	  outGBWT[z].sync = target ; outGBWT[z].count = 1 ; --z ; // again jNext is correct
	  outGBWT[z].count -= (sumCount + g->count - pre) ;
#ifdef DEBUG_ADD  
	  printf ("adding 2 out-blocks %d,%d - ", z, z+1) ;
#endif
	  cn->outG += 2 ;
	}
      else
	{ int z ; for (z = cn->inG+1 ; z > g - g0 + 2 ; --z) inGBWT[z] = inGBWT[z-2] ;
	  inGBWT[z].sync = inGBWT[z-2].sync ; inGBWT[z].count = sumCount + g->count - pre ; --z ;
	  inGBWT[z].sync = target ; inGBWT[z].count = 1 ; --z ;
	  inGBWT[z].count -= (sumCount + g->count - pre) ;
#ifdef DEBUG_ADD  
	  printf ("adding 2 in-blocks %d,%d - ", z, z+1) ;
#endif
	  cn->inG += 2 ;
	}
    }

#ifdef DEBUG_ADD
  int z ;
  if (cn->inG) { printf ("inG ") ; for (z = 0 ; z < cn->inG ; z++) printf ("(%d,%d)", inGBWT[z].sync, inGBWT[z].count) ; }
  if (cn->outG) { printf ("outG ") ; for (z = 0 ; z < cn->outG ; z++) printf ("(%d,%d)", outGBWT[z].sync, outGBWT[z].count) ; }
  for (z = 0 ; z < cn->inG ; z++) if (inGBWT[z].count < 0) die ("count < 0") ;
  for (z = 0 ; z < cn->outG ; z++) if (outGBWT[z].count < 0) die ("count < 0") ;
#endif

  // we are done
  if (!cn->sc) // it was split - we have to reform it
    { cn->sc = new (cn->inN + cn->outN + cn->inG + cn->outG + (cn->inN + cn->outN + 1)/2, SyncCount) ;
      SyncCount *sc = cn->sc ;
      memcpy (sc, inList, cn->inN*sizeof(SyncCount)) ;   sc += cn->inN ;
      memcpy (sc, outList, cn->outN*sizeof(SyncCount)) ; sc += cn->outN ;
      memcpy (sc, inGBWT, cn->inG*sizeof(SyncCount)) ;   sc += cn->inG ;
      memcpy (sc, outGBWT, cn->outG*sizeof(SyncCount)) ; sc += cn->outG ;
      I32 *ic = (I32*)sc ;
      memcpy (ic, inOffList, cn->inN*sizeof(I32)) ;      ic += cn->inN ;
      memcpy (ic, outOffList, cn->outN*sizeof(I32)) ;
      // call free() here because we already allowed for removal on creation
      free (inList) ; free (outList) ; free (inGBWT) ; free (outGBWT) ;
      free (inOffList) ; free(outOffList) ;
    }
  
#ifdef DEBUG_ADD  
  printf (" - jNext %d\n", jNext) ;
#endif
#ifdef TRACE_NODE  
  if (k == TRACE_NODE) { printf ("  jNext %d", jNext) ; nodePrint (n, s, sb->isFM) ; }
#endif

  return jNext ;
}

static I32 syngBWTnext (SyngBWT *sb, I32 k, I32 in, I32 inOff, I32 j, I32 *out, I32 *outOff)
{
  *out = 0 ; // default, indicating failure to match or extend

  bool isPositive = true ;
  if (k < 0)
    { isPositive = false ;
      k = -k ;
    } // don't swap in and out here

  if (k >= arrayMax(sb->node))
    die ("syngBWTnext: k %d >= arrayMax(sb->node) %lld", k, arrayMax(sb->node)) ;
  Node *n = arrayp(sb->node, k, Node) ;
  U8   *s = arrayp(sb->status, k, U8) ;
  if (!*s) die ("syngBWTnext %c %d hit an empty node", isPositive?'+':'-', k) ;
  if (*s & NODE_SIMPLE)
    { SimpleNode *sn = &(n->simple) ;
      if (isPositive && in == sn->in && inOff == sn->inOff && j < sn->inCount)
	{ *out = sn->out ; *outOff = sn->outOff ; return j ; }
      else if (!isPositive && in == -sn->out && inOff == sn->outOff && j < sn->outCount)
	{ *out = -sn->in ; *outOff = sn->inOff ; return j ; }
      else
	die ("syngBWTnext %c %d failed to match %d to %d in SIMPLE",
	     isPositive?'+':'-', k, in, isPositive?sn->in:-sn->out) ;
    }
  else if (*s & NODE_PACKED)
    nodeUnpack (n, s) ;
  
  // now node must complex
  ComplexNode *cn = &(n->complex) ;

  SyncCount *inList = cn->sc, *outList = inList + cn->inN ;
  SyncCount *inGBWT = outList + cn->outN, *outGBWT = inGBWT + cn->inG ;
  FM  *inFM = (FM*)inGBWT, *outFM = (FM*)outGBWT ;
  I32 *inOffList = (I32*) (outGBWT + cn->outG), *outOffList = inOffList + cn->inN ;
  I32  jNext ; // will hold return value

  if (isPositive) // find in within inList
    { int x = 0 ; // x is index of in within inList
      while (x < cn->inN && (inList[x].sync != in || inOffList[x] != inOff))
	j += inList[x++].count ;
      if (x == cn->inN) die ("syngBWTnext +%d failed to find %d in inList", k, in) ;

      if (sb->isFM) 
	{ int i0 = 0, i1 = cn->outG, i ;
	  while (i1 > i0+1) // binary search
	    { i = (i0+i1)/2 ;
	      if (outFM[i].cum <= j) i0 = i ;
	      else i1 = i ;
	    }
	  int index = fmIndex(outFM+i0) ;
	  *out = outList[index].sync ; *outOff = outOffList[index] ;
	  jNext = j - outFM[i0].cum + fmOcc(outFM+i0) ;
	}
      else
	{ SyncCount *g0 = outGBWT, *g = g0 ;  // start of the GBWT to search into
	  I32 *counts = new0 (cn->outN, I32) ;
	  int jStart = j ;
	  while (g->count <= j)                // find GBWT block for j
	    { j -= g->count ;
	      counts[g->sync] += g->count ;
	      ++g ;
	    }
	  if (g-g0 >= cn->outG) die ("error in nextNode") ;
	  *out = outList[g->sync].sync ;
	  *outOff = outOffList[g->sync] ;
	  jNext = j + counts[g->sync] ;
	  newFree (counts, cn->outN, I32) ;
	}
    }
  else
    { int x = 0 ; // x is index of in within outList
      while (x < cn->outN && (outList[x].sync != -in || outOffList[x] != inOff))
	j += outList[x++].count ;
      if (x == cn->outN) die ("syngBWTnext -%d failed to find %d in outList", k, in) ;

      if (sb->isFM)
	{ int i0 = 0, i1 = cn->inG, i ;
	  while (i1 > i0+1) // binary search
	    { i = (i0+i1)/2 ;
	      if (inFM[i].cum <= j) i0 = i ;
	      else i1 = i ;
	    }
	  int index = fmIndex(inFM+i0) ;
	  *out = -inList[index].sync ; *outOff = inOffList[index] ;
	  jNext = j - inFM[i0].cum + fmOcc(inFM+i0) ;
	}
      else
	{ SyncCount *g0 = inGBWT, *g = g0 ;  // start of the GBWT to insert into
	  I32 *counts = new0(cn->inN,I32) ;
	  int jStart = j ;
	  while (g->count <= j)                // find GBWT block for j
	    { j -= g->count ;
	      counts[g->sync] += g->count ;
	      ++g ;
	    }
	  if (g-g0 >= cn->inG) die ("error in nextNode") ;
	  *out = -inList[g->sync].sync ;
	  *outOff = inOffList[g->sync] ;
	  jNext = j + counts[g->sync] ;
	  newFree (counts, cn->outN, I32) ;
	}
    }

  return jNext ;
}

/*******************************************************/

#define CHECK_FM

void syngBWTtoFM (SyngBWT *sb) // convert BWT information to FM mode
{
  assert (sizeof(FM) == sizeof(SyncCount)) ; // we need this!
  if (sb->isFM) return ; // don't try to do it again!
  sb->isFM = true ;
  int i, j ;
  Array occ = arrayCreate (1024, I32) ;
  for (i = 1 ; i < arrayMax (sb->node) ; ++i)
    { Node n = arr(sb->node, i, Node) ;
      U8   s = arr(sb->status, i, U8) ;
      if (!s || (s & NODE_SIMPLE)) continue ; // nothing needed
      if (s & NODE_PACKED) nodeUnpack (&n, &s) ;
      ComplexNode *cn = &(n.complex) ;
 #ifdef CHECK_FM     
      if (cn->inN >= (1 << 16) || cn->outN >= (1 << 16))
	die ("inN %d %x outN %d %x", cn->inN, cn->inN, cn->outN, cn->outN) ;
 #endif

      SyncCount *inList = cn->sc, *outList = inList + cn->inN ;
      SyncCount *inGBWT = outList + cn->outN, *outGBWT = inGBWT + cn->inG ;
      FM        *inFM = (FM*) inGBWT, *outFM = (FM*) outGBWT ;  // will replace in place

      U32 cum = 0 ; array(occ,cn->inN-1,I32) = 0 ; memset (arrp(occ,0,I32), 0, cn->inN*sizeof(I32)) ;
      for (j = 0 ; j < cn->inG ; ++j, ++inGBWT, ++inFM)
	{ U32 k = inGBWT->sync, bit = inGBWT->count ;
	  if (k < 64) { inFM->is2 = false ; inFM->fm.index1 = k ; inFM->fm.occ1 = arr(occ,k,I32) ; }
	  else        { inFM->is2 = true  ; inFM->fm.index2 = k ; inFM->fm.occ2 = arr(occ,k,I32) ; }
	  inFM->cum = cum ; cum += bit ; arr(occ,k,I32) += bit ;
	}
 #ifdef CHECK_FM     
      for (j = 0 ; j < 64 && j < cn->inN ; ++j)
	if (arr(occ,j,I32) >= (1 << 26)) die ("occ[%d] %d %x", j, arr(occ,j,I32), arr(occ,j,I32)) ;
      for (j = 64 ; j < cn->inN ; ++j)
	if (arr(occ,j,I32) >= (1 << 16)) die ("occ[%d] %d %x", j, arr(occ,j,I32), arr(occ,j,I32)) ;
 #endif

      cum = 0 ; array(occ,cn->outN-1,I32) = 0 ; memset (arrp(occ,0,I32), 0, cn->outN*sizeof(I32)) ;
      for (j = 0 ; j < cn->outG ; ++j, ++outGBWT, ++outFM)
	{ U32 k = outGBWT->sync, bit = outGBWT->count ;
	  if (k < 64) { outFM->is2 = false ; outFM->fm.index1 = k ; outFM->fm.occ1 = arr(occ,k,I32) ; }
	  else        { outFM->is2 = true ;  outFM->fm.index2 = k ; outFM->fm.occ2 = arr(occ,k,I32) ; }
	  outFM->cum = cum ; cum += bit ; arr(occ,k,I32) += bit ;
	}
 #ifdef CHECK_FM     
      for (j = 0 ; j < 64 && j < cn->outN ; ++j)
	if (arr(occ,j,I32) >= (1 << 26)) die ("occ[%d] %d %x", j, arr(occ,j,I32), arr(occ,j,I32)) ;
      for (j = 64 ; j < cn->outN ; ++j)
	if (arr(occ,j,I32) >= (1 << 16)) die ("occ[%d] %d %x", j, arr(occ,j,I32), arr(occ,j,I32)) ;
 #endif
    }
  arrayDestroy (occ) ;
}

/*******************************************************/
/************* external interfaces *********************/

static SyngBWTpath *pathCreate (SyngBWT *sb, I32 startNode)
{
  SyngBWTpath *sbp = new0 (1, SyngBWTpath) ;
  sbp->sb = sb ;
  sbp->lastNode = 0 ;
  sbp->lastOff = 0 ;
  sbp->thisNode = startNode ;
  return sbp ;
}

void syngBWTpathDestroy (SyngBWTpath *sbp) { newFree (sbp, 1, SyngBWTpath) ; }

/************* interface to add a path *******************/

SyngBWTpath *syngBWTpathStartNew (SyngBWT *sb, I32 startNode)
{
  SyngBWTpath *sbp = pathCreate (sb, startNode) ;
  sbp->jLast = startCount (sb, startNode, true) ;
  return sbp ;
}

void syngBWTpathAdd (SyngBWTpath *sbp, I32 nextNode, I32 offset)
{
  sbp->jLast = syngBWTadd (sbp->sb, sbp->thisNode, sbp->lastNode, sbp->lastOff, sbp->jLast,
			     nextNode, offset) ;
  sbp->lastNode = sbp->thisNode ;
  sbp->lastOff = offset ;
  sbp->thisNode = nextNode ;
}

void syngBWTpathFinish (SyngBWTpath *sbp)
{ syngBWTadd (sbp->sb, sbp->thisNode, sbp->lastNode, sbp->lastOff, sbp->jLast, 0, 0) ;
  syngBWTpathDestroy (sbp) ;
}

/************* interface to follow an existing path *******************/

SyngBWTpath *syngBWTpathStartOld (SyngBWT *sb, I32 startNode, I32 count)
{
  SyngBWTpath *sbp = pathCreate (sb, startNode) ;
  sbp->jLast = count ;
  if (count >= startCount (sb, startNode, false))
    die ("syngBWTpathStartOld startNode %d count %d >= startCount %d",
	 startNode, count, startCount(sb,startNode,false)) ;
  return sbp ;
}

bool syngBWTpathNext (SyngBWTpath *sbp, I32 *nextNode, I32 *offset)
{
  int oldjLast = sbp->jLast ;
  sbp->jLast = syngBWTnext (sbp->sb, sbp->thisNode, sbp->lastNode, sbp->lastOff, sbp->jLast,
			      nextNode, offset) ;
  if (*nextNode)
    { sbp->lastNode = sbp->thisNode ;
      sbp->lastOff = *offset ;
      sbp->thisNode = *nextNode ;
      return true ;
    }
  else
    return false ;
}

/****************************************************************/
/********************* write the SyngBWT ************************/

void syngBWTwrite (OneFile *of, SyngBWT *sb)
{
  int i, j, bufSize = 0 ;
  I64 *buf = 0 ; 
  
  for (i = 1 ; i < arrayMax(sb->node) ; ++i)
    { Node n = arr(sb->node, i, Node) ;
      U8   s = arr(sb->status, i, U8) ;
      oneInt(of,0) = sb->fixedLen ; oneWriteLine (of, 'V', 0, 0) ;
      if (!s) continue ; // must come after we have written the node
      if (s & NODE_SIMPLE)
	{ SimpleNode *sn = &(n.simple) ;
	  oneInt(of,0) = sn->in ; oneInt(of,1) = sn->inOff ; oneInt(of,2) = sn->inCount ;
	  oneWriteLine(of, 'E', 0, 0) ; // + edge
	  oneInt(of,0) = sn->out ; oneInt(of,1) = sn->outOff ; oneInt(of,2) = sn->outCount ;
	  oneWriteLine(of, 'e', 0, 0) ; // - edge
	}
      else
	{ if (s & NODE_PACKED) nodeUnpack (&n, &s) ;
	  ComplexNode *cn = &(n.complex) ;
	  SyncCount *inList = cn->sc, *outList = inList + cn->inN ;
	  SyncCount *inGBWT = outList + cn->outN, *outGBWT = inGBWT + cn->inG ;
	  I32 *inOffList = (I32*) (outGBWT + cn->outG), *outOffList = inOffList + cn->inN ;
	  for (j = 0 ; j < cn->inN ; ++j) // incoming edges
	    { oneInt(of,0) = inList[j].sync; oneInt(of,1) = inOffList[j] ; 
	      oneInt(of,2) = inList[j].count ; oneWriteLine(of, 'E', 0, 0) ; // + edge
	    }
	  if (cn->outG > 1) // outgoing GBWT, which goes from the incoming edges
	    { if (cn->outG > bufSize)
		{ if (bufSize) newFree(buf, bufSize, I64) ;
		  bufSize = 2*cn->outG ; buf = new (bufSize, I64) ;
		}
	      for (j = 0 ; j < cn->outG ; ++j) buf[j] = outGBWT[j].sync ;
	      oneWriteLine (of, 'B', cn->outG, buf) ;
	      for (j = 0 ; j < cn->outG ; ++j) buf[j] = outGBWT[j].count ;
	      oneWriteLine (of, 'C', cn->outG, buf) ;
	    }
	  for (j = 0 ; j < cn->outN ; ++j) // outgoing edges
	    { oneInt(of,0) = outList[j].sync; oneInt(of,1) = outOffList[j] ; 
	      oneInt(of,2) = outList[j].count ; oneWriteLine(of, 'e', 0, 0) ; // - edge
	    }
	  if (cn->inG > 1) // incoming GBWT, which follows from the outgoing edges (in reverse)
	    { if (cn->inG > bufSize)
		{ if (bufSize) newFree(buf, bufSize, I64) ;
		  bufSize = 2*cn->inG ; buf = new (bufSize, I64) ;
		}
	      for (j = 0 ; j < cn->inG ; ++j) buf[j] = inGBWT[j].sync ;
	      oneWriteLine (of, 'b', cn->inG, buf) ;
	      for (j = 0 ; j < cn->inG ; ++j) buf[j] = inGBWT[j].count ;
	      oneWriteLine (of, 'c', cn->inG, buf) ;
	    }
	}
    }
  newFree (buf, bufSize, I64) ;
}

/****************** and read it - thread this for speed *******************/

typedef struct {
  OneFile *of ;
  SyngBWT *sb ;
  U32      v1, vn ; // first V line to read, and last
  I64      eTotal ;
} ReadThread ;

static void *threadRead (void *arg)
{
  ReadThread *rt = (ReadThread*)arg ;
  SyngBWT    *sb = rt->sb ;
  OneFile    *of = rt->of ;
  U32 i, j ;
  I64 EMax, eMax, BMax, bMax ;

  if (!oneStatsContains (of, 'V', 'E', &EMax, 0) ||  // max number of E lines
      !oneStatsContains (of, 'V', 'e', &eMax, 0) ||  // max number of e lines
      !oneStatsContains (of, 'V', 'B', 0, &BMax) ||  // longest sum of B lines
      !oneStatsContains (of, 'V', 'b', 0, &bMax))    // longest sum of b lines
    die ("syngBWTread: failed to find max list lengths for vertices") ;
  SyncCount *eIn = new (EMax, SyncCount), *eOut = new (eMax, SyncCount) ;
  SyncCount *gIn = new (bMax, SyncCount), *gOut = new (BMax, SyncCount) ;
  I32       *inOff = new (EMax, I32),     *outOff = new (eMax, I32) ;

  if (!oneGoto (of, 'V', rt->v1)) die ("failed to locate to V line %u", rt->v1) ;
  oneReadLine (of) ; // read the first V line
  rt->eTotal = 0 ;
  for (i = rt->v1 ; i <= rt->vn ; ++i)
    { if (oneInt(of,0) != sb->fixedLen)
	die ("syngBWTread node length %lld != expected %d", oneInt(of,0), sb->fixedLen) ;
      Node *n = arrayp(sb->node, i, Node) ;
      U8   *s = arrayp(sb->status, i, U8) ;
      int   inN = 0, outN = 0, inG = 0, outG = 0 ;
      I64  *oil ; // used for oneIntList()
      while (oneReadLine (of) && of->lineType != 'V')
	switch (of->lineType)
	  { 
	  case 'E':
	    eIn[inN].sync = oneInt(of,0) ; inOff[inN] = oneInt(of,1) ; eIn[inN].count = oneInt(of,2) ;
	    if (eIn[inN].sync == 0) startCountAdd (sb, i, eIn[inN].count) ; // start node
	    ++inN ; ++rt->eTotal ; break ;
	  case 'e':
	    eOut[outN].sync = oneInt(of,0) ; outOff[outN] = oneInt(of,1) ; eOut[outN].count = oneInt(of,2) ;
	    if (eOut[outN].sync == 0) startCountAdd (sb, -i, eOut[outN].count) ; // start node
	    ++outN ; ++rt->eTotal ; break ;
	  case 'B':
	    oil = oneIntList(of) ;
	    outG = oneLen(of) ; for (j = 0 ; j < outG ; ++j) gOut[j].sync = oil[j] ;
	    break ;
	  case 'b':
	    oil = oneIntList(of) ;
	    inG = oneLen(of) ; for (j = 0 ; j < inG ; ++j) gIn[j].sync = oil[j] ;
	    break ;
	  case 'C':
	    oil = oneIntList(of) ;
	    for (j = 0 ; j < outG ; ++j) gOut[j].count = oil[j] ;
	    break ;
	  case 'c':
	    oil = oneIntList(of) ;
	    for (j = 0 ; j < inG ; ++j) gIn[j].count = oil[j] ;
	    break ;
	  default:
	    die ("unrecognized linetype %c in vertex", of->lineType) ;
	  }
      // now we have read everything we need to construct the node
      if (inN == 0 && outN == 0)	// empty node
	*s = 0 ;
      else if (inN == 1 && outN == 1)	// simple node
	{ *s = NODE_SIMPLE ;
	  SimpleNode *sn = &(n->simple) ;
	  sn->in = eIn[0].sync ; sn->inOff = inOff[0] ; sn->inCount = eIn[0].count ;
	  sn->out = eOut[0].sync ; sn->outOff = outOff[0] ; sn->outCount = eOut[0].count ;
	}
      else				// complex node
	{ *s = NODE_COMPLEX ;
	  ComplexNode *cn = &(n->complex) ;
	  if (!inG) { inG = 1 ; gIn->sync = 0 ; assert (inN == 1) ; gIn->count = eIn->count ; }
	  if (!outG) { outG = 1 ; gOut->sync = 0 ; assert (outN == 1) ; gOut->count = eOut->count ; }
	  SyncCount *sc = cn->sc = new (inN + outN + inG + outG + (inN+outN+1)/2, SyncCount) ;
	  memcpy (sc, eIn, inN*sizeof(SyncCount)) ; cn->inN = inN ; sc += inN ;
	  memcpy (sc, eOut, outN*sizeof(SyncCount)) ; cn->outN = outN ; sc += outN ;
	  memcpy (sc, gIn, inG*sizeof(SyncCount)) ; cn->inG = inG ; sc += inG ;
	  memcpy (sc, gOut, outG*sizeof(SyncCount)) ; cn->outG = outG ; sc += outG ;
	  I32* ip = (I32*)sc ;
	  memcpy (ip, inOff, inN*sizeof(I32)) ; ip += inN ;
	  memcpy (ip, outOff, outN*sizeof(I32)) ;
	}
#ifdef TRACE_NODE
      if (i == TRACE_NODE) { printf ("syngBWTread %d", TRACE_NODE) ; nodePrint (n, s, false) ; }
#endif
    }
  
  newFree (eIn, EMax, SyncCount) ; newFree (eOut, eMax, SyncCount) ;
  newFree (gIn, bMax, SyncCount) ; newFree (gOut, BMax, SyncCount) ;
  newFree (inOff, EMax, I32) ; newFree (outOff, eMax, I32) ;
  return 0 ;
}


SyngBWT *syngBWTread (OneFile *of)
{
  I64 i, nv ;
  if (!oneStats (of, 'V', &nv, 0, 0) || !nv || !oneGoto (of, 'V', 1))
    die ("syngBWTread: no Vertex objects in .1gbwt file or can't locate to them") ;
  
  oneReadLine (of) ;	// we are just before the 1st V line
  int fixedLen = oneInt(of,0) ;
  SyngBWT *sb = syngBWTcreate (fixedLen, nv+1) ;
  if (!sb) die ("failed to create syngBWT of size %lld", nv+1) ;

  int         nThread = of->share ; if (!nThread) nThread = 1 ;
  pthread_t  *threads = new (nThread, pthread_t) ;
  ReadThread *rt = new0 (nThread, ReadThread) ;
  for (i = 0 ; i < nThread ; ++i)
    { rt[i].sb = sb ; rt[i].of = of + i ;
      rt[i].v1 = 1 + (i*nv) / nThread ;
      rt[i].vn = ((i+1)*nv) / nThread ;
      pthread_create (&threads[i], 0, threadRead, &rt[i]) ;
    }
  I64 eTotal = 0 ;
  for (i = 0 ; i < nThread ; ++i)
    { pthread_join (threads[i], 0) ; // wait for the threads to complete
      eTotal += rt[i].eTotal ;
    }
  fprintf (stdout, "read GBWT with %lld vertices and %lld edges\n", nv, eTotal) ;

  newFree (threads, nThread, pthread_t) ;
  newFree (rt, nThread, ReadThread) ;
  timeUpdate (stdout) ;
  return sb ;
}

/*****************************************************************************/

void syngBWTstat (SyngBWT *sb)
{
  int N = arrayMax(sb->node) ;
  int nSimple = 0, maxTotal = 0, nBlank ;
  Array eHist = arrayCreate (1024, int) ;
  Array cHist = arrayCreate (1024, int) ;
  Array gHist = arrayCreate (1024, int) ;
  int i, j ;
  I64 eSum = 0, eTotal = 0 ;
  for (i = 1 ; i < N ; ++i)
    { Node n = arr(sb->node, i, Node) ;
      U8   s = arr(sb->status, i, U8) ;
      if (!s) { ++nBlank ; continue ; }
      if (s & NODE_SIMPLE) ++nSimple ;
      else
	{ if (s & NODE_PACKED) nodeUnpack (&n, &s) ;
	  ComplexNode *cn = &(n.complex) ;
	  SyncCount *inList = cn->sc, *outList = inList + cn->inN ;
	  SyncCount *inGBWT = outList + cn->outN, *outGBWT = inGBWT + cn->inG ;

	  int total = 0 ;
	  ++array(eHist,cn->inN,int) ;
	  ++array(gHist,cn->inG,int) ;
	  for (j = 0 ; j < cn->inN ; ++j)
	    { ++array(cHist,inList[j].count,int) ;
	      total += inList[j].count ;
	      eSum += inList[j].count ;
	      eTotal += inList[j].count * (j+1) ;
	    }
	  if (total > maxTotal) maxTotal = total ;

	  total = 0 ;
	  ++array(eHist,cn->outN,int) ;
	  ++array(gHist,cn->outG,int) ;
	  for (j = 0 ; j < cn->outN ; ++j)
	    { ++array(cHist,outList[j].count,int) ;
	      total += outList[j].count ;
	      eSum += outList[j].count ;
	      eTotal += outList[j].count * (j+1) ;
	    }
	  if (total > maxTotal) maxTotal = total ;
	}
    }
  printf ("  %d nodes of which %d are blank, %d are simple, %d complex %.1f%%\n",
	  N, nBlank, nSimple, N - nBlank - nSimple, (N-nBlank-nSimple)/(0.01*N)) ;
  printf ("    max in or out edge list %d\n", (int)arrayMax(eHist)) ;
  printf ("    max count of any one edge %d\n", (int)arrayMax(cHist)) ;
  printf ("    max total paths through a node %d\n", maxTotal) ;
  printf ("    expected list search time %.2f\n", eTotal / (1.0*eSum)) ;
  printf ("    max size of GBWT\n", (int)arrayMax(gHist)) ;
  int bit = 0 ;
  printf ("edge list dbn\n") ;
  for (i = 1 ; i < 20 ; ++i)
    { printf ("    %2d %8d\n", i, arr(eHist, i, int)) ; bit += arr(eHist, i, int) ; }
  int left = 0 ; while (i < arrayMax(eHist)) left += arr(eHist, i++, int) ; 
  printf     ("   >20 %8d %8d\n", 2*(N-nBlank-nSimple)-bit, left) ;
  arrayDestroy (eHist) ; arrayDestroy (cHist) ;
}

/**************** basic node creation, packing and unpacking ******************/

static Node nodeCreate (I32 k, I32 in, I32 inOff, I32 out, I32 outOff)
{
  Node node ;
  SimpleNode *sn = &(node.simple) ;
  sn->in = in ; sn->inOff = inOff ;
  sn->out = out ; sn->outOff = outOff ;
  if (k > 0) { sn->inCount = 1 ; sn->outCount = 0 ; }
  else { sn->outCount = 1 ; sn->inCount = 0 ; }
  return node ;
}

static I32 nodePack (Node *n, U8 *s) // returns new size in bytes
{
  if (!*s) return 0 ; // nothing in this node
  if (*s & NODE_SIMPLE)
    return sizeof(SimpleNode) ;
  else if (*s & NODE_PACKED)
    return (I32) n->packed.size ;
  else // complex and edited
    { ComplexNode *cn = &(n->complex) ;
      int scSize = cn->inN + cn->outN + cn->inG + cn->outG + (cn->inN + cn->outN + 1)/2 ;
      int maxSize = sizeof(I32) + 5*(4 + scSize*2) ;
      U8 *u = new(maxSize,U8) ;
      U8 *u0 = u ;
      u += intPut (u, cn->inN) ;
      u += intPut (u, cn->outN) ;
      u += intPut (u, cn->inG) ;
      u += intPut (u, cn->outG) ;
      int i ;
      SyncCount *sc = cn->sc ;
      for (i = 0 ; i < scSize ; ++i, ++sc)
	{ u += intPut (u, sc->sync) ; u += intPut (u, sc->count) ; }
      newFree (cn->sc, scSize, SyncCount) ;
      PackedNode *pn = &(n->packed) ;
      pn->size = u - u0 ;
      pn->data = new(pn->size,U8) ;
      memcpy(pn->data,u0,(size_t)pn->size) ;
      newFree (u0, maxSize, U8) ;
      *s &= ~NODE_COMPLEX ;
      *s |= NODE_PACKED ;
      return pn->size ;
    }
}

static I32 nodeUnpack (Node *n, U8 *s) // returns new size in bytes (beyond the node)
{
  if (!*s) return 0 ;   // nothing in this node
  if (*s & NODE_SIMPLE) // don't do anything
    return 0 ;
  else if (*s & NODE_COMPLEX) // also don't do anything - must calculate size
    { ComplexNode *cn = &(n->complex) ;
      int scSize = cn->inN + cn->outN + cn->inG + cn->outG + (cn->inN + cn->outN + 1) / 2 ;
      return scSize * sizeof(SyncCount) ;
    }
  else if (*s & NODE_PACKED) // packed
    { I64 packedSize = n->packed.size ;
      U8 *u = n->packed.data, *u0 = u ;
      ComplexNode *cn = &(n->complex) ;
      *s &= ~NODE_PACKED ;
      *s |= NODE_COMPLEX ;
      u += intGet (u, &cn->inN) ;
      u += intGet (u, &cn->outN) ;
      u += cn->inN + cn->outN ; // offsets
      u += intGet (u, &cn->inG) ;
      u += intGet (u, &cn->outG) ;
      int scSize = cn->inN + cn->outN + cn->inG + cn->outG + (cn->inN + cn->outN + 1) / 2 ;
      SyncCount *sc = cn->sc = new (scSize, SyncCount) ;
      int i ;
      for (i = 0 ; i < scSize ; ++i, ++sc)
	{ u += intGet (u, &sc->sync) ; u += intGet (u, &sc->count) ; }
      newFree (u0, packedSize, U8) ;
      return scSize * sizeof(SyncCount) ;
    }
  else
    { die ("something wrong in nodeUnpack") ; return 0 ; } // need return 0 for compiler happiness
}

/************ low level integer packing/unpacking routines *************/

static inline int intGet (unsigned char *u, I32 *pval)
{
  switch (u[0] >> 5)
    {
    case 2: case 3: // single byte positive
      *pval = (I32) (u[0] & 0x3f) ; return 1 ;
    case 6: case 7: // single byte negative
      *pval =  (I32) u[0] | 0xffffff00 ; return 1 ;
    case 1: // two bytes positive
      *pval = (I32) (u[0] & 0x1f) << 8 | (I32)u[1] ; return 2 ;
      //     *pval = - ((I32) (u[0] & 0x1f) << 8 | (I32)u[1]) ; return 2 ;
    case 0:
      switch (u[0] & 0x07)
	{
	case 0: die ("int packing error") ; break ;
	case 1: *pval = *(I32*)(u+1) & 0x0000ffff ; return 3 ;
	case 2: *pval = *(I32*)(u+1) & 0x00ffffff ; return 4 ;
	case 3: *pval = *(I32*)(u+1) ; return 5 ;
	}
      break ;
    case 4:
      switch (u[0] & 0x07)
	{
	case 0: die ("int packing error") ; break ;
	case 1: *pval = *(I32*)(u+1) | 0xffff0000 ; return 3 ;
	case 2: *pval = *(I32*)(u+1) | 0xff000000 ; return 4 ;
	case 3: *pval = *(I32*)(u+1) ; return 5 ;
	}
      break ;
    }
  return 0 ; // shouldn't get here, but needed for compiler happiness
}

static inline int intPut (unsigned char *u, I32 val)
{
  if (val >= 0)
    { if (     !(val & 0xffffffc0)) { *u = val | 0x40 ;  return 1 ; } // up to 63
      else if (!(val & 0xffffe000)) { *u++ = (val >> 8) | 0x20 ; *u = val & 0xff ; return 2 ; } // up to 8191
      else if (!(val & 0xffff0000)) { *u++ = 1 ; *(I32*)u = val ; return 3 ; }
      else if (!(val & 0xff000000)) { *u++ = 2 ; *(I32*)u = val ; return 4 ; }
      else                          { *u++ = 3 ; *(I32*)u = val ; return 5 ; }
    }
  else
    { if (     !(~val & 0xffffffc0)) { *u = val | 0x40 ;  return 1 ; }
      //     else if (!(~val & 0xffffe000)) { *u++ = (val >> 8) | 0x20 ; *u = val & 0xff ; return 2 ; }
      else if (!(~val & 0xffff0000)) { *u++ = 0x81 ; *(I32*)u = val ; return 3 ; }
      else if (!(~val & 0xff000000)) { *u++ = 0x82 ; *(I32*)u = val ; return 4 ; }
      else                           { *u++ = 0x83 ; *(I32*)u = val ; return 5 ; }
    }
}

/************** main ****************/

#ifdef TEST

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  --argc ; ++argv ;
  if (argc != 2) die ("usage: syngbwt <XX.1khash> <XX.1syncseq>") ;

  OneFile *of = oneFileOpenRead (*argv, 0, "khash", 1) ;
  if (!of) die ("failed to open khash file %s", *argv) ;
  KmerHash *kh = kmerHashReadOneFile (of) ;
  if (!kh) die ("failed to read khash file %s", *argv) ;
  oneFileClose (of) ;
  timeUpdate (stdout) ;
  
  --argc ; ++argv ;
  of = oneFileOpenRead (*argv, 0, "path", 1) ;
  if (!of) die ("failed to open path file %s", *argv) ;
  I64 max, len, i, *x, total = 0 ;
  oneStats (of, 'P', 0, &max, 0) ;
  I32 *sync = new (max, I32), *off = new (max, I32), nPath = 0 ;
  Array aNode = arrayCreate (max+1, Node) ;
  Array aStatus = arrayCreate (max+1, U8) ;
  while (oneReadLine (of))
    if (of->lineType == 'P') // process it
      { // printf ("adding sequence %lld %lld\n", oneInt(of,0), oneInt(of,1)) ;
      }
    else if (of->lineType == 'z') // path = sync list
      { len = oneLen(of) ; x = oneIntList(of) ; total += len ;
	for (i = 0 ; i < len ; ++i) sync[i] = x[i] ;
      }
    else if (of->lineType == 'd') // directions
      { char *d = oneString (of) ;
	if (oneLen(of) != len) die ("length mismatch D len %lld != %lld", oneLen(of), len) ;
	for (i = 0 ; i < len ; ++i) if (d[i] == '-') sync[i] = -sync[i] ;
      }
    else if (of->lineType == 'o') // offsets (positions) in original read
      { x = oneIntList(of) ;
	if (oneLen(of) != len) die ("length mismatch P len %lld != %lld", oneLen(of), len) ;
	off[0] = x[0] ; for (i = 1 ; i < len ; ++i) off[i] = x[i] - x[i-1] ;
	// now process it
	I32 j = syngBWTadd (sync[0], aNode, aStatus, nSync++, 0, 0, sync[1], off[1]) ;
	sync[len] = 0 ; off[len] = 0 ;
	for (i = 1 ; i < len ; ++i)
	  j = syngBWTadd (sync[i], aNode, aStatus, j, sync[i-1], off[i], sync[i+1], off[i+1]) ;
      }
  oneFileClose (of) ;
  newFree (sync, max, I32) ; newFree (off, max, I32) ;
  printf ("read and processed %d paths with %lld total kmers\n", nSync, total) ;
  timeUpdate (stdout) ;

  // now write out the GBWT
  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
  of = oneFileOpenWriteNew ("out.1gbwt", schema, "gbwt", true, 1) ;
  syngBWTwrite (of, kh->len, aNode, aStatus) ;
  oneFileClose (of) ;
  arrayDestroy (aNode) ; arrayDestroy (aStatus) ;
  kmerHashDestroy (kh) ;
  
  timeTotal (stdout) ;
}

#endif // TEST

#ifdef INTGET_TEST

// this is a test of intPut(), intGet
//   all 4 billion I32 put and get in 26s on my Mac 240930 
  
  I32 i, inCount, outCount, x ;
  U8  buf[8] ;

  for (i = 0 ; i < ((U32)1 << 31) ; ++i)
    { inCount = intPut (buf, i) ;
      outCount = intGet (buf, &x) ;
      if (inCount != outCount || i != x)
	die ("inCount %d outCount %d i %d x %d", inCount, outCount, i, x) ;
      inCount = intPut (buf, -i) ;
      outCount = intGet (buf, &x) ;
      if (inCount != outCount || -i != x)
	die ("inCount %d outCount %d i %d x %d", inCount, outCount, -i, x) ;
    }

#endif // INTGET_TEST

/************** end of file **************/

#ifdef CONVERT_SERIAL

void updateNodePass2 (Node *node, U8 *status, I32 x, I64 off, bool isIn)
{
  if (*status & NODE_SIMPLE)
    { SimpleNode *sn = &(node->simple) ;
      if (isIn)
	{ if (!sn->inOff) { sn->inOff = off ; sn->inCount = 1 ; }
	  else if (sn->inOff == off) ++sn->inCount ;
	  else // convert to complex node
	    { SimpleNode old = *sn ;
	      ComplexNode *cn = &(node->complex) ;
	      *status &= ~NODE_SIMPLE ;
	      *status |= NODE_COMPLEX ;
	      cn->inN = 2 ; cn->outN = 1 ;
	      cn->sc = new0 (cn->inN + cn->outN + 2, SyncCount) ;
	      cn->sc[0].sync = old.in ; cn->sc[1].sync = old.in, cn->sc[2].sync = old.out ;
	      cn->sc[0].count = old.inCount ; cn->sc[1].count = 1 ; cn->sc[2].count = old.outCount ;
	      I32 *ofs = (I32*)&cn->sc[3] ; ofs[0] = old.inOff ; ofs[1] = off ; ofs[2] = old.outOff ;
	    }
	}
      else // isOut
	{ if (!sn->outOff) { sn->outOff = off ; sn->outCount = 1 ; }
	  else if (sn->outOff == off) ++sn->outCount ;
	  else // convert to complex node
	    { SimpleNode old = *sn ;
	      ComplexNode *cn = &(node->complex) ;
	      *status &= ~NODE_SIMPLE ;
	      *status |= NODE_COMPLEX ;
	      cn->inN = 1 ; cn->outN = 2 ;
	      cn->sc = new0 (cn->inN + cn->outN + 2, SyncCount) ;
	      cn->sc[0].sync = old.in ; cn->sc[1].sync = old.out, cn->sc[2].sync = old.out ;
	      cn->sc[0].count = old.inCount ; cn->sc[1].count = old.outCount ; cn->sc[2].count = 1 ;
	      I32 *ofs = (I32*)&cn->sc[3] ; ofs[0] = old.inOff ; ofs[1] = old.outOff ; ofs[2] = off ;
	    }
	}
    }
  else if (*status & NODE_COMPLEX)
    { ComplexNode *cn = &(node->complex) ;
      if (isIn)
	{ SyncCount *scIn = cn->sc ;
	  I32       *scInOff = (I32*)(scIn + cn->inN + cn->outN) ;
	  int i ;
	  for (i = 0 ; i < cn->inN ; ++i)
	    if (scIn[i].sync == x && (!scInOff[i] || scInOff[i] == off)) break ;
	  if (i < cn->inN) // found
	    { ++scIn[i].count ; scInOff[i] = off ; }
	  else // should be rare - have to extend the array
	    { SyncCount *oldSc = cn->sc, *old = oldSc ;
	      scIn = cn->sc = new0 (cn->inN+1 + cn->outN + (cn->inN+cn->outN+2)/2, SyncCount) ;
	      memcpy (scIn, old, cn->inN*sizeof(SyncCount)) ; scIn += cn->inN ; old += cn->inN ;
	      scIn->sync = x ; scIn->count = 1 ; ++scIn ;
	      memcpy (scIn, old, cn->outN*sizeof(SyncCount)) ; scIn += cn->outN ; old += cn->outN ;
	      I32 *iIn = (I32*)scIn, *iold = (I32*)old ;
	      memcpy (iIn, iold, cn->inN*sizeof(I32)) ; iIn += cn->inN ; iold += cn->inN ;
	      *iIn++ = off ;
	      memcpy (iIn, iold, cn->outN*sizeof(I32)) ;
	      newFree (oldSc, cn->inN + cn->outN + (cn->inN+cn->outN+1)/2, SyncCount) ;
	    }
	}
      else // isOut
	{ SyncCount *scOut = cn->sc + cn->inN ;
	  I32       *scOutOff = (I32*)(scIn + cn->inN + cn->outN) + cn->inN ;
	  int i ;
	  for (i = 0 ; i < cn->outN ; ++i)
	    if (scOut[i].sync == x && (!scOutOff[i] || scOutOff[i] == off)) break ;
	  if (i < cn->outN) // found
	    { ++scOut[i].count ; scOutOff[i] = off ; }
	  else // should be rare - have to extend the array
	    { SyncCount *oldSc = cn->sc, *old = oldSc ;
	      scOut = cn->sc = new0 (cn->inN + cn->outN+1 + (cn->inN+cn->outN+2)/2, SyncCount) ;
	      memcpy (scOut, old, (cn->inN+cn->outN)*sizeof(SyncCount)) ;
	      scOut += cn->inN + cn->outN ; old += cn->inN + cn->outN ;
	      scOut->sync = x ; scOut->count = 1 ; ++scOut ;
	      I32 *iOut = (I32*)sc, *iold = (I32*)old ;
	      memcpy (iOut, old, (cn->inN + cn->outN)*sizeof(I32)) ;
	      iOut += cn->inN + cn->outN ; iold += cn->inN + cn->outN ;
	      *iOut++ = off ;
	      newFree (oldSc, cn->inN + cn->outN + (cn->inN+cn->outN+1)/2, SyncCount) ;
	    }
	}
    }
}

typedef struct {
  I32 offset ;
  I32 count ;
} EdgeInfo ;

typedef struct {
  HashKey k ;
  I32 offset, count ;
} DoubleEdge ;

void addEdge (I32 to, I32 from, I32 off, Array count, Hash hash, Array info, Array doubleEdge)
{
  int i ;
  HashKey k = hashInt2 (to, from) ;
  if (hashAdd (hash, k, &i))
    { EdgeInfo *ei = arrayp(info,i,EdgeInfo) ;
      ei->offset = off ;
      ei->count = 1 ;
      ++array(count,to,int) ;
    }
  else
    { EdgeInfo *ei = arrp(info,i,EdgeInfo) ;
      if (ei->offset == off)
	++ei->count ;
      else if (ei->offset > 0) // a new double
	{ DoubleEdge *de = arrayp(doubleEdge, arrayMax(doubleEdge), DoubleEdge) ;
	  de->k = k ;
	  de->offset = ei->offset ;
	  de->count = ei->count ;
	  de = arrayp(doubleEdge, arrayMax(doubleEdge), DoubleEdge) ;
	  de->k = k ;
	  de->offset = off ;
	  de->count = 1 ;
	  ++array(count,to,int) ;
	  ei->offset = -1 ; // marker for this being a double edge
	}
      else // -1
	{ DoubleEdge *de = arrp(doubleEdge, 0, DoubleEdge) ;
	  int j ;
	  for (j = 0 ; j < arrayMax(doubleEdge) ; ++j, ++de)
	    if (de->k == k && de->offset == off) { ++de->count ; break ; }
	  if (j >= arrayMax(doubleEdge)) // no match to offset
	    { de = array(doubleEdge, arrayMax(doubleEdge), DoubleEdge) ; // redo because may be off end
	      de->k = k ;
	      de->offset = off ;
	      de->count = 1 ;
	      ++array(count,to,int) ;
	    }
	}
    }
}

SyngBWT *syngBWTconvert (I64 nNode, I64 nodeLen, OneFile *ofPath) // convert a .1path file to a SyngBWT
{
  int i, j ;
  
  I64 nPath, maxPathLen ;
  if (!ofPath || !ofPath->subType || strcmp(ofPath->subType, "path") ||
      !oneStats (ofPath, 'z', &nPath, &maxPathLen, 0) || !maxPathLen)
    die ("syngBWTconvert failed to find paths in path file") ;
  char *ofPathName = ofPath->fileName ; // for tidiness
  I64  *nodes = new(maxPathLen+1,64) ; // ensure enough space to read nodes[pathLen]
  oneUserBuffer (ofPath, 'z', nodes) ;
  I64  *offsets = new(maxPathLen+1,64) ; // same for offsets
  oneUserBuffer (ofPath, 'o', offsets) ;

  printf ("nNode %lld nodeLen %lld nPath %lld maxPathLen %lld\n", nNode, nodeLen, nPath, maxPathLen) ;

  SyngBWT *sb = syngBWTcreate (nodeLen, nNode + 1) ; // this will be the output

  // first pass through the path file to count edges into and out of each node

  Array edgeIn     = arrayCreate (nNode, int) ;  // count of incoming edges
  Array edgeOut    = arrayCreate (nNode, int) ; // count of outgoing edges
  Hash  edgeHash   = hashCreate (4*nNode) ;     // use this to store edges
  Array edgeInfo   = arrayCreate (4*nNode, EdgeInfo) ;
  Array doubleEdge = arrayCreate (1024, DoubleEdge) ; // use for cases where same edge has >1 offsets
  
  oneGoto (ofPath, 'P', 1) ; // can't fail since we have paths
  while (oneReadLine (ofPath) && ofPath->lineType != 'z') ; // get to path record
  int zCount = 0 ;
  while (ofPath->lineType == 'z') // for each path
    { ++zCount ;
      int pathLen = oneLen(ofPath) ;
      oneReadLine (ofPath) ;
      if (ofPath->lineType != 'o') die ("%s line %d is not o", ofPathName, ofPath->line) ;
      if (oneLen(ofPath) != pathLen) die ("%s line %d o len != z len", ofPathName, ofPath->line) ;
      if (pathLen) // not sure whether I need this, but safe
	{ I32 last = 0, lastOffset = 0 ;
	  for (i = 0 ; i <= pathLen ; ++i)
	    { I32 next = nodes[i] ;
	      I32 off = offsets[i] - lastOffset ;
	      if (next > 0) addEdge (next, last, off, edgeIn, edgeHash, edgeInfo, doubleEdge) ;
	      if (next < 0) addEdge (-next, -last, off, edgeOut, edgeHash, edgeInfo, doubleEdge) ;
	      if (last > 0) addEdge (last, next, off, edgeOut, edgeHash, edgeInfo, doubleEdge) ;
	      if (last < 0) addEdge (-last, -next, off, edgeIn, edgeHash, edgeInfo, doubleEdge) ;
	      last = next ; lastOffset = offsets[i] ;
	    }
	}
      while (oneReadLine (ofPath) && ofPath->lineType != 'z') ;
    }
  printf ("read %d paths: ", zCount) ; timeUpdate (stdout) ;

  // use the edge counts to set up the node objects
  int nSimple = 0, nComplex = 0, nEdge = 0 ;
  int n1 = 0, n2 = 0, n4 = 0, n8 = 0, n16 = 0, n32 = 0 ; // number of bits needed for an edge
  for (i = 1 ; i < nNode ; ++i)
    { int inN = arr(edgeIn,i,int), outN = arr(edgeOut,i,int) ;
      if (inN < 1 || outN < 1) continue ; // printf ("node %d in %d out %d\n", i, inN, outN) ;
      else if (inN == 1 && outN == 1) // make a simple node
	{ arr(sb->status,i,U8) = NODE_SIMPLE ;
	  ++nSimple ;
	}
      else
	{ arr(sb->status,i,U8) = NODE_COMPLEX ;
	  ComplexNode *cn = &(arr(sb->node,i,Node).complex) ;
	  cn->inN = inN ; cn->outN = outN ;
	  nEdge += inN + outN ;
	  int nSC = inN + outN + (inN + outN + 1)/2 ;
	  cn->sc = new0 (nSC, SyncCount) ;
	  ++nComplex ;
	  if (inN == 2) ++n1 ; else if (inN <= 4) ++n2 ;
	  else if (inN <= 16) ++n4 ; else if (inN <= 256) ++n8 ;
	  else if (inN <= 65536) ++n16 ; else ++n32 ;
	  if (outN == 2) ++n1 ; else if (outN <= 4) ++n2 ;
	  else if (outN <= 16) ++n4 ; else if (outN <= 256) ++n8 ;
	  else if (outN <= 65536) ++n16 ; else ++n32 ;
	}
    }
  printf ("BWT convert pass 1\n") ;
  printf ("  %d nodes %d zero %d simple %d complex with %d edges\n",
	  (int)nNode,(int) (nNode - nSimple - nComplex), nSimple, nComplex, nEdge) ;
  printf ("  edge complexity %d 1bit, %d 2bit, %d 4bit, %d 8bit, %d 16bit, %d 32bit\n",
	  n1, n2, n4, n8, n16, n32) ;
  timeUpdate (stdout) ;

  // second pass - count occupancy of each edge and add offsets, splitting if necessary

  oneGoto (ofPath, 'P', 1) ; // can't fail since we have paths
  while (oneReadLine (ofPath) && ofPath->lineType != 'z') ; // get to path record
  while (ofPath->lineType == 'z') // for each path
    { int pathLen = oneLen(ofPath) ;
      oneReadLine (ofPath) ; // read the 'o' line - we checked this worked in pass 1
      if (pathLen) // not sure whether I need this, but safe
	{ I32   last = 0 ;
	  U8   *lastStatus = arrp(sb->status,0,Status) ;
	  Node *lastNode = arrp(sb->node,0,Node) ;
	  for (i = 0 ; i <= pathLen ; ++i)
	    { I32   next = nodes[i] ;
	      U8   *nextStatus = arrp(sb->status,next>0?next:-next,U8) ;
	      Node *nextNode = arrp(sb->node,next>0?next:-next,Node) ;
	      if (next < 0) updateNodePass2 (nextNode, nextStatus, last, offsets[i], true) ;
	      if (next > 0) updateNodePass2 (nextNode, nextStatus, -last, offsets[i], false) ;
	      if (last > 0) updateNodePass2 (lastNode, lastStatus, next, offsets[i], false) ;
	      if (last < 0) updateNodePass2 (lastNode, lastStatus, -next, offsets[i], true) ;
	      last = next ;
	      lastStatus = nextStatus ;
	      lastNode = nextNode ;
	    }
	}
      while (oneReadLine (ofPath) && ofPath->lineType != 'z') ;
    }
  printf ("BWT convert pass 2\n") ;
  printf (" read %d paths\n", zCount) ;
  timeUpdate (stdout) ;
	  
  hashDestroy (edgeHash) ;
  
  arrayDestroy (edgeIn) ;
  arrayDestroy (edgeOut) ;
  newFree (nodes, maxPathLen+1, I64) ;

  return sb ;
}

int main (int argc, char *argv[])
{
  timeUpdate(0) ;
  if (argc != 3) die ("usage: syngbwtconvert <.1khash file> <.1path file>") ;

  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;

  OneFile *ofK = oneFileOpenRead (argv[1], schema, "khash", 1) ;
  if (!ofK || strcmp (ofK->fileType, "khash") || !oneGoto (ofK, 't', 1) || !oneReadLine (ofK))
    die ("syngBWTconvert failed to extract info from khash file") ;
  I64 nNode = oneInt(ofK,0), nodeLen = oneInt(ofK,1) ;
  oneFileClose (ofK) ;
  
  OneFile *ofPath = oneFileOpenRead (argv[2], schema, "path", 1) ;
  SyngBWT *sb = syngBWTconvert (nNode, nodeLen, ofPath) ;
  oneFileClose (ofPath) ;
  
  syngBWTdestroy (sb) ;
  printf ("FINISHED: ") ; timeTotal (stdout) ;
}

#endif
