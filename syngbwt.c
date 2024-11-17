/*  File: syngbwt.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 11 01:51 2024 (rd109)
 * Created: Mon Sep  9 11:34:51 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "syng.h"

typedef struct {
  I32 in, out ;           // syncs
  I32 inCount, outCount ; // inCount is number of paths entering in +ve direction, outCount in -ve
  U8  inOff, outOff ;     // sync offsets
} SimpleNode ;

typedef struct {
  I32 sync ;
  I32 count ;
} SyncCount ;

typedef struct {
  I32        inN, outN ;   // how many distinct syncs coming in, going out
  I32        inG, outG ;   // how big are the run-length encoded GBWTs
  I32        packOff ;     // offset in packed to find the sync offsets
  SyncCount *sc ;          // inList, then outList, then inGBWT, then outGBWT
  U8        *packed ;      // original packed version to replace this by if not edited
} UnpackedNode ;

typedef union {
  SimpleNode   *simple ;
  UnpackedNode *unpacked ;
  U8           *packed ;
} Node ;

#define NODE_SIMPLE     0x01
#define NODE_PACKED     0x02
#define NODE_UNPACKED   0x02
#define NODE_EDITED     0x08 // only meaningful if not SIMPLE and not PACKED

static Node nodeCreate (I32 k, I32 in, I32 out, I32 inOff,I32 outOff) ;
static I32  nodePack (Node *n, U8 *status) ;   // returns new size in bytes
static I32  nodeUnpack (Node *n, U8 *status) ; // returns new size in bytes

static inline int intPut (U8 *u, I32 val) ;   // at end of file for integer byte packing
static inline int intGet (U8 *u, I32 *pval) ; // at end of file for integer byte unpacking

/*********************** create and destroy ****************************/

SyngBWT *syngBWTcreate (int fixedLen)
{
  if (!fixedLen) die ("syngBWT does not yet support variable length sequences") ;
  SyngBWT *sb = new0 (1, SyngBWT) ;
  sb->fixedLen = fixedLen ;
  I64 max = 1 << 20 ;
  sb->node = arrayCreate (max, Node) ;
  sb->status = arrayCreate (max, U8) ;
  return sb ;
}

void syngBWTdestroy (SyngBWT *sb)
{
  arrayDestroy (sb->node) ; // should reall destroy the contents
  arrayDestroy (sb->status) ;
  if (sb->length) arrayDestroy (sb->length) ;
  if (sb->startHash) hashDestroy (sb->startHash) ;
  if (sb->startHashCount) arrayDestroy (sb->startHashCount) ;
  newFree (sb, 1, SyngBWT) ;
}

/*********************** add - the main function ***********************/

// first a couple of utilities

static inline void scSplit (UnpackedNode *un, SyncCount **inList, SyncCount **outList,
			    SyncCount **inGBWT, SyncCount **outGBWT)
{
  if (!un->sc) return ; // don't do this more than once
  SyncCount *t ;
  // don't need to reallocated *inList, because it takes over un->sc
  t = new(un->outN+1,SyncCount) ; memcpy(t,*outList,un->outN*sizeof(SyncCount)) ; *outList = t ;
  t = new(un->inG+2,SyncCount) ;  memcpy(t,*inGBWT,un->inG*sizeof(SyncCount)) ;   *inGBWT = t ;
  t = new(un->outG+2,SyncCount) ; memcpy(t,*outGBWT,un->outG*sizeof(SyncCount)) ; *outGBWT = t ;
  un->sc = 0 ;
  // do this now since we know size
  newFree (0, un->inN + 2*(un->outN + un->inG + un->outG) + 5, SyncCount) ;
}

static inline void packSplit (UnpackedNode *un, U8 **inOffList, U8 **outOffList)
{
  if (!un->packed) return ; // don't do this more than once
  U8 *t = new(un->outN+1,U8) ; memcpy(t,*outOffList,un->outN*sizeof(U8)) ; *outOffList = t ;
  // don't need to reallocated *inOffList, because it takes over un->packed
  un->packed = 0 ;
  newFree (0, un->inN + 2*un->outN + 1, char) ; // do this now since we know size
}

static U32 startCount (SyngBWT *sb, I32 k)
{
  if (!sb->startHash)
    { sb->startHash = hashCreate (8192) ;
      sb->startHashCount = arrayCreate (8192, U32) ;
    }
  int index ;
  hashAdd (sb->startHash, HASH_INT(k), &index) ;
  return array(sb->startHashCount, index, U32)++ ;
}

// this is the core operation to build the GBWT, by inserting path (in, k, out) at k
// j is the offset in the list of paths from in to k; returns next j

//#define DEBUG_ADD // use for detailed debugging - very verbose!

static I32 syngBWTadd (SyngBWT *sb, I32 k, I32 in, I32 inOff, I32 j, I32 out, I32 outOff) 
{
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
  if (!*s)               // make a new simple node
    { assert (j == 0) ;
      *n = nodeCreate (isPositive ? k : -k, in, inOff, out, outOff) ;
      *s = NODE_SIMPLE ;
#ifdef DEBUG_ADD  
      printf ("new simple node - jNext 0\n") ;
#endif
      return 0 ;
    }
  if (*s & NODE_SIMPLE)    // try to update an old simple node
    { SimpleNode *sn = n->simple ;
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
	  return j ;
	}
      else                  // expand out to Unpacked then fall through to code below
	{ UnpackedNode *un = n->unpacked = new(1,UnpackedNode) ;
	  un->inN = un->outN = 1 ;
	  un->inG = un->outG = 1 ;
	  SyncCount *sc = un->sc = new (4, SyncCount) ;
	  sc->sync = sn->in ; sc->count = sn->inCount ; ++sc ;    // inList
	  sc->sync = sn->out ; sc->count = sn->outCount ; ++sc ;  // outList
	  sc->sync = 0 ; sc->count = sn->outCount ; ++sc ;        // inGBWT
	  sc->sync = 0 ; sc->count = sn->inCount ;                // outGBWT
	  *s = NODE_UNPACKED | NODE_EDITED ; // it will be edited, and this makes ->packed simpler
	  un->packed = new(2,U8) ; un->packed[0] = sn->inOff ; un->packed[1] = sn->outOff ;
#ifdef DEBUG_ADD  
	  printf ("unpacking simple node - ") ;
#endif
	}
    }
  else if (*s & NODE_PACKED)
    { nodeUnpack (n, s) ;
#ifdef DEBUG_ADD  
      printf ("unpacking packed node - ") ;
#endif
    }
  
  // now node must be unpacked
  UnpackedNode *un = n->unpacked ;

  // first check whether it is in edit status, and if not put it there
  if (!(*s & NODE_EDITED)) // must change packed to just contain the offsets
    { U8 *p = un->packed ;
      U8 *u = p ;
      I32 oldSize = *(I32*)u ;
      u += sizeof(I32) ;
      I32 dummy ;
      u += intGet (u, &dummy) ;
      u += intGet (u, &dummy) ;
      U8 *newPacked = new(un->inN + un->outN, U8) ;
      memcpy (newPacked, u, un->inN + un->outN) ;
      un->packed = newPacked ;
      newFree (p, oldSize, U8) ;
      *s |= NODE_EDITED ;
#ifdef DEBUG_ADD  
      printf ("setting edit mode - ") ;
#endif
    }

  // strategy is to conceptually separate out the SyncCount lists as below
  // update them in place if we can
  // if we need to extend one or more of them scSplit() will reassign them and set un->sc to 0
  SyncCount *inList = un->sc, *outList = inList + un->inN ;
  SyncCount *inGBWT = outList + un->outN, *outGBWT = inGBWT + un->inG ;
  // similarly for the offset lists, with packSplit() reassigning and setting un->packed to 0
  U8        *inOffList = un->packed, *outOffList = un->packed + un->inN ;

  // next find in within inList, or add to the end of the list if necessary
  int inK = 0, inPre = 0 ; // inK is index of in within inList, and inPre is sum(count) before inK
  while (inK < un->inN && (inList[inK].sync != in || inOffList[inK] != inOff))
    inPre += inList[inK++].count ;
  if (inK == un->inN)      // add in to inList
    { scSplit (un, &inList, &outList, &inGBWT, &outGBWT) ;
      packSplit (un, &inOffList, &outOffList) ;
      inList[inK].sync = in ; inList[inK].count = 0 ; inOffList[inK] = inOff ;
      ++un->inN ;
#ifdef DEBUG_ADD  
      printf ("adding %d to inList - ", in) ;
#endif
    }

  // now do the same for out
  int outK = 0, outPre = 0 ;
  while (outK < un->outN && (outList[outK].sync != out || outOffList[outK] != outOff))
    outPre += outList[outK++].count ;
  if (outK == un->outN) // add a sync to outList
    { scSplit (un, &inList, &outList, &inGBWT, &outGBWT) ;
      packSplit (un, &inOffList, &outOffList) ;
      outList[outK].sync = out ; outList[outK].count = 0 ; outOffList[outK] = outOff ;
      ++un->outN ;
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
  if (g-g0 >= (isPositive ? un->outG : un->inG)) die ("error in addNode") ;
  if (g->sync == target)                                    // within matching block - easy
    { jNext -= sumCount - pre ;
      g->count++ ;
#ifdef DEBUG_ADD  
      printf ("increment g block %d - ", (int)(g-g0)) ;
#endif
    }
  else if (sumCount + g->count == pre)                      // at the end of this GBWT block
    { if (++g - g0 < (isPositive?un->outG:un->inG) && g->sync == target) // can increment next block
	{ g->count++ ;                                      // jNext is already correct
#ifdef DEBUG_ADD  
	  printf ("increment %s-block %d - ", isPositive?"out":"in", (int)(g-g0)) ;
#endif
	}
      else                                                  // add one row to the GBWT
	{ scSplit (un, &inList, &outList, &inGBWT, &outGBWT) ;
	  if (isPositive)                                   // remember we incremented g!
	    { int z ; for (z = un->outG ; z > g - g0 ; --z) outGBWT[z] = outGBWT[z-1] ;
	      outGBWT[z].sync = target ; outGBWT[z].count = 1 ; // again jNext is correct
#ifdef DEBUG_ADD  
	      printf ("adding out-block %d - ", z) ;
#endif
	      ++un->outG ;
	    }
	  else
	    { int z ; for (z = un->inG ; z > g - g0 ; --z) inGBWT[z] = inGBWT[z-1] ;
	      inGBWT[z].sync = target ; inGBWT[z].count = 1 ;
#ifdef DEBUG_ADD  
	      printf ("adding in-block %d - ", z) ;
#endif
	      ++un->inG ;
	    }
	}
    }
  else // we are in the middle of a GBWT block - we need to add two rows: new one and reversal
    { scSplit (un, &inList, &outList, &inGBWT, &outGBWT) ;
      if (isPositive)
	{ int z ; for (z = un->outG+1 ; z > g - g0 + 2 ; --z) outGBWT[z] = outGBWT[z-2] ;
	  outGBWT[z].sync = outGBWT[z-2].sync ; outGBWT[z].count = sumCount + g->count - pre ; --z ;
	  outGBWT[z].sync = target ; outGBWT[z].count = 1 ; --z ; // again jNext is correct
	  outGBWT[z].count -= (sumCount + g->count - pre) ;
#ifdef DEBUG_ADD  
	  printf ("adding 2 out-blocks %d,%d - ", z, z+1) ;
#endif
	  un->outG += 2 ;
	}
      else
	{ int z ; for (z = un->inG+1 ; z > g - g0 + 2 ; --z) inGBWT[z] = inGBWT[z-2] ;
	  inGBWT[z].sync = inGBWT[z-2].sync ; inGBWT[z].count = sumCount + g->count - pre ; --z ;
	  inGBWT[z].sync = target ; inGBWT[z].count = 1 ; --z ;
	  inGBWT[z].count -= (sumCount + g->count - pre) ;
#ifdef DEBUG_ADD  
	  printf ("adding 2 in-blocks %d,%d - ", z, z+1) ;
#endif
	  un->inG += 2 ;
	}
    }

#ifdef DEBUG_ADD
  int z ;
  if (un->inG) { printf ("inG ") ; for (z = 0 ; z < un->inG ; z++) printf ("(%d,%d)", inGBWT[z].sync, inGBWT[z].count) ; }
  if (un->outG) { printf ("outG ") ; for (z = 0 ; z < un->outG ; z++) printf ("(%d,%d)", outGBWT[z].sync, outGBWT[z].count) ; }
  for (z = 0 ; z < un->inG ; z++) if (inGBWT[z].count < 0) die ("count < 0") ;
  for (z = 0 ; z < un->outG ; z++) if (outGBWT[z].count < 0) die ("count < 0") ;
#endif

  // we are done
  if (!un->sc) // have to reform it
    { un->sc = new (un->inN + un->outN + un->inG + un->outG, SyncCount) ;
      memcpy (un->sc, inList, un->inN*sizeof(SyncCount)) ;
      memcpy (un->sc + un->inN, outList, un->outN*sizeof(SyncCount)) ;
      memcpy (un->sc + un->inN + un->outN, inGBWT, un->inG*sizeof(SyncCount)) ;
      memcpy (un->sc + un->inN + un->outN + un->inG, outGBWT, un->outG*sizeof(SyncCount)) ;
      // call free() here because we already allowed for removal on creation
      free (inList) ; free (outList) ; free (inGBWT) ; free (outGBWT) ;
    }
  if (!un->packed) // have to reform it
    { un->packed = new (un->inN + un->outN, U8) ;
      memcpy (un->packed, inOffList, un->inN) ;
      memcpy (un->packed + un->inN, outOffList, un->outN) ;
      // call free() here because we already allowed for removal on creation
      free (inOffList) ; free (outOffList) ;
    }
  
#ifdef DEBUG_ADD  
  printf (" - jNext %d\n", jNext) ;
#endif
    
  return jNext ;
}

SyngBWTpath *syngBWTpathStart (SyngBWT *sb, I32 startNode)
{
  SyngBWTpath *sbp = new0 (1, SyngBWTpath) ;
  sbp->sb = sb ;
  sbp->lastNode = 0 ;
  sbp->lastOff = 0 ;
  sbp->thisNode = sbp->startNode = startNode ;
  sbp->inCount = sbp->startCount = startCount (sb, startNode) ;
  return sbp ;
}

void syngBWTpathAdd (SyngBWTpath *sbp, I32 nextNode, I32 offset)
{
  sbp->inCount = syngBWTadd (sbp->sb, sbp->thisNode, sbp->lastNode, sbp->lastOff, sbp->inCount,
			     nextNode, offset) ;
  sbp->lastNode = sbp->thisNode ;
  sbp->lastOff = offset ;
  sbp->thisNode = nextNode ;
}

void syngBWTpathFinish (SyngBWTpath *sbp)
{ syngBWTadd (sbp->sb, sbp->thisNode, sbp->lastNode, sbp->lastOff, sbp->inCount, 0, 0) ;
  newFree (sbp, 1, SyngBWTpath) ;
}

/***************** read and write the SyngBWT *******************/

void syngBWTwrite (OneFile *of, SyngBWT *sb)
{
  int i, j, bufSize = 0 ;
  I64 *buf = 0 ; 
  
  for (i = 1 ; i < arrayMax(sb->node) ; ++i)
    { Node n = arr(sb->node, i, Node) ;
      U8   s = arr(sb->status, i, U8) ;
      if (!s) continue ;
      oneInt(of,0) = sb->fixedLen ; oneWriteLine (of, 'V', 0, 0) ;
      if (s & NODE_SIMPLE)
	{ SimpleNode *sn = n.simple ;
	  oneInt(of,1) = sn->in ; oneInt(of,2) = sn->inOff ; oneInt(of,3) = sn->inCount ;
	  oneChar(of,0) = '+' ; oneWriteLine(of, 'E', 0, 0) ; // edge
	  oneInt(of,1) = sn->out ; oneInt(of,2) = sn->outOff ; oneInt(of,3) = sn->outCount ;
	  oneChar(of,0) = '-' ; oneWriteLine(of, 'E', 0, 0) ; // edge
	}
      else
	{ if (s & NODE_PACKED) nodeUnpack (&n, &s) ;
	  UnpackedNode *un = n.unpacked ;
	  SyncCount *inList = un->sc, *outList = inList + un->inN ;
	  SyncCount *inGBWT = outList + un->outN, *outGBWT = inGBWT + un->inG ;
	  U8        *inOffList = un->packed + ((s & NODE_EDITED) ? 0 : un->packOff) ;
	  U8        *outOffList = inOffList + un->inN ;
	  for (j = 0 ; j < un->inN ; ++j) // incoming edges
	    { oneChar(of,0) = '+'; oneInt(of,1) = inList[j].sync; oneInt(of,2) = inOffList[j] ; 
	      oneInt(of,3) = inList[j].count ; oneWriteLine(of, 'E', 0, 0) ; // edge
	    }
	  if (un->outG > 1) // outgoing GBWT, which goes from the incoming edges
	    { if (un->outG > bufSize)
		{ if (bufSize) newFree(buf, bufSize, I64) ;
		  bufSize = 2*un->outG ; buf = new (bufSize, I64) ;
		}
	      for (j = 0 ; j < un->outG ; ++j) buf[j] = outGBWT[j].sync ;
	      oneChar(of,0) = '+' ; oneWriteLine (of, 'B', un->outG, buf) ;
	      for (j = 0 ; j < un->outG ; ++j) buf[j] = outGBWT[j].count ;
	      oneChar(of,0) = '+' ; oneWriteLine (of, 'C', un->outG, buf) ;
	    }
	  for (j = 0 ; j < un->outN ; ++j) // outgoing edges
	    { oneChar(of,0) = '-'; oneInt(of,1) = outList[j].sync; oneInt(of,2) = outOffList[j] ; 
	      oneInt(of,3) = outList[j].count ; oneWriteLine(of, 'E', 0, 0) ; // edge
	    }
	  if (un->inG > 1) // incoming GBWT, which follows from the outgoing edges (in reverse)
	    { if (un->inG > bufSize)
		{ if (bufSize) newFree(buf, bufSize, I64) ;
		  bufSize = 2*un->inG ; buf = new (bufSize, I64) ;
		}
	      for (j = 0 ; j < un->inG ; ++j) buf[j] = inGBWT[j].sync ;
	      oneChar(of,0) = '-' ; oneWriteLine (of, 'B', un->inG, buf) ;
	      for (j = 0 ; j < un->inG ; ++j) buf[j] = inGBWT[j].count ;
	      oneChar(of,0) = '-' ; oneWriteLine (of, 'C', un->inG, buf) ;
	    }
	}
    }
  newFree (buf, bufSize, I64) ;
}

/**************** basic node creation, packing and unpacking ******************/

static Node nodeCreate (I32 k, I32 in, I32 inOff, I32 out, I32 outOff)
{
  Node node ;
  SimpleNode *sn = node.simple = new0 (1, SimpleNode) ;
  sn->in = in ; sn->inOff = inOff ;
  sn->out = out ; sn->outOff = outOff ;
  if (k > 0) sn->inCount = 1 ; else sn->outCount = 1 ;
  return node ;
}

static I32 nodePack (Node *n, U8 *s) // returns new size in bytes
{
  if (!*s) return 0 ; // nothing in this node
  if (*s & NODE_SIMPLE)
    return sizeof(SimpleNode) ;
  else if (*s & NODE_PACKED)
    return *(I32*)(n->packed) ;
  else if (!(*s & NODE_EDITED))
    { UnpackedNode *un = n->unpacked ;
      int scG = un->inN + un->outN + un->inG + un->outG ;
      n->packed = un->packed ;
      *s &= ~NODE_UNPACKED ;
      *s |= NODE_PACKED ;
      newFree (un->sc, scG, SyncCount) ;
      newFree (un, 5*sizeof(I32) + 2*sizeof(void*), char) ;
      return *(I32*)(n->packed) ;
    }
  else // complex and edited
    { UnpackedNode *un = n->unpacked ;
      int scSize = un->inN + un->outN + un->inG + un->outG ;
      int maxSize = sizeof(I32) + 5*(4 + scSize*2) + un->inN + un->outN ;
      U8 *u = new(maxSize,U8) ;
      U8 *u0 = u ;
      u += sizeof(I32) ; // space to record the size of u
      u += intPut (u, un->inN) ;
      u += intPut (u, un->outN) ;
      memcpy (u, un->packed, (size_t) (un->inN + un->outN)) ; // sync offsets
      u += un->inN + un->outN ;
      u += intPut (u, un->inG) ;
      u += intPut (u, un->outG) ;
      int i ;
      SyncCount *sc = un->sc ;
      for (i = 0 ; i < scSize ; ++i, ++sc)
	{ u += intPut (u, sc->sync) ; u += intPut (u, sc->count) ; }
      I32 size = u - u0 ;
      *(I32*)u0 = size ; // store the true size at the start of u
      n->packed = new(size,U8) ;
      memcpy(n->packed,u0,(size_t)size) ;
      newFree (u0, maxSize, U8) ;
      *s &= ~NODE_UNPACKED ;
      *s |= NODE_PACKED ;
      newFree (un->sc, scSize, SyncCount) ;
      newFree (un->packed, un->inN + un->outN, U8) ; // held offsets
      newFree (un, 5*sizeof(I32) + 2*sizeof(void*), char) ;
      return *(I32*)(n->packed) ;
    }
}

static I32 nodeUnpack (Node *n, U8 *s) // returns new size in bytes
{
  if (!*s) return 0 ;   // nothing in this node
  if (*s & NODE_SIMPLE) // don't do anything
    return sizeof(SimpleNode) ;
  else if (*s & NODE_UNPACKED) // also don't do anything - must calculate size
    { UnpackedNode *un = n->unpacked ;
      int scSize = un->inN + un->outN + un->inG + un->outG ;
      I32 size = 5*sizeof(I32) + 2*sizeof(void*) + scSize*sizeof(SyncCount) ;
      if (*s & NODE_EDITED)
	size += un->inN + un->outN ;
      else
	size += *(I32*)un->packed ;
      return size ;
    }
  else if (*s & NODE_PACKED) // packed
    { U8 *u0 = n->packed ;
      U8 *u = u0 ;
      UnpackedNode *un = new (1, UnpackedNode) ;
      un->packed = n->packed ;
      n->unpacked = un ;
      *s &= ~NODE_PACKED ;
      *s |= NODE_UNPACKED ;
      *s &= ~NODE_EDITED ;
      u += sizeof(I32) ; // skip past the initial I32 size
      u += intGet (u, &un->inN) ;
      u += intGet (u, &un->outN) ;
      un->packOff = u - u0 ;
      u += un->inN + un->outN ;
      u += intGet (u, &un->inG) ;
      u += intGet (u, &un->outG) ;
      int scSize = un->inN + un->outN + un->inG + un->outG ;
      un->sc = new (scSize, SyncCount) ;
      int i ;
      SyncCount *sc = un->sc ;
      for (i = 0 ; i < scSize ; ++i, ++sc)
	{ u += intPut (u, sc->sync) ; u += intPut (u, sc->count) ; }
      
      return sizeof(UnpackedNode) +
	5*sizeof(I32) + 2*sizeof(void*) + scSize*sizeof(SyncCount) + *(I32*)un->packed ;
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
