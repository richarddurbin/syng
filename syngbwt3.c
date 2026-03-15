/*  File: syngbwt.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 15 16:23 2026 (rd109)
 * * Nov 23 01:15 2025 (rd109): converted to skipList (balanced tree)
 * Created: Mon Sep  9 11:34:51 2024 (rd109)
 *-------------------------------------------------------------------
 */

static int DEBUG = 0 ; // set to the sync you want to debug, -1 for general, or 0 for none
static int PATH_DEBUG = 0 ;
int pathCount = 0 ; // global for debugging

#include "syng.h"
#include "rskip.h"

typedef union {
  struct {
    I32 sync ;
    U16 offset ; // requires maximum offset (1>>16)-1
    U16 count ;  // if count >= 1 >> 16 then use rs
  } ;
  Rskip rs ; // internally a 64-bit pointer of some sort
} NodeSide ; // 64 bits

typedef struct {
  NodeSide in ;
  NodeSide out ;
} Node ;     // 128 bits

#define NODE_EXISTS      0x01
#define NODE_SIMPLE_IN   0x02
#define NODE_SIMPLE_OUT  0x04
#define NODE_SIMPLE      0x07

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
    { if (!(arr(sb->status,i,U8) & NODE_SIMPLE_IN))  rsDestroy (arr(sb->node,i,Node).in.rs) ;
      if (!(arr(sb->status,i,U8) & NODE_SIMPLE_OUT)) rsDestroy (arr(sb->node,i,Node).out.rs) ;
    }
  arrayDestroy (sb->node) ;
  arrayDestroy (sb->status) ;
  if (sb->length) arrayDestroy (sb->length) ;
  if (sb->startHash) hashDestroy (sb->startHash) ;
  if (sb->startHashCount) arrayDestroy (sb->startHashCount) ;
  newFree (sb, 1, SyngBWT) ;
}

/*********************** simple node creation *************************/

static Node nodeCreate (I32 in, U32 inOff, I32 out, U32 outOff)
{
  Node node ;
  node.in.sync = in ;   node.in.offset = inOff ; node.in.count = 0 ;
  node.out.sync = out ; node.out.offset = outOff ; node.out.count = 0 ;
  return node ;
}

static inline void nodePrint (Node *n, U8 s)
{
  if (s & NODE_SIMPLE_IN)
    printf   ("  in sync %d offset %u count %u\n", n->in.sync, n->in.offset, n->in.count) ;
  else
    { printf ("  in rs ") ; rsPrint (n->in.rs) ; }
  if (s & NODE_SIMPLE_OUT)
    printf   ("  out sync %d offset %u count %u\n", n->out.sync, n->out.offset, n->out.count) ;
  else
    { printf ("  out rs ") ; rsPrint (n->out.rs) ; }
}

static inline bool nodeCheck (Node *n, U8 s)
{
  bool isSimpleIn = (s & NODE_SIMPLE_IN), isSimpleOut = (s & NODE_SIMPLE_OUT) ;
  bool isBad = false ;

  if (!s) return true ;

  if (isSimpleIn && !isSimpleOut && n->in.count != rsLength(n->out.rs))
    { isBad = true ;
      warn ("nodeCheck simple in.count %u != length out.rs %d", n->in.count, rsLength(n->out.rs)) ;
    }
  else if (isSimpleOut && !isSimpleIn && n->out.count != rsLength(n->in.rs))
    { isBad = true ;
      warn ("nodeCheck simple out.count %u != length in.rs %d", n->out.count, rsLength(n->in.rs)) ;
    }
  else if (!isSimpleIn && !isSimpleOut)
    { if (rsDirSum (n->in.rs) != rsLength(n->out.rs))
	{ isBad = true ;
	  warn ("nodeCheck dirSum in %u != length out %d", rsDirSum(n->in.rs), rsLength(n->out.rs)) ;
	}
      if (rsDirSum (n->out.rs) != rsLength(n->in.rs))
	{ isBad = true ;
	  warn ("nodeCheck dirSum out %u != length in %d", rsDirSum(n->out.rs), rsLength(n->in.rs)) ;
	}
    }
  if (!isSimpleIn && !rsCheck (n->in.rs)) isBad = true ;
  if (!isSimpleOut && !rsCheck (n->out.rs)) isBad = true ;
  return isBad ;
}

/******************** manage the starts *********************/

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

/*********************** main functions to add and follow *************************/

static inline void convertToRskip (NodeSide *ns, I32 sync, U32 off, U32 sum)
{
  I32 symbol[2] = { ns->sync, sync } ;
  U32 offset[2] = { ns->offset, off } ;
  U32 oldCount = ns->count ;
  if (sum)
    { I64 iSym = 0, runLen = sum ;
      ns->rs = rsBuildDynamicSyng (2, symbol, offset, 1, &iSym, &runLen) ;
    }
  else
    ns->rs = rsBuildDynamicSyng (2, symbol, offset, 0, 0, 0) ;
  if (oldCount) rsDirSetCount (ns->rs, 0, oldCount) ;
}

static inline void convertBig (NodeSide *ns, U32 sum)
{
  I32 symbol = ns->sync ;
  U32 offset = ns->offset ;
  U32 oldCount = ns->count ;
  if (sum)
    { I64 iSym = 0, runLen = sum ;
      ns->rs = rsBuildDynamicSyng (1, &symbol, &offset, 1, &iSym, &runLen) ;
    }
  else
    ns->rs = rsBuildDynamicSyng (1, &symbol, &offset, 0, 0, 0) ;
  if (oldCount) rsDirSetCount (ns->rs, 0, oldCount) ;
}

static U32 syngBWTadd (SyngBWT *sb, I32 k, I32 in, U32 inOff, U32 j, I32 out, U32 outOff)
{ // returns the next value of j
  bool  isPositive = (k >= 0) ;
  I32   kPos = isPositive ? k : -k ;
  Node *n = arrayp(sb->node,   kPos, Node) ;
  U8   *s = arrayp(sb->status, kPos, U8) ;
  static int MAX = 65000 ;

  static U32 callCount = 0 ; ++callCount ;
  if (DEBUG == kPos)
    printf ("syngBWTadd %u k %d s %u in %d inOff %u j %u out %d outOff %u\n",
	    callCount, k, *s, in, inOff, j, out, outOff) ;
  
  U32 jOrig = j ;
  if ((I32)jOrig < 0) die ("syngBWTadd negative input j %d callCount %d k %d", (I32)j, callCount, k) ;

  if (callCount == 0) // set to bad callCount for debugging
    printf ("!!wierdness here\n") ;
  
  if (!*s)               // empty - make a new simple node
    { assert (j == 0) ;
      if (isPositive)
	{ *n = nodeCreate (in, inOff, out, outOff) ;
	  n->in.count = 1 ;
	}
      else
	{ *n = nodeCreate (-out, outOff, -in, inOff) ;
	  n->out.count = 1 ;
	}
      *s = NODE_SIMPLE ;
     return 0 ;
    }

  int retVal = 0 ;

  if (isPositive)
    { if ((*s & NODE_SIMPLE_IN) && n->in.count >= MAX)
	{ U32 outSum = (*s & NODE_SIMPLE_OUT) ? n->out.count : rsDirSum (n->out.rs) ;
	  convertBig (&n->in, outSum) ; *s &= ~NODE_SIMPLE_IN ;
	}
      if ((*s & NODE_SIMPLE_IN) && n->in.sync == in && n->in.offset == inOff)
        ++n->in.count ;
      else if (*s & NODE_SIMPLE_IN) // convert to Rskip
	{ assert (j == 0) ; // a new in->inOff for this node, so must be 0
	  j = n->in.count ;
	  U32 outSum = (*s & NODE_SIMPLE_OUT) ? n->out.count : rsDirSum (n->out.rs) ;
	  convertToRskip (&n->in, in, inOff, outSum) ;
	  *s &= ~NODE_SIMPLE_IN ;
	  rsDirSetCount (n->in.rs, 1, 1) ;
	}
      else // already an Rskip
	j += rsDirAddSyng (&n->in.rs, in, inOff) ; // returns rank

      if ((*s & NODE_SIMPLE_OUT) && n->out.sync == out && n->out.offset == outOff)
        retVal = j ;
      else
	{ if (*s & NODE_SIMPLE_OUT) // convert to Rskip
	    { U32 inSum = (*s & NODE_SIMPLE_IN) ? n->in.count : rsDirSum (n->in.rs) ;
	      convertToRskip (&n->out, out, outOff, inSum-1) ; // -1 because we already incremented
	      *s &= ~NODE_SIMPLE_OUT ;
	    }
	  retVal = rsAddSyng (&n->out.rs, j, out, outOff) ;
	}
    }
  else // have to come from .out == -in and leave on .in == -out
    { if ((*s & NODE_SIMPLE_OUT) && n->out.count >= MAX)
	{ U32 inSum = (*s & NODE_SIMPLE_IN) ? n->in.count : rsDirSum (n->in.rs) ;
	  convertBig (&n->out, inSum) ; *s &= ~NODE_SIMPLE_OUT ;
	}
      if ((*s & NODE_SIMPLE_OUT) && n->out.sync == -in && n->out.offset == inOff)
	++n->out.count ;
      else if (*s & NODE_SIMPLE_OUT) // convert to Rskip
	{ assert (j == 0) ; // a new in->inOff for this node, so must be 0
	  j = n->out.count ;
	  U32 inSum = (*s & NODE_SIMPLE_IN) ? n->in.count : rsDirSum (n->in.rs) ;
	  convertToRskip (&n->out, -in, inOff, inSum) ;
	  *s &= ~NODE_SIMPLE_OUT ;
	  rsDirSetCount (n->out.rs, 1, 1) ;
	}
      else // already an Rskip
	j += rsDirAddSyng (&n->out.rs, -in, inOff) ; // returns rank
      
      if ((*s & NODE_SIMPLE_IN) && n->in.sync == -out && n->in.offset == outOff)
	retVal = j ;
      else
	{ if (*s & NODE_SIMPLE_IN) // convert to Rskip
	    { U32 outSum = (*s & NODE_SIMPLE_OUT) ? n->out.count : rsDirSum (n->out.rs) ;
	      convertToRskip (&n->in, -out, outOff, outSum-1) ; // -1 because we already incremented
	      *s &= ~NODE_SIMPLE_IN ;
	    }
	  retVal = rsAddSyng (&n->in.rs, j, -out, outOff) ;
	}
    }

  if (DEBUG && nodeCheck (n, *s))
    { nodePrint (n, *s) ;
      printf ("syngBWTadd %u k %d s %u in %d inOff %u j %u out %d outOff %u\n",
	      callCount, k, *s, in, inOff, j, out, outOff) ;
      die ("syngBWTadd nodeCheck fail callCount %d k retVal %d", callCount, retVal) ;
    }
  
  if (DEBUG == kPos) nodePrint (n, *s) ;
  if (retVal < 0) die ("syngBWTadd -ve output retVal %d callCount %d k %d", retVal, callCount, k) ;

  return (U32)retVal ;
}

static U32 syngBWTnext (SyngBWT *sb, I32 k, I32 in, U32 inOff, U32 j, I32 *out, U32 *outOff)
{
  bool isPositive = (k >= 0) ;
  I32  kPos = isPositive ? k : -k ;
  if (kPos >= arrayMax(sb->node))
    die ("syngBWTnext: k %d >= arrayMax(sb->node) %lld", k, arrayMax(sb->node)) ;
  Node n = arr(sb->node,   kPos, Node) ;
  U8   s = arr(sb->status, kPos, U8) ;
  if (!s) die ("syngBWTnext %d hit an empty node", k) ;

  if (pathCount == PATH_DEBUG)
    { printf ("++path %d k %d in %d j %u\n", PATH_DEBUG, k, in, j) ; nodePrint (&n, s) ; }

  if (isPositive)
    { if (s & NODE_SIMPLE_IN)
	{ if (n.in.sync != in || n.in.offset != inOff || j >= n.in.count) 
	    die ("syngBWTnext SIMPLE_IN mismatch k %d", k) ;
	}
      else
	j += rsDirRankSyng (n.in.rs, in, inOff) ; // returns rank - dies if not found

      if (s & NODE_SIMPLE_OUT)
	{ *out = n.out.sync ; *outOff = n.out.offset ;
	  return j ;
	}
      else // already an Rskip
	return rsFindSyng (n.out.rs, j, out, outOff) ;
    }
  else // have to come from .out == -in and leave on .in == -out
    { if (s & NODE_SIMPLE_OUT)
	{ if (n.out.sync != -in || n.out.offset != inOff || j >= n.out.count) 
	    die ("syngBWTnext SIMPLE_OUT mismatch k %d", k) ;
	}
      else // already an Rskip
	j += rsDirRankSyng (n.out.rs, -in, inOff) ; // returns rank - dies if not found

      if (s & NODE_SIMPLE_IN)
	{ *out = -n.in.sync ; *outOff = n.in.offset ;
	  return j ;
	}
      else // already an Rskip
	{ U32 jOut = rsFindSyng (n.in.rs, j, out, outOff) ;
	  *out = -*out ;
	  return jOut ;
	}
    }
}

static inline bool syngBWTmatch (SyngBWT *sb, I32 k, I32 out, U32 off, U32 *low, U32 *high)
{
  bool isPositive = (k >= 0) ;
  I32  kPos = isPositive ? k : -k ;
  if (kPos >= arrayMax(sb->node))
    die ("syngBWTnext: k %d >= arrayMax(sb->node) %lld", k, arrayMax(sb->node)) ;
  Node n = arr(sb->node,   kPos, Node) ;
  U8   s = arr(sb->status, kPos, U8) ;
  if (!s) { printf ("MATCH at %d node empty!\n", k) ; return false ; }

  if (pathCount == PATH_DEBUG)
    { printf ("++match %d k %d out %d low %u high %u\n", PATH_DEBUG, k, out, *low, *high) ;
      nodePrint (&n, s) ;
    }
  
  if (isPositive)
    { if (s & NODE_SIMPLE_OUT)
	if (out == n.out.sync && off == n.out.offset)
	  { if (*high > n.out.count) die ("syngBWTmatch error") ; }
	else return false ;
      else
	{ U32 newLow = rsRankSyng (n.out.rs, *low, out, off) ; // returns rank - 0 if not found
	  U32 newHigh = rsRankSyng (n.out.rs, *high, out, off) ;
	  if (newHigh == newLow) return false ;
	  *high = newHigh ; *low = newLow ;
	}
      if (*high > *low)
	{ if (out > 0 && !(arr(sb->status, out, U8) & NODE_SIMPLE_IN))
	    { int dr = rsDirRankSyng (arrp(sb->node, out, Node)->in.rs, k, off) ;
	      *high += dr ; *low += dr ;
	    }
	  else if (out < 0 && !(arr(sb->status, -out, U8) & NODE_SIMPLE_OUT))
	    { int dr = rsDirRankSyng (arrp(sb->node, -out, Node)->out.rs, -k, off) ;
	      *high += dr ; *low += dr ;
	    }
	}
      }
  else
    { if (s & NODE_SIMPLE_IN)
	if (-out == n.in.sync && off == n.in.offset)
	  { if (*high > n.in.count) die ("syngBWTmatch error") ; }
	else return false ;
      else
	{ U32 newLow = rsRankSyng (n.in.rs, *low, -out, off) ; // returns rank - 0 if not found
	  U32 newHigh = rsRankSyng (n.in.rs, *high, -out, off) ;
	  if (newHigh == newLow) return false ;
	  *high = newHigh ; *low = newLow ;
	}
      if (*high > *low)
	{ if (out > 0 && !(arr(sb->status, out, U8) & NODE_SIMPLE_IN))
	    { int dr = rsDirRankSyng (arrp(sb->node, out, Node)->in.rs, k, off) ;
	      *high += dr ; *low += dr ;
	    }
	  else if (out < 0 && !(arr(sb->status, -out, U8) & NODE_SIMPLE_OUT))
	    { int dr = rsDirRankSyng (arrp(sb->node, -out, Node)->out.rs, -k, off) ;
	      *high += dr ; *low += dr ;
	    }
	}
    }
  return (*high > *low) ;
}

/*******************************************************/
/************* external interfaces *********************/

static SyngBWTpath *pathCreate (SyngBWT *sb, I32 startNode)
{
  SyngBWTpath *sbp = new0 (1, SyngBWTpath) ;
  sbp->sb = sb ;
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

void syngBWTpathAdd (SyngBWTpath *sbp, I32 nextNode, U32 offset)
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

SyngBWTpath *syngBWTpathStartOld (SyngBWT *sb, I32 startNode, U32 count)
{
  SyngBWTpath *sbp = pathCreate (sb, startNode) ;
  sbp->jLast = count ;
  if (count >= startCount (sb, startNode, false))
    die ("syngBWTpathStartOld startNode %d count %u >= startCount %u",
	 startNode, count, startCount(sb,startNode,false)) ;
  return sbp ;
}

bool syngBWTpathNext (SyngBWTpath *sbp, I32 *nextNode, U32 *offset)
{
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

/************* interface to match a new sequence to the BWT, e.g. to find MEMs *******************/

SyngBWTpath *syngBWTmatchStart (SyngBWT *sb, I32 startNode, U32 *high)
{
  SyngBWTpath *sbp = pathCreate (sb, startNode) ;
  if (high)
    { bool  isPositive = (startNode >= 0) ;
      I32   kPos = isPositive ? startNode : -startNode ;
      Node *n = arrp(sb->node,  kPos, Node) ;
      U8    s = arr(sb->status, kPos, U8) ;
      if (isPositive)
	if (s & NODE_SIMPLE_IN) *high = n->in.count ;
	else *high = rsDirSum (n->in.rs) ; // could be rsLength (n.out.rs) ;
      else
	if (s & NODE_SIMPLE_OUT) *high = n->out.count ;
	else *high = rsDirSum (n->out.rs) ;
    }
  return sbp ;
}

bool syngBWTmatchNext (SyngBWTpath *sbp, I32 nextNode, U32 nextOff, U32 *low, U32 *high)
{
  if (syngBWTmatch (sbp->sb, sbp->thisNode, nextNode, nextOff, low, high))
    { sbp->thisNode = nextNode ;
      return true ;
    }
  else
    return false ;
}

/****************************************************************/
/********************* write the SyngBWT ************************/

void rsWrite (Rskip rs, OneFile *of, char eChar, char bChar, char cChar)
{
  I64 nRun ;
  int i, nSym = rsNsym(rs) ;
  for (i = 0 ; i < nSym ; ++i)
    { rsDirSyng (rs, i, &oneInt(of,0), &oneInt(of,1), &oneInt (of, 2)) ;
      oneWriteLine (of, eChar, 0, 0) ;
    }
  static I64 *iSym, *runLen ;
  static int size = 0 ;
  if (rsSize(rs, 0, 0) > size)
    { if (size) { newFree (iSym, size, I64) ; newFree (runLen, size, I64) ; }
      size = rsSize(rs, 0, 0) ; iSym = new (size, I64) ; runLen = new (size, I64) ;
    }
  rsLinearise (rs, &nRun, iSym, runLen) ;
  oneWriteLine (of, bChar, nRun, iSym) ;
  oneWriteLine (of, cChar, nRun, runLen) ;
}


void syngBWTwrite (OneFile *of, SyngBWT *sb)
{
  for (int i = 1 ; i < arrayMax(sb->node) ; ++i)
    { Node n = arr(sb->node, i, Node) ;
      U8   s = arr(sb->status, i, U8) ;
      oneInt(of,0) = sb->fixedLen ; oneWriteLine (of, 'V', 0, 0) ;
      if (!s) continue ; // must come after we have written the node
      if (s & NODE_SIMPLE_IN)
	{ oneInt(of,0) = n.in.sync ; oneInt(of,1) = n.in.offset ; oneInt(of,2) = n.in.count ;
	  oneWriteLine (of, 'E', 0, 0) ; // + edge
	}
      else
	rsWrite (n.in.rs, of, 'E', 'B', 'C') ;
      if (s & NODE_SIMPLE_OUT)
	{ oneInt(of,0) = n.out.sync ; oneInt(of,1) = n.out.offset ; oneInt(of,2) = n.out.count ;
	  oneWriteLine (of, 'e', 0, 0) ; // + edge
	}
      else
	rsWrite (n.out.rs, of, 'e', 'b', 'c') ;
    }
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
  I32 *inSync   = new (EMax, I32),  *outSync   = new (eMax, I32) ;
  U32 *inOffset = new (EMax, U32),  *outOffset = new (eMax, U32) ;
  U32 *inSum    = new (EMax, U32),  *outSum    = new (eMax, U32) ;
  I64 *inSym, *inRunLen, *outSym, *outRunLen ;

  if (!oneGoto (of, 'V', rt->v1)) die ("failed to locate to V line %u", rt->v1) ;
  oneReadLine (of) ; // read the first V line
  rt->eTotal = 0 ;
  for (i = rt->v1 ; i <= rt->vn ; ++i)
    { if (oneInt(of,0) != sb->fixedLen)
	die ("syngBWTread node length %lld != expected %d", oneInt(of,0), sb->fixedLen) ;
      Node *n = arrp(sb->node, i, Node) ;
      U8   *s = arrp(sb->status, i, U8) ;
      int   inN = 0, outN = 0, inNrun = 0, outNrun = 0 ;
      while (oneReadLine (of) && of->lineType != 'V')
	switch (of->lineType)
	  { 
	  case 'E':
	    inSync[inN] = oneInt(of,0) ; inOffset[inN] = oneInt(of,1) ; inSum[inN] = oneInt(of,2) ;
	    if (inSync[inN] == 0) startCountAdd (sb, i, inSum[inN]) ; // start node
	    ++inN ; ++rt->eTotal ; break ;
	  case 'e':
	    outSync[outN] = oneInt(of,0) ; outOffset[outN] = oneInt(of,1) ; outSum[outN] = oneInt(of,2) ;
	    if (outSync[outN] == 0) startCountAdd (sb, -i, outSum[outN]) ; // start node
	    ++outN ; ++rt->eTotal ; break ;
	  case 'B': inNrun = oneLen(of) ; inSym = oneIntList(of) ; break ;
	  case 'b': outNrun = oneLen(of) ; outSym = oneIntList(of) ; break ;
	  case 'C': inRunLen = oneIntList(of) ; break ;
	  case 'c': outRunLen = oneIntList(of) ; break ;
	  default:
	    die ("unrecognized linetype %c in vertex", of->lineType) ;
	  }
      // now we have read everything we need to construct the node
      *s = 0 ; // default for empty node
      if (inN > 0 || outN > 0)
	{ *s = NODE_EXISTS ;
	  if (inN == 1)
	    { *s |= NODE_SIMPLE_IN ;
	      n->in.sync = inSync[0] ; n->in.offset = inOffset[0] ; n->in.count = inSum[0] ;
	    }
	  else
	    n->in.rs = rsBuildFixedSyng (inN, inSync, inOffset, inNrun, inSym, inRunLen) ;
	  if (outN == 1)
	    { *s |= NODE_SIMPLE_OUT ;
	      n->out.sync = outSync[0] ; n->out.offset = outOffset[0] ; n->out.count = outSum[0] ;
	    }
	  else
	    n->out.rs = rsBuildFixedSyng (outN, outSync, outOffset, outNrun, outSym, outRunLen) ;
	}
      if (DEBUG) nodeCheck (n, *s) ;
    }
  
  newFree (inSync, EMax, I32) ;   newFree (outSync, eMax, I32) ;
  newFree (inOffset, EMax, U32) ; newFree (outOffset, eMax, U32) ;
  newFree (inSum, EMax, U32) ;    newFree (outSum, eMax, I32) ;
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
  printf ("read GBWT with %lld vertices and %lld edges from %s\n", nv, eTotal, oneFileName(of)) ;

  arrayMax (sb->node) = nv+1 ; // need to set these because we filled with arrp() for thread safety
  arrayMax (sb->status) = nv+1 ;

  newFree (threads, nThread, pthread_t) ;
  newFree (rt, nThread, ReadThread) ;
  timeUpdate (stdout) ;
  return sb ;
}

/*****************************************************************************/

void syngBWTstat (SyngBWT *sb)
{
  I64 N = arrayMax(sb->node) ;
  I64 nSimple = 0, maxTotal = 0, nBlank = 0 ;
  Array eHist = arrayCreate (1024, I64) ;
  Array cHist = arrayCreate (1024, I64) ;
  //  Array gHist = arrayCreate (1024, I64) ;
  int i, j ;
  I64 eSum = 0, eTotal = 0 ;
  for (i = 1 ; i < N ; ++i)
    { Node n = arr(sb->node, i, Node) ;
      U8   s = arr(sb->status, i, U8) ;

      if (!s) { ++nBlank ; continue ; }

      if (s & NODE_SIMPLE_IN)
	++nSimple ;
      else
	{ int nSym = rsNsym(n.in.rs), total = 0 ;
	  I64 count ;
	  ++array(eHist,nSym,I64) ;
	  for (j = 0 ; j < nSym ; ++j)
	    { rsDirSyng (n.in.rs, j, 0, 0, &count) ;
	      ++array(cHist,count,I64) ;
	      total += count ;
	      eSum += count ;
	      eTotal += count * (j+1) ;
	    }
	  if (total > maxTotal) maxTotal = total ;
	}

      if (s & NODE_SIMPLE_OUT)
	++nSimple ;
      else
	{ int nSym = rsNsym(n.out.rs), total = 0 ;
	  I64 count ;
	  ++array(eHist,nSym,I64) ;
	  for (j = 0 ; j < nSym ; ++j)
	    { rsDirSyng (n.out.rs, j, 0, 0, &count) ;
	      ++array(cHist,count,I64) ;
	      total += count ;
	      eSum += count ;
	      eTotal += count * (j+1) ;
	    }
	  if (total > maxTotal) maxTotal = total ;
	}
    }
  printf ("  %'lld nodes of which %'lld are blank, %'lld sides are simple, %'lld sides are complex %.1f%%\n",
	  N, nBlank, nSimple, 2*N - 2*nBlank - nSimple, (2*N-2*nBlank-nSimple)/(0.02*N)) ;
  printf ("    max in or out edge list %'lld\n", arrayMax(eHist)) ;
  printf ("    max count of any one edge %'lld\n", arrayMax(cHist)) ;
  printf ("    max total paths through a node %'lld\n", maxTotal) ;
  printf ("    expected list search time %.2f\n", eTotal / (1.0*eSum)) ;
  I64 bit = 0 ;
  printf ("edge list dbn per side\n") ;
  printf     ("    %2d %'12lld\n", 1, nSimple) ; bit += nSimple ;
  for (i = 2 ; i < 20 ; ++i)
    { printf ("    %2d %'12lld\n", i, arr(eHist, i, I64)) ; bit += arr(eHist, i, I64) ; }
  I64 left = 0 ; while (i < arrayMax(eHist)) left += arr(eHist, i++, I64) ; 
  printf     ("   >20 %'12lld\n", left) ;
  arrayDestroy (eHist) ; arrayDestroy (cHist) ;
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
