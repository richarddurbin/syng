/*  File: syng.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 11 01:17 2024 (rd109)
 * Created: Mon May 29 08:19:18 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "kmerhash.h"
#include "ONElib.h"

typedef struct {
  int   fixedLen ;
  Array node ;   // of Node - only defined in syngbwt.c
  Array status ; // of U8 bitFlags
  Array length ; // of I32, if fixedLen == 0
  // the next annoying set of properties are to manage start counts
  Hash  startHash ; // for nodes with >= 255 starts in at least one direction
  Array startHashCount ; // of I32
  Array startPlusCount, startMinusCount ; // of U8, property of a ReadSet - don't read/write
} SyngBWT ;

typedef struct {
  SyngBWT *sb ;
  I32      startNode, startCount ; // this is all that is needed to reconstruct the node path
  U32      lastNode, lastOff, thisNode, inCount ;
} SyngBWTpath ;

// in syngbwt.c

SyngBWT       *syngBWTcreate (int fixedLen) ;
void           syngBWTdestroy (SyngBWT *sb) ;
void           syngBWTwrite (OneFile *of, SyngBWT *sb) ;
bool           syngBWTread  (OneFile *of, SyngBWT *sb ) ;
SyngBWTpath   *syngBWTpathStart (SyngBWT *sb, I32 startNode) ;
void           syngBWTpathFinish (SyngBWTpath *sbp) ;
void           syngBWTpathAdd (SyngBWTpath *sbp, I32 nextNode, I32 offset) ;

static char *syngSchemaText =
  "1 3 def 1 0               schema for syng\n"
  ".\n"
  "P 3 seq                   SEQUENCE\n"
  "S 4 path                  contains P (path) objects = syncmer sequences\n"
  "S 3 gfa                   sequence graph - contains V (vertex) objects, probably with E lines\n"
  "S 4 gbwt                  gbwt: a gfa with B, C, Z lines\n"
  ".\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  ".\n"
  "O S 1 3 DNA               sequence of the node\n"
  ".\n"
  "O V 1 3 INT               graph node (vertex): length\n"
  "D K 1 3 INT               coverage count of node\n"
  "D E 4 4 CHAR 3 INT 3 INT 3 INT   edge: +/-, adjacent node (- if reversed), offset, count\n"
  "D B 2 4 CHAR 8 INT_LIST          GBWT: +/-, list of node indices from opposite-signed E lines\n"
  "D C 2 4 CHAR 8 INT_LIST          GBWT: +/-, list of run-length counts\n"
  ".\n"
  "O P 2 3 INT 3 INT         path: reference file number and sequence number in file\n"
  "D Z 2 3 INT 3 INT         GBWT path: starting node, starting count, length in nodes\n"
  "D z 1 8 INT_LIST          alternative explicit list of node ids (-ve if reversed)\n"
  "D d 1 6 STRING            directions (orientations) of nodes in path (optional): 1:1 with P\n"
  "D o 1 8 INT_LIST          offsets of the nodes from start of sequence, 1:1 with P\n"
  "D X 1 3 DNA               prefix before first node - required to fully reconstruct\n"
  "D Y 1 3 DNA               suffix after last node - required to fully reconstruct\n"
  ".\n"
  "P 5 khash                 KMER HASH\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  "D t 3 3 INT 3 INT 3 INT   max, len, dim for KmerHash table\n"
  "O S 1 3 DNA               packed sequences aligned to 64-bit boundaries\n" 
  "D L 1 8 INT_LIST          locations in the table\n"
  "D C 1 8 INT_LIST          kmer counts\n"
  ".\n"
  ;

/****************** end of file ********************/
