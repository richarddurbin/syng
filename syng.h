/*  File: syng.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  6 11:47 2025 (rd109)
 * Created: Mon May 29 08:19:18 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "kmerhash.h"
#include "ONElib.h"

#define SYNG_VERSION  "2.0"

typedef struct {
  int   fixedLen ;
  Array node ;   // of Node - only defined in syngbwt.c
  Array status ; // of U8 bitFlags
  Array length ; // of I32, if fixedLen == 0
  // the next annoying set of properties are to manage start counts
  Hash  startHash ;      // for nodes with starts - good if relatively few starts compared to nodes
  Array startHashCount ; // of I32, indexed by hash value
} SyngBWT ;

typedef struct {
  SyngBWT *sb ;
  U32      lastNode, lastOff, thisNode, inCount ;
} SyngBWTpath ;

// in syngbwt.c

SyngBWT       *syngBWTcreate (int fixedLen, I64 max) ;
void           syngBWTdestroy (SyngBWT *sb) ;
void           syngBWTwrite (OneFile *of, SyngBWT *sb) ;
SyngBWT       *syngBWTread  (OneFile *of) ;
SyngBWTpath   *syngBWTpathStartNew (SyngBWT *sb, I32 startNode) ;
void           syngBWTpathAdd (SyngBWTpath *sbp, I32 nextNode, I32 offset) ; // add to a new path
void           syngBWTpathFinish (SyngBWTpath *sbp) ; // use when creating a new path
SyngBWTpath   *syngBWTpathStartOld (SyngBWT *sb, I32 startNode, I32 count) ; // follow an existing path
bool           syngBWTpathNext (SyngBWTpath *sbp, I32 *nextNode, I32 *nextPos) ;
void           syngBWTpathDestroy (SyngBWTpath *sbp) ;

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
  "O S 1 3 DNA               sequence of the node\n" // for the future, for general GFA
  ".\n"
  "O V 1 3 INT               graph node (vertex): length\n"
  "D K 1 3 INT               coverage count of node\n"
  "D E 3 3 INT 3 INT 3 INT   edge +: adjacent node (- if reversed), offset, count\n"
  "D e 3 3 INT 3 INT 3 INT   edge -: adjacent node (- if reversed), offset, count\n"
  "D B 1 8 INT_LIST          GBWT +: list of node indices from opposite-signed E lines\n"
  "D b 1 8 INT_LIST          GBWT -: list of node indices from opposite-signed E lines\n"
  "D C 1 8 INT_LIST          GBWT +: list of run-length counts\n"
  "D c 1 8 INT_LIST          GBWT -: list of run-length counts\n"
  ".\n"
  "O P 3 3 INT 3 INT 3 INT   path: length in bp, source file number, sequence number in file\n"
  "D Z 4 3 INT 3 INT 3 INT 3 INT   GBWT path: starting node, pos, count, then length in nodes\n"
  "D z 1 8 INT_LIST          alternative explicit list of node ids (-ve if reversed)\n"
  "D o 1 8 INT_LIST          if z, then offsets of the nodes from start of sequence, 1:1 with z\n"
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
  "P 3 map                   MAP\n"
  "P 3 ref                   REFERENCE INFORMATION\n"
  "O N 1 6 STRING            name of sequence, e.g. chr1\n"
  "O I 1 8 INT_LIST          list of indexes into name table for each sync: 0 if missing, -ve if RC\n"
  "D P 1 8 INT_LIST          positions in sequence\n"
  ;

/****************** end of file ********************/
