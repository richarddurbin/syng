/*  File: seqio.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 19 02:07 2024 (rd109)
 * * Dec 20 23:23 2022 (rd109): fully revised SeqPack and QualPack
 * Created: Sat Nov 10 08:51:49 2018 (rd109)
 *-------------------------------------------------------------------
 */

#ifndef SEQIO_DEFINED
#define SEQIO_DEFINED

#include "utils.h"
#include <zlib.h>

/* first the packages to pack and unpack sequences and qualities, then the IO package */

/* SeqPack interconverts to 2-bit packed binary (4 bases per byte) and works on these */
/* maps acgt, ACGT and 0123 encodings.  N and n are mapped to 0 = a */

typedef struct {
  char unconv[5] ;		/* 5 so we can strcpy() into it for coding simplicity */
  char unconvC[5] ;
  U32  byteExpand[256] ;
  U32  byteExpandC[256] ;
} SeqPack ;

SeqPack *seqPackCreate (char unpackA) ; 
	/* unpackA = 0 maps to 0123, a to acgt, A to ACGT and 1 to 4-bit 1248 */
#define  seqPackDestroy(sp) newFree(sp,1,SeqPack)
U8*      seqPack (SeqPack *sp, char *s, U8 *u, U64 len) ;
U8*      seqPackRevComp (SeqPack *sp, char *s, U8 *u, U64 len) ; /* packs the RC of the sequence s */
	/* compress s into u if non-zero else allocates memory, returns packed array */
char*    seqUnpack (SeqPack *sp, U8 *u, char *s, U64 i, U64 len) ;
char*    seqUnpackRevComp (SeqPack *sp, U8 *u, char *s, U64 i, U64 len) ;
	/* uncompress u into s, memory/return like seqPack(), len and offset i in base coords */
U8*      seqRevCompPacked (U8 *u, U8 *rc, U64 len) ; /* reverse complements 2-bit packed binary */
U64      seqMatchPacked (U8 *a, U64 ia, U8 *b, U64 ib, U64 len) ;
	/* returns 0 if match, index+1 of first mismatching site if mismatch */

/* QualPack is similar for 1-bit qualities, mapping q < qualThresh to 0, q >= qualThresh to 1 */

typedef struct {
  int  qualThresh ;
  U64  qualExpand[256] ;
} QualPack ;

QualPack *qualPackCreate (int qualThresh) ;
#define   qualPackDestroy(sp) free(sp)
U8*       qualPack (QualPack *qp, char *q, U8 *u, U64 len) ; /* compress q into (len+7)/8 u */
char*     qualUnpack (QualPack *qp, U8 *u, char *q, U64 len) ; /* uncompress (len+7)/8 u into q */

/* SeqIO for flexible and efficient IO - I implement my own buffering */

typedef enum { UNKNOWN, FASTA, FASTQ, BINARY, ONE, BAM } SeqIOtype ;
extern char* seqIOtypeName[] ; // = { "unknown", "fasta", "fastq", "binary", "onecode", "bam" } ;

typedef struct {
  SeqIOtype type ;
  bool  isWrite ;
  U64   nSeq, totIdLen, totDescLen, totSeqLen, maxIdLen, maxDescLen, maxSeqLen ;
  U64   idLen, descLen, seqLen ;
  U64   idStart, descStart, seqStart, qualStart ;
  bool  isQual ;       		/* if set then convert qualities by subtracting 33 (FASTQ) */
  int   qualThresh ;		/* used for binary representation of qualities */
  /* below here private */
  U64   bufSize ;
  U64   nb ;			/* nb is how many characters left to read in the buffer */
  U64   line, recStart ;	/* recStart is the offset for the current record */
  int   fd ;			/* file descriptor, if gzf is not set */
  gzFile gzf ;
  char *buf, *b ;		/* b is current pointer in buf */
  int  *convert ;
  char *seqBuf, *qualBuf ;	/* used in modes BINARY, VGP, BAM */
  void *handle;			/* used for ONEseq, BAM */
  SeqPack  *seqPack ;
  QualPack *qualPack ;
} SeqIO ;

/* Reads/writes FASTA or FASTQ, gzipped or not, ONEseq, SAM/BAM/CRAM and a custom packed binary. */
/* Philosophy here is to read blocks of 8Mb and provide direct access into the buffer. */
/* So the user does not own the pointers. */
/* Add 0 terminators to ids.  Convert sequences in place if convert != 0, and quals if isQual. */

SeqIO  *seqIOopenRead (char *filename, int* convert, bool isQual) ; /* can use "-" for stdin */
bool    seqIOread (SeqIO *si) ;
#define sqioId(si)   ((si)->buf+(si)->idStart)
#define sqioDesc(si) ((si)->buf+(si)->descStart)
#define sqioSeq(si)  ((si)->type >= BINARY ? (si)->seqBuf : (si)->buf+(si)->seqStart)
#define sqioQual(si) ((si)->type >= BINARY ? (si)->qualBuf : (si)->buf+(si)->qualStart)

void    seqIOreferenceFileName (char *refFileName) ; /* resets this (globally) for CRAM */

SeqIO  *seqIOopenWrite (char *filename, SeqIOtype type, int* convert, int qualThresh) ;
void    seqIOwrite (SeqIO *si, char *id, char *desc, U64 seqLen, char *seq, char *qual) ;
void    seqIOflush (SeqIO *si) ;	/* NB writes are buffered, so need this to ensure in file */

void    seqIOclose (SeqIO *si) ;	/* will flush file opened for writing */

/* For ONEcode files, instead of seqio opening the file you can pass the OneCode handle. */
/* The handle must have primary type seq and support at least this schema (it can contain more). */

static char *seqioSchemaText =
  "1 3 def 1 0  schema for seqio\n"
  ".\n"
  "P 3 seq SEQUENCE\n"
  "O S 1 3 DNA               sequence: the DNA string\n"
  "D I 1 6 STRING            id: (optional) sequence identifier\n"
  "D Q 1 6 STRING            quality: Q values (ascii string = q+33)\n"
  "D N 3 3 INT 4 CHAR 3 INT  non-acgt base: pos (0-indexed), base, number\n" ;

SeqIO *seqIOadoptOneFile (void *handle, int* convert, int qualThresh) ;

/* utility */

char* seqRevComp (char *s, U64 len) ; /* reverse complements both index and text encodings */

/* standard converters - instantiated in seqio.c */

extern int dna2textConv[] ;
extern int dna2textAmbigConv[] ;
extern int dna2textAmbig2NConv[] ;
extern int dna2indexConv[] ;
extern int dna2index4Conv[] ;
extern int dna2binaryConv[] ;
extern int dna2binaryAmbigConv[] ;
static const char index2char[] = "acgtn" ;
static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;
extern int aa2textConv[] ;
extern int aa2indexConv[] ;
static const char index2aa[] = "ACDEFGHIKLMNPQRSTVWYX*" ;
extern int noConv[] ;
extern int complementBase[] ; /* complements both index and text with ambiguity, not binary */
extern int acgtCheck[] ;      /* 1 if acgtACGT, else 0 */

#endif	/* SEQIO_DEFINED */

/******************************************************************/
