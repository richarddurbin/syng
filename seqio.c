/*  File: seqio.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: buffered package to read arbitrary sequence files - much faster than readseq
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  6 11:22 2025 (rd109)
 * * Dec 15 09:45 2022 (rd109): separated out 2bit packing/unpacking into SeqPack
 * Created: Fri Nov  9 00:21:21 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <fcntl.h>
#include <unistd.h>

#ifdef ONEIO
#include "ONElib.h"
#endif

#ifdef BAMIO
bool bamFileOpenRead (char* filename, SeqIO *si) ;
bool bamRead (SeqIO *si) ;
void bamFileClose (SeqIO *si) ;
#endif

// global
char* seqIOtypeName[] = { "unknown", "fasta", "fastq", "binary", "onecode", "bam" } ;

SeqIO *seqIOopenRead (char *filename, int* convert, bool isQual)
{
  SeqIO *si = new0 (1, SeqIO) ;
  if (!strcmp (filename, "-")) si->gzf = gzdopen (fileno (stdin), "r") ;
  else si->gzf = gzopen (filename, "r") ;
  if (!si->gzf) { free(si) ; return 0 ; }
  si->bufSize = 1<<24 ; // 16 MB
  si->b = si->buf = new (si->bufSize, char) ;
  si->convert = convert ;
  si->isQual = isQual ;
  si->nb = gzread (si->gzf, si->buf, si->bufSize) ;
  if (!si->nb)
    { fprintf (stderr, "sequence file %s unreadable or empty\n", filename) ;
      seqIOclose (si) ;
      return 0 ;
    }
  si->line = 1 ;
  if (*si->buf == '>')
    { si->type = FASTA ; si->isQual = false ;
      if (!si->convert) si->convert = dna2textAmbigConv ; /* default: need to remove whitespace */
    }
  else if (*si->buf == '@')
    {
#ifdef BAMIO // problem: SAM file headers can start with this...
      if (si->buf[3] == '\t' &&
	  ( (si->buf[1] == 'H' && si->buf[2] == 'D') ||
	    (si->buf[1] == 'S' && si->buf[2] == 'Q') ||
	    (si->buf[1] == 'R' && si->buf[2] == 'G') ||
	    (si->buf[1] == 'P' && si->buf[2] == 'G') ||
	    (si->buf[1] == 'C' && si->buf[2] == 'O'))) // then almost certainly a SAM file
	{ gzclose (si->gzf) ; si->gzf = 0 ;
	  if (!bamFileOpenRead (filename, si))
	    { fprintf (stderr, "failed to open file %s as SAM/BAM/CRAM\n", filename) ;
	      seqIOclose (si) ;
	      return 0 ;
	    }
	  si->type = BAM ; // important that this is after successful open
	}
      else
	si->type = FASTQ ;
#else
      si->type = FASTQ ;
#endif
    }
  else if (*si->buf == 'b')
    { si->type = BINARY ;
      if (!si->convert) si->convert = dna2textConv ;
      if (si->nb <= 64) die ("binary file too short\n") ;
      si->qualThresh = si->buf[1] ; si->b += 8 ;
      si->nSeq = *(U64*)si->b ; si->b += 8 ;
      si->totIdLen = *(U64*)si->b ; si->b += 8 ;
      si->totDescLen = *(U64*)si->b ; si->b += 8 ;
      si->totSeqLen = *(U64*)si->b ; si->b += 8 ;
      si->maxIdLen = *(U64*)si->b ; si->b += 8 ;
      si->maxDescLen = *(U64*)si->b ; si->b += 8 ;
      si->maxSeqLen = *(U64*)si->b ; si->b += 8 ;
      si->nb -= 64 ;
      si->seqPack = seqPackCreate (si->convert['a']) ;
      si->qualPack = qualPackCreate (si->qualThresh) ;
      si->seqBuf = new0 (si->maxSeqLen+1, char) ;
      if (si->isQual) si->qualBuf = new0 (si->maxSeqLen+1, char) ;
      { U64 maxBufSize = 3*sizeof(int) + 5 + si->maxIdLen + si->maxDescLen + si->maxSeqLen / 4 ;
	if (si->qualThresh) maxBufSize += (si->maxSeqLen / 8) ; /* 5 = 1 + 1 + 3 for pad */
	if (maxBufSize > si->nb)
	  { maxBufSize = ((maxBufSize >> 20) + 1) << 20 ; /* so a clean number of megabytes */
	    char *newBuf = new (maxBufSize, char) ; memcpy (newBuf, si->b, si->nb) ;
	    si->b = si->buf = newBuf ; si->bufSize = maxBufSize ;
	    si->nb += gzread (si->gzf, si->b + si->nb, si->bufSize - si->nb) ;
	  }
      }
    }
#ifdef ONEIO
  else if (*si->buf == '1')
    { gzclose (si->gzf) ; si->gzf = 0 ;
      OneFile *vf = oneFileOpenRead (filename, 0, "seq", 1) ;
      if (!vf)
	{ fprintf (stderr, "failed to open ONE seq file %s\n", filename) ;
	  seqIOclose (si) ;
	  return 0 ;
	}
      si->type = ONE ; // important that this is after successful open
      si->handle = vf ;
      if (vf->info['S']->given.count)
	{ si->nSeq = vf->info['S']->given.count ;
	  si->totSeqLen = vf->info['S']->given.total ;
	  si->maxSeqLen = vf->info['S']->given.max ;
	  si->seqBuf = new0 (si->maxSeqLen+1,char) ;
	}
      if (si->isQual && vf->info['Q'] && vf->info['Q']->given.count)
	si->qualBuf = new0 (si->maxSeqLen+1,char) ;
      else
	si->isQual = false ;
      while (oneReadLine (vf) && vf->lineType != 'S') ; // move up to first sequence line
      si->seqStart = 0 ;
    }
#endif
#ifdef BAMIO
  else
    { gzclose (si->gzf) ; si->gzf = 0 ;
      if (!bamFileOpenRead (filename, si))
	{ fprintf (stderr, "failed to open file %s as SAM/BAM/CRAM\n", filename) ;
	  seqIOclose (si) ;
	  return 0 ;
	}
      si->type = BAM ; // important that this is after successful open
    }
#else
  else
    { fprintf (stderr, "sequence file %s is unknown type\n", filename) ;
      seqIOclose (si) ;
      return 0 ;
    }
#endif
  return si ;
}

void seqIOclose (SeqIO *si)
{ if (si->isWrite)
    { if (si->type <= BINARY)
	seqIOflush (si) ;
      if (si->type == BINARY)	/* write header */
	{ if (lseek (si->fd, 0, SEEK_SET)) die ("failed to seek to start of binary file") ;
	  si->b = si->buf; 
	  *si->b++ = 'b' ; *si->b++ = si->qualThresh ; si->b += 6 ;
	  *(U64*)si->b = si->nSeq  ; si->b += 8 ;
	  *(U64*)si->b = si->totIdLen ; si->b += 8 ;
	  *(U64*)si->b = si->totDescLen ; si->b += 8 ;
	  *(U64*)si->b = si->totSeqLen ; si->b += 8 ;
	  *(U64*)si->b = si->maxIdLen ; si->b += 8 ;
	  *(U64*)si->b = si->maxDescLen ; si->b += 8 ;
	  *(U64*)si->b = si->maxSeqLen ; si->b += 8 ;
	  seqIOflush (si) ;
	}
    }
  free (si->buf) ;
  if (si->seqBuf) free (si->seqBuf) ;
  if (si->qualBuf) free (si->qualBuf) ;
  if (si->gzf) gzclose (si->gzf) ;
  if (si->fd) close (si->fd) ;
#ifdef ONEIO
  if (si->type == ONE)
    oneFileClose ((OneFile*)si->handle) ;
#endif
#ifdef BAMIO
  if (si->type == BAM)
    bamFileClose (si) ;
#endif
  free (si) ;
}

/********** local routines for seqIOread() ***********/
 
static void bufRefill (SeqIO *si)
{ si->b -= si->recStart ;		/* will be position after move */
  memmove (si->buf, si->buf + si->recStart, si->b - si->buf) ;
  si->idStart -= si->recStart ; si->descStart -= si->recStart ; /* adjust all the offsets */
  si->seqStart -= si->recStart ; si->qualStart -= si->recStart ;
  si->recStart = 0 ;
  si->nb = gzread (si->gzf, si->b, si->buf + si->bufSize - si->b) ;
}

static void bufDouble (SeqIO *si)
{
  char *newbuf = new (si->bufSize*2, char) ;
  memcpy (newbuf, si->buf, si->bufSize) ;
  si->b = newbuf + si->bufSize ; si->nb = si->bufSize ; /* rely on being at end of old buf */
  free (si->buf) ; si->buf = newbuf ;
  si->nb = gzread (si->gzf, si->b, si->bufSize) ;
  si->bufSize *= 2 ;
}

#define bufAdvanceEndRecord(si) \
  { ++si->b ; \
    if (!--si->nb) \
      { if (si->recStart) bufRefill (si) ; \
        else if (si->b == si->buf + si->bufSize) bufDouble (si) ; \
      } \
  }

#define bufAdvanceInRecord(si) \
  { bufAdvanceEndRecord(si) ; \
    if (!si->nb) \
      { fprintf (stderr, "incomplete sequence record line %llu\n", si->line) ; return false ; } \
  } 

static void bufHardRefill (SeqIO *si, U64 n) /* like bufRefill() but for bufConfirmNbytes() */
{					     /* NB buf should be big enough because of header */
  si->b -= si->recStart ;		/* will be position after move */
  memmove (si->buf, si->buf + si->recStart, si->b - si->buf) ;
  si->recStart = 0 ; si->b = si->buf ;
  si->nb += gzread (si->gzf, si->b + si->nb, si->bufSize - si->nb) ;
  if (si->nb < n) die ("incomplete sequence record %llu", si->line) ;
}

#define bufConfirmNbytes(si, n) { if (si->nb < n) bufHardRefill (si, n) ; }

#include <ctype.h>

bool seqIOread (SeqIO *si)
{
#ifdef ONEIO
  if (si->type == ONE)
    { OneFile *vf = (OneFile*) si->handle ;
      if (vf->lineType != 'S') return false ; // at end of file
      si->seqLen = oneLen(vf) ; // otherwise we are at an 'S' line
      if (si->seqLen > si->maxSeqLen)
	{ if (si->maxSeqLen) { free (si->seqBuf) ; if (si->isQual) free (si->qualBuf) ; }
	  si->maxSeqLen = si->seqLen ;
	  si->seqBuf = new0 (si->maxSeqLen+1, char) ;
	  if (si->isQual) si->qualBuf = new0 (si->maxSeqLen+1, char) ;
	}
      if (si->convert)
	{ char *s = si->seqBuf, *e = s + si->seqLen, *sv = oneString(vf) ;
	  while (s < e) *s++ = si->convert[(int)*sv++] ;
	}
      else
	memcpy (si->seqBuf, oneString(vf), si->seqLen) ;
      si->idLen = 0 ; si->idStart = 0 ;
      si->descLen = 0 ; si->descStart = 0 ;
      *si->buf = 0 ;
      while (oneReadLine (vf) && vf->lineType != 'S')
	{ if (vf->lineType == 'Q' && si->isQual)
	    { char *q = si->qualBuf, *e = q + si->seqLen, *qv = oneString(vf) ;
	      while (q < e) *q++ = *qv++ - 33 ;
	    }
	  if (vf->lineType == 'I')
	    { char *desc = oneReadComment (vf) ;
	      si->idLen = oneLen(vf) ;
	      if (desc) si->descLen = strlen(desc) ; else si->descLen = 0 ;
	      while (si->idLen + si->descLen + 2 > si->bufSize) bufDouble (si) ;
	      si->idStart = 0 ; strcpy (si->buf, oneString(vf)) ;
	      si->descStart = si->idLen+1 ;
	      if (desc) strcpy (si->buf+si->descStart, desc) ; else si->buf[si->descStart] = 0 ;
	    }
	  if (vf->lineType == 'N')
	    { int n = oneInt(vf,2) ;
	      char base = si->convert ? si->convert[(int)oneChar(vf,1)] : oneChar(vf,1) ;
	      while (n--) si->seqBuf[oneInt(vf,0)+n] = base ;
	    }
	}
      return true ;
    }
#endif
#ifdef BAMIO
  if (si->type == BAM) return bamRead (si) ;
#endif

  if (!si->nb) return false ;
  si->recStart = si->b - si->buf ;
  
  if (si->type == BINARY)
    { if (si->line > si->nSeq) return false ; /* have already read all the sequences */
      bufConfirmNbytes (si, (U64)(3*sizeof(int))) ;
      si->idLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->descLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->seqLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->nb -= 3*sizeof(int) ;
      U64 nBytes = si->idLen + 1 + si->descLen + 1 + (si->seqLen+3)/4 ;
      if (si->qualThresh) nBytes += (si->seqLen+7)/8 ; /* NB qualThresh (from file) not isQual */
      nBytes = 4*((nBytes+3)/4) ;
      bufConfirmNbytes (si, nBytes) ;
      si->idStart = si->b - si->buf ;
      si->descStart = si->idStart + si->idLen + 1 ;
      si->seqStart = si->descStart + si->descLen + 1 ;
      seqUnpack (si->seqPack, (U8*)(si->buf+si->seqStart), si->seqBuf, 0, si->seqLen) ;
      if (si->isQual)
	{ si->qualStart = si->seqStart + (si->seqLen + 3) / 4 ;
	  qualUnpack (si->qualPack, (U8*)(si->buf+si->qualStart), si->qualBuf, si->seqLen) ;
	}
      ++si->line ;
      si->b += nBytes ; si->nb -= nBytes ;
      return true ;
    }
  
  /* if get to here then this is a text file, FASTA or FASTQ */
  
  if (si->type == FASTA)
    { if (*si->b != '>') die ("no initial > for FASTA record line %llu", si->line) ; }
  else if (si->type == FASTQ)
    { if (*si->b != '@') die ("no initial @ for FASTQ record line %llu", si->line) ; }
  bufAdvanceInRecord(si) ; si->idStart = si->b - si->buf ;
  while (!isspace(*si->b)) bufAdvanceInRecord(si) ;
  si->idLen = si->b - sqioId(si) ;
  if (*si->b != '\n') /* a space or tab - whatever follows on this line is description */
    { *si->b = 0 ; bufAdvanceInRecord(si) ;
      si->descStart = si->b - si->buf ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ;
      si->descLen = si->b - sqioDesc(si) ;
    }
  else { si->descLen = si->descStart = 0 ; }
  *si->b = 0 ;
  ++si->line ; bufAdvanceInRecord(si) ;	              /* line 2 */
  si->seqStart = si->b - si->buf ;
  if (si->type == FASTA)
    { while (si->nb && *si->b != '>')
	{ while (*si->b != '\n') bufAdvanceInRecord(si) ;
	  ++si->line ; bufAdvanceEndRecord(si) ;
	}
      char *s = sqioSeq(si), *t = s ;
      while (s < si->b) if ((*t++ = si->convert[(int)*s++]) < 0) --t ;
      si->seqLen = t - sqioSeq(si) ;
    }
  else if (si->type == FASTQ)
    { while (*si->b != '\n') bufAdvanceInRecord(si) ;
      si->seqLen = si->b - sqioSeq(si) ;
      if (si->convert)
	{ char *s = sqioSeq(si) ;
	  while (s < si->b) { *s = si->convert[(int)*s] ; ++s ; }
	}
      ++si->line ; bufAdvanceInRecord(si) ; 	      /* line 3 */
      if (*si->b != '+') die ("missing + FASTQ line %llu", si->line) ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ; /* ignore remainder of + line */
      ++si->line ; bufAdvanceInRecord(si) ;	      /* line 4 */
      si->qualStart = si->b - si->buf ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ;
      if (si->b - si->buf - si->qualStart != si->seqLen)
	die ("qual not same length as seq line %llu", si->line) ;
      if (si->isQual) { char *q = sqioQual(si), *e = q + si->seqLen ; while (q < e) *q++ -= 33 ; }
      ++si->line ; bufAdvanceEndRecord(si) ;
    }

  ++si->nSeq ;
  return true ;
}

/*********************** open for writing ***********************/

SeqIO *seqIOopenWrite (char *filename, SeqIOtype type, int* convert, int qualThresh)
{
  SeqIO *si = new0 (1, SeqIO) ;

  int nameLen = strlen (filename) ;
  char *ext ; // extension
  for (ext = filename + nameLen ; ext > filename && *--ext != '.' ;) { ; } ++ext ;
  bool isGzip = !strcmp (ext, "gz") ;
  
  si->isWrite = true ;
  si->convert = convert ;
  si->isQual = (qualThresh > 0) ;
  si->qualThresh = qualThresh ;
  if (si->type == FASTA && si->isQual)
    { fprintf (stderr, "warning : can't write qualities to FASTA file %s\n", filename) ;
      si->isQual = false ;
    }
  if (type)
    si->type = type ;
  else
    { if (isGzip)
	{ if ((nameLen > 9 && !strcmp (filename+nameLen-9, ".fasta.gz")) ||
	      (nameLen > 6 && !strcmp (filename+nameLen-6, ".fa.gz"))) si->type = FASTA ;
	  if ((nameLen > 9 && !strcmp (filename+nameLen-9, ".fastq.gz")) ||
	      (nameLen > 6 && !strcmp (filename+nameLen-6, ".fq.gz"))) si->type = FASTQ ;
	}
      else if (!strcmp (ext, "fasta") || !strcmp (ext, "fa")) si->type = FASTA ;
      else if (!strcmp (ext, "fastq") || !strcmp (ext, "fq")) si->type = FASTQ ;
      else if (*ext == '1') si->type = ONE ;
      else { warn ("can't determine sequence file type for %s", filename) ; free(si) ; return 0 ; }
    }
  
  if (si->type == ONE) 
#ifdef ONEIO
    { OneSchema *schema = oneSchemaCreateFromText (seqioSchemaText) ;
      char *oneType = "seq" ;
      //      if (*ext == '1') ++ext ; oneType = ext ; // don't think I want this
      OneFile *vf = oneFileOpenWriteNew (filename, schema, oneType, true, 1) ;
      oneSchemaDestroy (schema) ;
      if (!vf) { free (si) ; return 0 ; }
      si->handle = vf ;
      char *commandLine = getCommandLine() ;
      if (!commandLine) commandLine = "-" ;
      oneAddProvenance (vf, "seqio", "1.0", commandLine) ;
      return si ;
    }
#else
    { warn ("sorry, seqio not compiled with ONElib so can't write 1seq") ;
      free (si) ; return 0 ;
    }
#endif

  if (si->type == BAM)
    { warn ("sorry, seqio can't write BAM") ;
      free (si) ; return 0 ;
    }
  
  if (!strcmp (filename, "-"))
    { si->fd = fileno (stdout) ;
      if (si->fd == -1) { warn ("failed to write to stdout") ; free (si) ; return 0 ; }
    }
  else if (!strcmp (filename, "-z"))
    { si->gzf = gzdopen (fileno(stdout), "w") ;
      if (!si->gzf) { warn ("failed to gzdopen stdout", filename) ; free (si) ; return 0 ; }
    }
  else if (isGzip)
    { si->gzf = gzopen (filename, "w") ;
      if (!si->gzf) { warn ("failed to gzopen %s", filename) ; free (si) ; return 0 ; }
    }
  else
    { si->fd = open (filename, O_CREAT | O_TRUNC | O_WRONLY, 00644) ;
      if (si->fd == -1) { warn ("failed to open %s", filename) ; free (si) ; return 0 ; }
    }
  if (si->type == BINARY && si->gzf)
    { fprintf (stderr, "can't write a gzipped binary file\n") ; seqIOclose (si) ; return 0 ; }
  
  si->nb = si->bufSize = 1<<24 ;
  si->b = si->buf = new (si->bufSize, char) ;

  if (si->type == BINARY)
    { si->b += 64 ; si->nb -= 64 ; /* make space for header - written on file close */
      if (!si->convert) si->convert = dna2textConv ;
      si->seqPack = seqPackCreate (si->convert['a']) ;
      si->qualPack = qualPackCreate (si->qualThresh) ;
    }

  return si ;
}

SeqIO *seqIOadoptOneFile (void *handle, int* convert, int qualThresh)
{
  if (!handle) return 0 ;
  
 #ifndef ONEIO
  warn ("sorry, seqio not compiled with ONElib so can't write 1seq") ;
  return 0 ;
#else
  if (!oneFileCheckSchemaText ((OneFile*)handle, seqioSchemaText)) return 0 ;
  
  SeqIO *si = new0 (1, SeqIO) ;
  si->isWrite = true ;
  si->convert = convert ;
  si->isQual = (qualThresh > 0) ;
  si->qualThresh = qualThresh ;
  si->type = ONE ;
  si->handle = handle ;
  return si ;
#endif
}

void seqIOflush (SeqIO *si)	/* writes buffer to file and resets to 0 */
{
  if (!si->isWrite) return ;
  U64 retVal, nBytes = si->b - si->buf ;
  if (si->gzf) retVal = gzwrite (si->gzf, si->buf, nBytes) ;
  else retVal = write (si->fd, si->buf, nBytes) ;
  if (retVal != nBytes) die ("seqio write error %llu not %llu bytes written", retVal, nBytes) ;
  si->b = si->buf ;
  si->nb = si->bufSize ;
}

static void writeExtend (SeqIO *si, U64 len)
{
  while (si->bufSize < len) si->bufSize <<= 1 ;
  free (si->buf) ;
  si->b = si->buf = new (si->bufSize, char) ;
  si->nb = si->bufSize ;
}
  
void seqIOwrite (SeqIO *si, char *id, char *desc, U64 seqLen, char *seq, char *qual)
{
  U64 len, pad = 0 ;
  assert (si->isWrite) ;

  ++si->nSeq ;

  { static char buf[24] ;
    if (si->type != ONE && !id) { id = buf ; sprintf (buf, "%lld", si->nSeq) ; }
  }

  si->idLen = id ? strlen(id) : 0 ;
  si->totIdLen += si->idLen ; if (si->idLen > si->maxIdLen) si->maxIdLen = si->idLen ;
  si->descLen = desc ? strlen(desc) : 0 ;
  si->totDescLen += si->descLen ; if (si->descLen > si->maxDescLen) si->maxDescLen = si->descLen ;
  si->seqLen = seqLen ;
  si->totSeqLen += si->seqLen ; if (si->seqLen > si->maxSeqLen) si->maxSeqLen = si->seqLen ;

#ifdef ONEIO
  if (si->type == ONE)
    { OneFile *vf = (OneFile*)(si->handle) ;
      static U64 bufLen = 0 ;
      static char *buf = 0 ;
      I64 i ;
      if (seqLen > bufLen)
	{ if (buf) free (buf) ;
	  bufLen = seqLen ;
	  buf = new(bufLen+1,char) ;
	}
      if (si->convert)
	{ for (i = 0 ; i < seqLen ; ++i) buf[i] = si->convert[(int)(seq[i])] ;
	  seq = buf ;
	}
      oneWriteLine (vf, 'S', seqLen, seq) ;
      if (id)
	{ oneWriteLine (vf, 'I', strlen (id), id) ;
	  if (desc) oneWriteComment (vf, "%s", desc) ;
	}
      if (qual && si->isQual)
	{ for (i = 0 ; i < seqLen ; ++i) buf[i] = qual[i] + 33 ;
	  oneWriteLine (vf, 'Q', seqLen, buf) ;
	}
      for (i = 0 ; i < seqLen ; ++i) // write exceptions for non-ACGT characters
	if (!acgtCheck[(int)seq[i]])
	  { oneInt(vf,0) = i ;
	    char c = oneChar(vf,1) = seq[i] ;
	    I64 n = 1 ;
	    while (++i < seqLen && seq[i] == c) n++ ;
	    oneInt(vf,2) = n ;
	    oneWriteLine (vf, 'N', 0, 0) ;
	  }
      return ;
    }
#endif
  
  if (si->type == FASTA)
    len = 3 + si->idLen + (desc ? 1+si->descLen : 0) + seqLen ;   /* >id[ desc]\nseq\n */
  else if (si->type == FASTQ)
    len = 6 + si->idLen + (desc ? 1+si->descLen : 0) + 2*seqLen ; /* @id[ desc]\nseq\n+\nqual\n */
  else	/* binary */
    { U64 nBytes = si->idLen + si->descLen + 2 + (si->seqLen+3)/4 ;
      if (si->isQual) nBytes += (si->seqLen+7)/8 ;
      pad = 3 - ((nBytes+3) % 4) ; /* 0,1,2,3 maps to 0,3,2,1 */
      len = 3*sizeof(int) + nBytes + pad ;
    }
  if (len > si->nb) seqIOflush (si) ;
  if (len > si->nb) writeExtend (si, len) ;

  if (si->type == FASTA)
    { *si->b++ = '>' ;
      if (id) { strcpy (si->b, id) ; si->b += si->idLen ; }
      if (desc) { *si->b++ = ' ' ; strcpy (si->b, desc) ; si->b += si->descLen ; }
      *si->b++ = '\n' ;
      memcpy (si->b, seq, seqLen) ;
      if (si->convert) while (seqLen--) { *si->b = si->convert[(int)*si->b] ; ++si->b ; }
      else si->b += seqLen ;
      *si->b++ = '\n' ;
    }
  else if (si->type == FASTQ)
    { *si->b++ = '@' ;
      if (id) { strcpy (si->b, id) ; si->b += si->idLen ; }
      if (desc) { *si->b++ = ' ' ; strcpy (si->b, desc) ; si->b += si->descLen ; }
      *si->b++ = '\n' ;
      memcpy (si->b, seq, seqLen) ;
      if (si->convert) while (seqLen--) { *si->b = si->convert[(int)*si->b] ; ++si->b ; }
      else si->b += seqLen ;
      *si->b++ = '\n' ;
      *si->b++ = '+' ;
      *si->b++ = '\n' ;
      seqLen = si->seqLen ; while (seqLen--) *si->b++ = qual ? *qual++ + 33 : 33 ;
      *si->b++ = '\n' ;
    }
  else				/* binary */
    { int *ib = (int*)si->b ; si->b += 3 * sizeof(int) ;
      *ib++ = si->idLen ; *ib++ = si->descLen ; *ib++ = seqLen ;
      if (si->idLen) { strcpy (si->b, id) ; si->b += si->idLen ; } *si->b++ = 0 ;
      if (si->descLen) { strcpy (si->b, desc) ; si->b += si->descLen ; } *si->b++ = 0 ;
      seqPack (si->seqPack, seq, (U8*)si->b, si->seqLen) ; si->b += (si->seqLen + 3)/4 ;
      if (si->isQual)
	{ qualPack (si->qualPack, qual, (U8*)si->b, si->seqLen) ;
	  si->b += (si->seqLen + 7)/8 ;
	}
      si->b += pad ;
    }
  si->nb -= len ;
}

/********* reverse complement sequences ********/

char* seqRevComp (char* s, U64 len) // index and text (including ambig) 
{
  char *r = new(len,char) ;
  r += len ;
  while (len--) *--r = complementBase[(int)*s++] ;
  return r ;
}

/********** some routines to pack sequence and qualities for binary representation ***********/

SeqPack *seqPackCreate (char unpackA)
{
  SeqPack *sp = new (1, SeqPack) ;
  switch (unpackA)
    {
    case 0:
      sp->unconv[0] = 0; sp->unconv[1] = 1; sp->unconv[2] = 2; sp->unconv[3] = 3;
      sp->unconvC[0] = 3; sp->unconvC[1] = 2; sp->unconvC[2] = 1; sp->unconv[3] = 0;
      break ;
    case 1:
      sp->unconv[0] = 1; sp->unconv[1] = 2; sp->unconv[2] = 4; sp->unconv[3] = 8;
      sp->unconvC[0] = 8; sp->unconvC[1] = 4; sp->unconv[2] = 2; sp->unconv[3] = 1;
      break ;
    case 'a': strcpy (sp->unconv, "acgt") ; strcpy (sp->unconvC, "tgca") ; break ;
    case 'A': strcpy (sp->unconv, "ACGT") ; strcpy (sp->unconvC, "TGCA") ; break ;
    default: 
      die ("seqPackCreate: unrecognised unpackA character %d = %c - must be one of a, A, 0, 1",
	   unpackA, unpackA) ;
    }
  int i, j ;
  for (i = 0 ; i < 256 ; ++i)
    { U8 u ;
      char *s = (char*)&sp->byteExpand[i] ;
      u = i ; for (j = 0 ; j < 4 ; ++j) { s[j] = sp->unconv[u & 0x03] ; u >>= 2 ; }
      s = (char*)&sp->byteExpandC[i] ;
      u = i ; for (j = 4 ; j-- ; ) { s[j] = sp->unconvC[u & 0x03] ; u >>= 2 ; }
    }
  return sp ;
}

static U8 pack[] = {    // sends N (indeed any non-CGT) to A, except 0,1,2,3 are maintained
   0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
} ;

static U8 packC[] = {   // same but send to the complement
   3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
} ;

U8* seqPack (SeqPack *sp, char *s, U8 *u, U64 len) /* compress s into (len+3)/4 u */
{
  if (!u) u = new((len+3)/4,U8) ;
  U8 *u0 = u ;
  while (len >= 4)
    { *u++ = pack[(int)s[0]] | (pack[(int)s[1]] << 2) |
	(pack[(int)s[2]] << 4) | (pack[(int)s[3]] << 6) ;
      len -= 4 ; s += 4 ;
    }
  switch (len)
    { case 3: *u++ = pack[(int)s[0]] | (pack[(int)s[1]] << 2) |	(pack[(int)s[2]] << 4) ; break ;
      case 2: *u++ = pack[(int)s[0]] | (pack[(int)s[1]] << 2) ; break ;
      case 1: *u++ = pack[(int)s[0]] ; break ;
      case 0: break ;
    }
  return u0 ;
}

U8* seqPackRevComp (SeqPack *sp, char *s, U8 *u, U64 len) /* packs the RC of the sequence s */
{
  if (!u) u = new((len+3)/4,U8) ;
  U8 *u0 = u ;
  s += len-4 ;
  while (len >= 4)
    { *u++ = packC[(int)s[3]] | (packC[(int)s[2]] << 2) |
	(packC[(int)s[1]] << 4) | (packC[(int)s[0]] << 6) ;
      len -= 4 ; s-= 4 ;
    }
  switch (len)
    { case 3: *u++ = packC[(int)s[3]] | (packC[(int)s[2]] << 2) | (packC[(int)s[1]] << 4) ; break ;
      case 2: *u++ = packC[(int)s[3]] | (packC[(int)s[2]] << 2) ; break ;
      case 1: *u++ = packC[(int)s[3]] ; break ;
      default: break ;
    }
  return u0 ;
}

char* seqUnpack (SeqPack *sp, U8 *u, char *s, U64 i, U64 len)
{
  if (!s) s = (*sp->unconv >= 'A') ? new0 (len+1,char) : new(len,char) ;
  char *s0 = s ;
  u += (i >> 2) ; i = (i & 3) ; // first offset and go to byte (U8) boundary
  if (i)
    { U8 uu = *u >> 2*i ;
      while ((i++ < 4) && len) { *s++ = sp->unconv[uu & 3] ; uu >>= 2 ; len-- ; }
      ++u ;
    }
  if (len)
    { while (len >= 4)
	{ *(U32*)s = sp->byteExpand[*u] ;
	  ++u ; s += 4 ; len -= 4 ;
	}
      switch (len) // note fall through
	{
	case 3: s[2] = sp->unconv[(*u >> 4) & 3] ;
	case 2: s[1] = sp->unconv[(*u >> 2) & 3] ;
	case 1: s[0] = sp->unconv[*u & 3] ; 
	case 0: break ;
	}
    }
  return s0 ;
}

char* seqUnpackRevComp (SeqPack *sp, U8 *u, char *s, U64 i, U64 len)
{
  if (!s) s = (*sp->unconv >= 'A') ? new0 (len+1,char) : new(len,char) ;
  char *s0 = s ;
  s += len ;
  u += (i >> 2) ; i = (i & 3) ; // first offset and go to byte (U8) boundary
  if (i)
    { U8 uu = *u >> 2*i ;
      while ((i++ < 4) && len) { *--s = sp->unconvC[uu & 3] ; uu >>= 2 ; len-- ; }
      ++u ;
    }
  if (len)
    { while (len >= 4)
	{ s -= 4 ;
	  *(U32*)s = sp->byteExpandC[*u] ;
	  ++u ; len -= 4 ;
	}
      switch (len) // note fall through
	{
	case 3: s[-3] = sp->unconvC[(*u >> 4) & 3] ;
	case 2: s[-2] = sp->unconvC[(*u >> 2) & 3] ;
	case 1: s[-1] = sp->unconvC[*u & 3] ; 
	case 0: break ;
	}
    }
  return s0 ;
}

U8 *seqRevCompPacked (U8* u, U8 *rc, U64 len)
{
  if (!rc) rc = new(len,U8) ;
  U64 i ;
  static U8 rcByte[256] ;
  static bool isFirst = true ;
  if (isFirst)
    { for (i = 0 ; i < 256 ; ++i)
	rcByte[255-i] = ((i & 3) << 6) | ((i & 12) << 2) | ((i & 48) >> 2) | ((i & 192) >> 6) ;
      isFirst = false ;
    }

  U64 blen = (len+3)/4 ;
  if (!rc) rc = new((len+3)/4,U8) ;
  for (i = blen ; i-- ;) rc[i] = rcByte[*u++] ;
  switch (len % 4) // now shift back
    {
    case 0:
      break ;
    case 1:
      for (i = 0 ; i < blen-1 ; ++i) rc[i] = (rc[i] >> 6) | (rc[i+1] << 2) ;
      rc[blen-1] = rc[blen-1] >> 6 ;
      break ;
    case 2:
      for (i = 0 ; i < blen-1 ; ++i) rc[i] = (rc[i] >> 4) | (rc[i+1] << 4) ;
      rc[blen-1] = rc[blen-1] >> 4 ;
      break ;
    case 3:
      for (i = 0 ; i < blen-1 ; ++i) rc[i] = (rc[i] >> 2) | (rc[i+1] << 6) ;
      rc[blen-1] = rc[blen-1] >> 2 ;
      break ;
    }
  
  return rc ;
}

U64 seqMatchPacked2 (U8 *a, U64 ia, U8 *b, U64 ib, U64 len) /* returns 0 or distance to mismatch */
{
  static SeqPack *sp = 0 ;
  if (sp == 0) sp = seqPackCreate ('a') ;
  U64 i = 0 ;
  char *sa = seqUnpack (sp, a, 0, ia, len), *sb = seqUnpack (sp, b, 0, ib, len) ;
  while (i < len && sa[i] == sb[i]) i++ ;
  free (sa) ; free (sb) ;
  if (i == len) return 0 ; else return i+1 ;
}

U64 seqMatchPacked (U8 *a, U64 ia, U8 *b, U64 ib, U64 len)
{
  static SeqPack *sp ;
  static U64 mask[33] ; // need 1 to 32
  static bool isFirst = true ;
  int i, chunk = 0 ;
  if (isFirst) // initialise statics
    { sp = seqPackCreate ('a') ;
      union { U64 u ; U8 b[8] ; } z ;
      z.u = 0 ;
      for (i = 0 ; i < 8 ; i++)
	{ z.b[i] = 0x03 ; mask[4*i+1] = z.u ;
	  z.b[i] = 0x0f ; mask[4*i+2] = z.u ;
	  z.b[i] = 0x3f ; mask[4*i+3] = z.u ;
	  z.b[i] = 0xff ; mask[4*i+4] = z.u ;
	}
      isFirst = false ;
    }

  // advance sequences to the 64-bit words in which the start points begin
  // recall that a U64 holds 32 bases, (z & 0x1f) == z % 32, and (z >> 5) == z / 32
  U64 *ua = &(((U64*)a)[ia>>5]), *ub = &(((U64*)b)[ib>>5]) ;
  ia &= 0x1f ; ib &= 0x1f ;

  // debug
  int iia = ia, iib = ib, ilen = len ;
  char sa[16], sb[16] ; sa[15] = sb[15] = 0 ;
  int l = (len < 15) ? len : 15 ;
  a = (U8*) ua ; b = (U8*) ub ;
  seqUnpack (sp, a, sa, ia, l) ; seqUnpack (sp, b, sb, ib, l) ;
  U64 d2 = seqMatchPacked2 (a, ia, b, ib, len) ;
  char ca, cb ; if (d2) { seqUnpack (sp, a, &ca, d2-1, 1) ; seqUnpack (sp, b, &cb, d2-1, 1) ; }

  U64 uua = 0, uub = 0, d = 0 ; // distance of match
  if (len > 32) // can work in 64-bit words
    { if (ia == ib)
	{ if (((*ua++ >> (2*ia)) & mask[32-ia]) != ((*ub++ >> (2*ib)) & mask[32-ib]))
	    { d = 32 - ia ;
	      goto end ;
	    }
	  d += (32 - ia) ; len -= (32 - ia) ;
	  while (len > 32)
	    { d += 32 ; len -= 32 ;
	      uua = *ua ; uub = *ub ; // DEBUG
	      if ((*ua++ & mask[32]) != (*ub++ & mask[32])) goto end ;
	    }
	}
      else
	{ if (ib < ia) { U64 *t = ua ; ua = ub ; ub = t ; U64 it = ia ; ia = ib ; ib = it ; }
	  uua = *ua >> (2*ia) ; uub = *ub >> (2*ib) ;
	  while (len > 32)
	    { chunk = 32-ib ; d += chunk ; len -= chunk ;
	      if ((uua & mask[chunk]) != (uub & mask[chunk])) goto end ;
	      uua >>= 2*chunk ; ia += chunk ; uub = *++ub ; ib = 0 ;
	      chunk = 32-ia ; d += chunk ; len -= chunk ;
	      if ((uua & mask[chunk]) != (uub & mask[chunk])) goto end ;
	      uub >>= 2*chunk ; ib += chunk ; uua = *++ua ; ia = 0 ; 
	    }
	}
    }
  
  // final bases which may not fit into a 32 bit word
  char abuf[33], bbuf[33] ;
  seqUnpack (sp, (U8*)ua, abuf, ia, len) ; seqUnpack (sp, (U8*)ub, bbuf, ib, len) ;
  for (i = 0 ; i < len ; ++i)
    { ++d ;
      if (abuf[i] != bbuf[i])
	{ abuf[len] = bbuf[len] = 0 ;
	  if (d2 != d)
	    printf ("SEQMATCH ia %d ib %d len %d d %lld d2 %lld abuf %s bbuf %s i %d\n",
		    iia, iib, ilen, d, d2, abuf, bbuf, i) ;
	  return d ;
	}
    }
  if (d2)
    printf ("SEQMATCH ia %d ib %d len %d d %lld d2 %lld abuf %s bbuf %s i %d ca %c cb %c\n",
	    iia, iib, ilen, d, d2, abuf, bbuf, i, ca, cb) ;
  d = 0 ; // if we get here then they all match

 end:
  if ((!d2 && d) || (d != d2 && d+32 < d2))
    { printf ("seqMatch ia %d ib %d len %d d %lld d2 %lld a %s b %s\n",
	      iia, iib, ilen, d, d2, sa, sb) ;
      printf ("   local d %lld  uua %0llx uub %0llx chunk %d mask[chunk] %0llx\n",
	      d, uua, uub, chunk, mask[chunk]) ;
    }
  return d ;
}

/************ QualPack package ***************/

QualPack *qualPackCreate (int qualThresh)
{
  QualPack *qp = new (1, QualPack) ;
  if (qualThresh < 0 || qualThresh >= 128)
    die ("seqPackCreate: qualThresh %d is out of range 0..127", qualThresh) ;
  qp->qualThresh = qualThresh ;

  int i, j ;
  for (i = 0 ; i < 256 ; ++i)
    { U8 u ;
      char *q = (char*)&qp->qualExpand[i] ;
      u = i ; for (j = 0 ; j < 8  ; ++j) { q[j] = (u & 0x01) ? qualThresh : 0 ; u >>= 1 ; }
    }

  return qp ;
}

U8* qualPack (QualPack *sp, char *q, U8 *u, U64 len) /* compress q into (len+7)/8 u */
{
  if (!u) u = new((len+7)/8,U8) ;
  U8 *u0 = u ;
  int i ;
  while (len >= 8)
    { *u = 0 ; for (i = 0 ; i < 8 ; ++i) { if (*q++ >= sp->qualThresh) *u |= (1 << i) ; }
      len -= 8 ; ++u ;
    }
  if (len)
    { *u = 0 ; for (i = 0 ; i < len ; ++i) { if (*q++ >= sp->qualThresh) *u |= (1 << i) ; }
      ++u ;
    }
  return u0 ;
}

char* qualUnpack (QualPack *sp, U8 *u, char *q, U64 len) /* uncompress (len+7)/8 u into q */
{
  if (!q) q = new(len,char) ;
  char *q0 = q ;
  int i ;
  while (len >= 8)
    { *(U64*)q = sp->qualExpand[*u] ;
      ++u ; q += 8 ; len -= 8 ;
    }
  if (len) for (i = len ; --i ; ) { q[i] = (*u & 0x1) ? sp->qualThresh : 0 ; *u >>= 1 ; }
  return q0 ;
}

/*********** standard conversion tables **************/

int dna2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A',  -2, 'C',  -2,  -2,  -2, 'G',  -2,  -2,  -2,  -2,  -2,  -2, 'N',  -2,
  -2,  -2,  -2,  -2, 'T',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'a',  -2, 'c',  -2,  -2,  -2, 'g',  -2,  -2,  -2,  -2,  -2,  -2, 'n',  -2,
  -2,  -2,  -2,  -2, 't',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, '-',  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'a', 'b', 'c', 'd',  -2,  -2, 'g', 'h',  -2,  -2, 'k',  -2, 'm', 'n',  -2,
  -2,  -2, 'r', 's', 't',  -2, 'v', 'w',  -2, 'y',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbig2NConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'a', 'n', 'c', 'n',  -2,  -2, 'g', 'n',  -2,  -2, 'n',  -2, 'n', 'n',  -2,
  -2,  -2, 'n', 'n', 't',  -2, 'n', 'n',  -2, 'n',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2index4Conv[] = {    /* sends N,n to A - NB also 0,1,2,3 maintained */
   0,   1,   2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int aa2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
} ;

int aa2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
} ;

int noConv[] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
} ;

int acgtCheck[] = {
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
} ;



// complements index and text DNA bases
// ambiguity codes and complements are at: https://droog.gs.washington.edu/parc/images/iupac.html

int complementBase[] = {
   3,   2,   1,   0,   4,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, '-',  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'T', 'V', 'G', 'H',  -2,  -2, 'C', 'D',  -2,  -2, 'M',  -2, 'K', 'N',  -2,
  -2,  -2, 'Y', 'W', 'A',  -2, 'B', 'S',  -2, 'R',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 't', 'v', 'g', 'h',  -2,  -2, 'c', 'd',  -2,  -2, 'm',  -2, 'k', 'n',  -2,
  -2,  -2, 'y', 'w', 'a',  -2, 'b', 's',  -2, 'r',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

/************ BAM access code ****************/

static char *referenceFileName = 0 ;
void seqIOreferenceFileName (char *refFileName) { referenceFileName = strdup (refFileName) ; }

#ifdef BAMIO

#include "sam.h"

typedef struct {
  samFile *f ;
  sam_hdr_t *h ;
  bam1_t *b ;
} BamFile ;

bool bamFileOpenRead (char* filename, SeqIO *si)
{
  BamFile *bf = new0 (1, BamFile) ;

  bf->f = sam_open (filename, "r") ;
  if (!bf->f) return false ;

  uint32_t rf = SAM_FLAG | SAM_SEQ ;
  if (si->isQual) rf |= SAM_QUAL ;
  if (hts_set_opt (bf->f, CRAM_OPT_REQUIRED_FIELDS, rf))
    { fprintf (stderr, "BamFileOpen failed to set CRAM_OPT_REQUIRED_FIELDS\n") ;
      bamFileClose (si) ;
      return false ;
    }
  if (hts_set_opt (bf->f, CRAM_OPT_DECODE_MD, 0))
    { fprintf (stderr, "BamFileOpen failed to set CRAM_OPT_DECODE_MD\n") ;
      bamFileClose (si) ;
      return false ;
    }
  if (referenceFileName && hts_set_fai_filename (bf->f, referenceFileName) != 0)
    { fprintf (stderr, "BamFileOpen failed to set reference genome from %s\n",
	       referenceFileName) ;
      bamFileClose (si) ;
      return false ;
    }

  if (!(bf->h = sam_hdr_read (bf->f)))
    { fprintf(stderr, "BamFileOpen failed to read header\n") ;
      bamFileClose (si) ;
      return false ;
    }

  if (!(bf->b = bam_init1 ()))
    { fprintf(stderr, "BamFileOpen failed to create record buffer\n") ;
      bamFileClose (si) ;
      return false ;
    }

  si->handle = (void*) bf ;
  return true ;
}

static const char binaryAmbigComplement[16] =
  { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 } ;
static const char binaryAmbig2text[] = "=ACMGRSVTWYHKDBN" ;

bool bamRead (SeqIO *si)
{
  int i ;
  BamFile *bf = (BamFile*) si->handle ;

  int res = sam_read1 (bf->f, bf->h, bf->b) ;
  if (res < -1)
    die ("BamFileRead failed to read bam record %d\n", si->nSeq) ;
  if (res == -1) // end of file
    return false ;
  
  si->seqLen = bf->b->core.l_qseq ;
  if (si->seqLen > si->maxSeqLen)
    { if (si->maxSeqLen) { free (si->seqBuf) ; if (si->isQual) free (si->qualBuf) ; }
      si->maxSeqLen = si->seqLen ;
      si->seqBuf = new0 (si->maxSeqLen+1, char) ;
      if (si->isQual) si->qualBuf = new0 (si->maxSeqLen+1, char) ;
    }

  char *bseq = (char*) bam_get_seq (bf->b) ;
  char *s = si->seqBuf ;
  //  if (bf->b->core.flag & BAM_FREVERSE)
  //    for (i = si->seqLen ; i-- ; )
  //      *s++ = binaryAmbig2text[(int)binaryAmbigComplement[(int)bam_seqi(bseq,i)]] ;
  //  else
    for (i = 0 ; i < si->seqLen ; ++i)
      *s++ = binaryAmbig2text[bam_seqi(bseq,i)] ;
  if (si->convert)
    while (s-- > si->seqBuf) 		// NB relies on s == seqBuf+seqLen at start
      *s = si->convert[(int)*s] ;
  
  if (si->isQual)
    { char *bq = (char*) bam_get_qual (bf->b) ;
      if (*bq == '\xff') bzero (si->qualBuf, si->seqLen) ;
      else if (bf->b->core.flag & BAM_FREVERSE)
	{ char *q = si->qualBuf + si->seqLen ;
	  while (q-- > si->qualBuf) *q = *bq ;
	}
      else
	memcpy (si->qualBuf, bq, si->seqLen) ;
    }

  // NB bam_get_qname(bf->b) returns the sequence name

  ++si->nSeq ;
  return true ;
}

void bamFileClose (SeqIO *si)
{ BamFile *bf = (BamFile*) si->handle ;
  if (bf->b) bam_destroy1 (bf->b) ;
  if (bf->h) sam_hdr_destroy (bf->h) ;
  if (bf->f) sam_close (bf->f) ;
  free (bf) ;
}

#endif

/******** end of file ********/
