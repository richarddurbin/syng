/*  File: syncmerset.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 27 20:24 2025 (rd109)
 * Created: Sun Mar 16 17:20:12 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "syncmerset.h"

/******************* package to handle syncmer parameters *************/

SyncmerParams syncmerParamsDefault (void)
{ SyncmerParams p ;
  p.k = 8 ;
  p.w = 55 ;
  p.seed = 7 ;
  return p ;
}

void syncmerParamsWrite (OneFile *of, SyncmerParams p)
{
  oneInt(of,0) = p.k ;
  oneInt(of,1) = p.w ;
  oneInt(of,2) = p.seed ;
  oneWriteLine (of, 'h', 0, 0) ;
}

SyncmerParams syncmerParamsRead (OneFile *of)
{ SyncmerParams p ;
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType != 'h') die ("sync file %s has no 'h' parameters record", oneFileName(of)) ;
  p.k = oneInt(of,0) ; p.w = oneInt(of,1) ; p.seed = oneInt(of,2) ;
  fprintf (stdout, "read syncmer parameters k %d w %d (size %d) seed %d\n",
	   p.k, p.w, p.k + p.w, p.seed) ;
  return p ;
}

void syncmerParamsCheck (OneFile *of, SyncmerParams p)
{
  if (of->lineType != 'h') die ("syncmerParamsCheck called on a %c line not 'h'", of->lineType) ;
  if (oneInt(of,0) != p.k ||
      oneInt(of,1) != p.w ||
      oneInt(of,2) != p.seed)
    die ("hash parameters mismatch: (k,w,s) file (%d,%d,%d) != code (%d,%d,%d)",
	 oneInt(of,0), oneInt(of,1), oneInt(of,2), p.k, p.w, p.seed) ;
}

/*************** syncmerHash package built on kmerhash *****************/

SyncmerSet *syncmerSetCreate (SyncmerParams params, U64 initialSize)
{
  SyncmerSet *sms = new0 (1, SyncmerSet) ;
  sms->params = params ;
  sms->kh = kmerHashCreate (initialSize, params.w + params.k) ;
  I64 aSize = (1 << (sms->kh->dim-2)) + 1 ;
  sms->count = arrayCreate (aSize, I64) ;
  sms->thisCount = arrayCreate (aSize, char) ;
  sms->maxCount = arrayCreate (aSize, char) ;
  return sms ;
}

void syncmerSetDestroy (SyncmerSet *sms)
{
  kmerHashDestroy (sms->kh) ;
  arrayDestroy (sms->count) ;
  arrayDestroy (sms->thisCount) ;
  arrayDestroy (sms->maxCount) ;
  newFree (sms, 1, SyncmerSet) ;
}

void syncmerUpdateMaxCount (SyncmerSet *sms)
{
  I64 i ;
  if (arrayMax(sms->thisCount) > arrayMax(sms->maxCount))
    array(sms->maxCount, arrayMax(sms->thisCount)-1, char) = 0 ;
  for (i = 1 ; i < arrayMax(sms->thisCount) ; ++i)
    if (arr(sms->thisCount,i,char) > arr(sms->maxCount,i,char))
      arr(sms->maxCount,i,char) = arr(sms->thisCount,i,char) ;
  memset (arrayp(sms->thisCount,0,char), 0, arrayMax(sms->thisCount)) ;
}

// syncmerSetWrite takes *of not filename so can record sources in upstream file

bool syncmerSetWrite (SyncmerSet *sms, OneFile *of) 
{
  syncmerParamsWrite (of, sms->params) ;

  KmerHash *kh = sms->kh ;
  if (!kmerHashWriteOneFile (kh, of)) { oneFileClose (of) ; return false ; }

  I64 chunk = (I64)1 << 22 ;              // chunk size in 8-byte words, giving 128 Mb chunks
  I64 total = kh->max ;
  
  I64 *counts = arrp (sms->count, 0, I64) ;        // write the counts
  while (total > chunk) { oneWriteLine (of, 'C', chunk, counts) ; counts += chunk ; total -= chunk ; }
  if (total) oneWriteLine (of, 'C', total, counts) ;
  
  char *mc = arrp (sms->maxCount, 0, char) ;        // write maxCount
  total = kh->max ;
  while (total > chunk) { oneWriteLine (of, 'M', chunk, mc) ; mc += chunk ; total -= chunk ; }
  if (total) oneWriteLine (of, 'M', total, mc) ;

  fprintf (stdout, "wrote %llu syncmers to file %s\n", kmerHashMax(sms->kh), oneFileName(of)) ;
  timeUpdate (stdout) ;
  return true ;
}

static char *schemaText =
  "1 3 def 1 0               schema for SyncmerSet\n"
  ".\n"
  "P 5 khash                 KMER HASH\n"
  "S 7 syncset               SYNCMER SET\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n"
  "O t 3 3 INT 3 INT 3 INT   max, len, dim for KmerHash table\n"
  "D S 1 3 DNA               packed sequences aligned to 64-bit boundaries\n" 
  "D L 1 8 INT_LIST          locations in the table\n"
  "D C 1 8 INT_LIST          kmer counts\n"
  "D M 1 6 STRING            maximum count in any input - (1..127)\n"
  ;

SyncmerSet *syncmerSetRead (char *filename)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (filename, schema, "syncset", 1) ;
  oneSchemaDestroy (schema) ;
  if (!of) die ("failed to open syncmer file %s to read", filename) ;
  
  SyncmerParams params = syncmerParamsRead (of) ;
  
  KmerHash *kh = kmerHashReadOneFile (of) ; // dies if fails
  SyncmerSet *sms = syncmerSetCreate (params, kh->max) ;
  kmerHashDestroy (sms->kh) ; sms->kh = kh ;
  if (sms->kh->len != params.w + params.k)
    die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;

  // read the counts
  array (sms->count, kh->max, I64) = 0 ; // ensure array is big enough
  I64 n = 0 ;
  while (of->lineType == 'C')
    { memcpy (arrp(sms->count, n, I64), oneIntList(of), oneLen(of)*sizeof(I64)) ;
      n += oneLen(of) ;
      oneReadLine(of) ;
    }
  if (n != kh->max) die ("wrong number of C counts %llx read for syncmerset %s size %llx",
			 n, filename, kh->max) ;

  // read maxCount
  array (sms->maxCount, kh->max, char) = 0 ; // ensure array is big enough
  n = 0 ;
  while (of->lineType == 'M')
    { memcpy (arrp(sms->maxCount, n, char), oneString(of), oneLen(of)) ;
      n += oneLen(of) ;
      oneReadLine(of) ;
    }
  if (n != kh->max) die ("wrong number of M counts %llx read for syncmerset %s size %llx",
			 n, filename, kh->max) ;

  oneFileClose (of) ; // close the old file

  U64 i, totCount = 0 ;
  for (i = 1 ; i <= kmerHashMax(kh) ; ++i) totCount += arr(sms->count,i,I64) ;
  fprintf (stdout, "read %llu syncmers from %s with total count %llu\n",
	   kmerHashMax(kh), filename, totCount) ;

  timeUpdate (stdout) ;

  return sms ;
}

/******************* end of file *****************/
