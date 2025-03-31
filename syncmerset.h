/*  File: syncmerset.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 27 19:49 2025 (rd109)
 * Created: Sun Mar 16 17:21:03 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "ONElib.h"
#include "kmerhash.h"

typedef struct {
  int w, k, seed ;
} SyncmerParams ;

SyncmerParams syncmerParamsDefault (void) ;
void          syncmerParamsWrite   (OneFile *of, SyncmerParams p) ;
SyncmerParams syncmerParamsRead    (OneFile *of) ;
void          syncmerParamsCheck   (OneFile *of, SyncmerParams p) ;

typedef struct {
  SyncmerParams params ;
  KmerHash *kh ;
  Array     count ;   // total count across all inputs; of I64 (fast ONEcode read/write)
  I64       totCount ;
  Array     thisCount ; // maximum count in this input; of char - use 1..127
  Array     maxCount ;  // maximum count in any input; of char - use 1..127
} SyncmerSet ;

SyncmerSet *syncmerSetCreate (SyncmerParams params, U64 initialSize) ;
void        syncmerSetDestroy (SyncmerSet *sms) ;
bool        syncmerSetWrite (SyncmerSet *sms, OneFile *of) ;
SyncmerSet *syncmerSetRead (char *filename) ;
void        syncmerUpdateMaxCount (SyncmerSet *sms) ;

static inline void syncmerAdd (SyncmerSet *sms, char *s, I64 *index)
{ bool added = kmerHashAdd (sms->kh, s, index) ;
  if (index)
    { I64 i = (*index < 0) ? -*index : *index ;
      if (added)
	{ array(sms->count, i, I64)++ ;
	  array(sms->thisCount, i, char) = 1 ;
	}
      else
	{ arr(sms->count, i, I64)++ ;
	  if (++arr(sms->thisCount, i, I64) & 0x80) arr(sms->thisCount,i,I64) = 0x7f ;
	}
    }
}

static inline void syncmerCount (SyncmerSet *sms, I32 sync)
{ if (sync < 0) sync = -sync ;
  ++arr(sms->count,sync,I64) ;
  if (++arr(sms->thisCount,sync,I64) & 0x80) arr(sms->thisCount,sync,I64) = 0x7f ;
}

/*************** end of file ******************/
