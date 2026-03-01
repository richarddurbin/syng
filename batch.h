/*  File: batch.h
 *  Batch helpers for syng: fill, post-process and output thread batches.
 *  Included from syng.c after type definitions (ThreadInfo, SeqInfo, SyncPos, OutType).
 */

// Fill ~100Mb of sequence per thread from sio.  Returns total bytes filled (0 = EOF).
static U64 fillBatch (SeqIO *sio, ThreadInfo *ti, int nThread, U64 *totSeqOut)
{
  U64 filled = 0 ;
  int i ;
  for (i = 0 ; i < nThread ; ++i)
    { arrayMax(ti[i].seq) = 0 ;
      arrayMax(ti[i].seqInfo) = 0 ;
      int seqStart = 0 ;
      while (arrayMax(ti[i].seq) < 100<<20 && seqIOread (sio))
        { arrayp(ti[i].seqInfo, arrayMax(ti[i].seqInfo), SeqInfo)->len = sio->seqLen ;
          array(ti[i].seq, seqStart+sio->seqLen, char) = 0 ;
          memcpy (arrp(ti[i].seq, seqStart, char), sqioSeq(sio), sio->seqLen) ;
          seqStart += sio->seqLen ;
          *totSeqOut += sio->seqLen ;
        }
      filled += arrayMax(ti[i].seq) ;
    }
  seqIOReleaseRead (sio) ;
  return filled ;
}

// After threads join: resize hash, grow count arrays, count syncmers.
static void postProcessBatch (ThreadInfo *ti, int nThread,
    SyncmerSet *sms, SyncmerSet *smsNew, bool isAdd, bool isLast)
{
  int i ;
  I64 j, k ;
  if (isAdd) kmerHashResize (sms->kh) ;
  if (isAdd)
    { I64 curMax = kmerHashMax(sms->kh) ;
      if (curMax >= arrayMax(sms->count))
        { array(sms->count, curMax, I64) = 0 ;
          array(sms->thisCount, curMax, char) = 0 ;
        }
    }
  for (i = 0 ; i < nThread ; ++i)
    { SyncPos *sp = arrp(ti[i].syncPos, 0, SyncPos) ;
      char *seq = arrp(ti[i].seq, 0, char) ;
      for (j = 0 ; j < arrayMax(ti[i].seqInfo) ; ++j)
        { for (k = 0 ; k < arrp(ti[i].seqInfo, j, SeqInfo)->nSync ; ++k, ++sp)
            if (sp->sync) syncmerCount (sms, sp->sync) ;
            else if (smsNew) syncmerAdd (smsNew, seq + sp->pos, 0) ;
          seq += arrp(ti[i].seqInfo, j, SeqInfo)->len ;
        }
    }
  if (isLast) syncmerUpdateMaxCount (sms) ;
}

// Write output for a completed batch: sequences, paths, or GBWT.
// outType is accessed as a file-scope static from syng.c.
static void outputBatch (ThreadInfo *ti, int nThread,
    OneFile *ofOut, SyngBWT *gbwtOut, int khLen, bool isOutputEnds,
    U64 *nSeqP, U64 *nSeq0P, int *nSourceP, U64 *totSyncP)
{
  int i ;
  I64 j, k ;
  for (i = 0 ; i < nThread ; ++i)
    { char *seq = arrp(ti[i].seq, 0, char) ;
      SyncPos *sp = arrp(ti[i].syncPos, 0, SyncPos) ;
      for (j = 0 ; j < arrayMax(ti[i].seqInfo) ; ++j, ++*nSeqP)
        { I64 nSync = arrp(ti[i].seqInfo, j, SeqInfo)->nSync ;
          *totSyncP += nSync ;
          if (arrp(ti[i].seqInfo, j, SeqInfo)->inSource == 1)
            { ++*nSourceP ; *nSeq0P = *nSeqP ; }
          if (outType == SEQ)
            oneWriteLine (ofOut, 'S', arrp(ti[i].seqInfo, j, SeqInfo)->len, seq) ;
          else if (outType == PATH || outType == GBWT)
            { oneInt(ofOut, 0) = arrp(ti[i].seqInfo, j, SeqInfo)->len ;
              oneInt(ofOut, 1) = *nSourceP ;
              oneInt(ofOut, 2) = *nSeqP - *nSeq0P + 1 ;
              oneWriteLine (ofOut, 'P', 0, 0) ;
              if (nSync && outType == GBWT)
                { SyngBWTpath *sbp = syngBWTpathStartNew (gbwtOut, sp->sync) ;
                  oneInt(ofOut, 0) = sp->sync ;
                  oneInt(ofOut, 1) = sp->pos ;
                  oneInt(ofOut, 2) = sbp->jLast ;
                  oneInt(ofOut, 3) = nSync ;
                  oneWriteLine (ofOut, 'Z', 0, 0) ;
                  for (k = 1 ; k < nSync ; ++k)
                    syngBWTpathAdd (sbp, sp[k].sync, sp[k].pos - sp[k-1].pos) ;
                  syngBWTpathFinish (sbp) ;
                  sbp = syngBWTpathStartNew (gbwtOut, -sp[nSync-1].sync) ;
                  for (k = nSync-2 ; k >= 0 ; --k)
                    syngBWTpathAdd (sbp, -sp[k].sync, sp[k+1].pos - sp[k].pos) ;
                  syngBWTpathFinish (sbp) ;
                }
              else if (nSync && outType == PATH)
                { static I64 *x = 0 ;
                  static size_t xSize = 0 ;
                  if (!x)
                    { xSize = nSync ; x = new (xSize, I64) ; }
                  else if (xSize < nSync)
                    { newFree (x,xSize,I64) ; xSize = nSync ; x = new (xSize, I64) ; }
                  for (k = 0 ; k < nSync ; ++k) x[k] = sp[k].sync ;
                  oneWriteLine (ofOut, 'z', nSync, x) ;
                  for (k = 0 ; k < nSync ; ++k) x[k] = sp[k].pos ;
                  oneWriteLine (ofOut, 'o', nSync, x) ;
                }
              if (isOutputEnds)
                { I64 len = arrp(ti[i].seqInfo, j, SeqInfo)->len ;
                  if (nSync)
                    { oneWriteLine (ofOut, 'X', sp->pos, seq) ;
                      I64 endOff = sp[nSync-1].pos + khLen ;
                      oneWriteLine (ofOut, 'Y', len - endOff, seq + endOff) ;
                    }
                  else
                    { oneWriteLine (ofOut, 'X', len/2, seq) ;
                      oneWriteLine (ofOut, 'Y', len - len/2, seq + len/2) ;
                    }
                }
              sp += nSync ;
            } // PATH or GBWT
          seq += arrp(ti[i].seqInfo, j, SeqInfo)->len ;
        } // j sequences
    } // i threads
}

/******************** end of file ********************/
