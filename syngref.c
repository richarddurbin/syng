/*  File: syngref.c
 *  Author: Anant Maheshwari (am3320@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: determine locations of syncmers from a <syng>.1khash within a reference (.fa) genome
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 10 16:08 2025 (am3320)
 * Created: Mon Jan 6 14:51  2025 (am3320)
 *-------------------------------------------------------------------
 */

#include "syng.h"
#include "seqhash.h"
#include "seqio.h"
#include "ONElib.h"

/******************* handle syncmer parameters, data structures to store reference positions based on syncmer matches *************/

typedef struct {
  int w, k, seed ;
} Params ;

typedef struct {
  I32	    sync ;	// initially 0 if not in kh; -ve if if reverse orientation
  I32       pos ;       // offset in the current sequence
  bool      forward;    // stores orientation
} SyncPos ;

static int PARAMS_K_DEFAULT = 8 ;
static int PARAMS_W_DEFAULT = 55 ;
static int PARAMS_SEED_DEFAULT = 7 ;

static void readParams (OneFile *of, Params *p)
{ while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType != 'h') die ("sync file %s has no 'h' parameters record", oneFileName(of)) ;
  p->k = oneInt(of,0) ; p->w = oneInt(of,1) ; p->seed = oneInt(of,2) ;
  fprintf (stdout, "read syncmer parameters k %d w %d (size %d) seed %d\n",
	   p->k, p->w, p->k + p->w, p->seed) ;
}

static void checkParams (OneFile *of, Params *p)
{
  if (oneInt(of,0) != p->k ||
      oneInt(of,1) != p->w ||
      oneInt(of,2) != p->seed)
    die ("hash parameters mismatch: (k,w,s) file (%d,%d,%d) != code (%d,%d,%d)",
	 oneInt(of,0), oneInt(of,1), oneInt(of,2), p->k, p->w, p->seed) ;
}

/**************** main program ********************/

static char *usage = 
  "Usage: syngref [options]* <.1khash file> <reference .fa file>\n"
  "possible options are:\n"
  "  -o <outfile prefix>    : default is syngref, file type is .1ref\n";

int main (int argc, char *argv[])
{
    storeCommandLine (argc--, argv++) ;
    timeUpdate (0) ;

    if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }
    char *outPrefix = "syngref" ;
    int nThread = 1 ;
    printf("arguments: %d %s %s \n", argc, argv[0], argv[1]);
    if (argc != 2) die ("missing the two required arguments\n%s", usage) ;

    OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
    if (!schema) die ("problem creating schema - needs debugging/recompilation") ;
    
    // open all the files here
    OneFile *ofK = oneFileOpenRead (argv[0], schema, "khash", 1) ;
    if (!ofK) die ("failed to open .1khash file %s", argv[0]) ;

    SeqIO   *sio = seqIOopenRead (argv[1], dna2indexConv, false) ;
    if (!sio) die ("failed to open sequence file %s", argv[1]) ;
    
    // set up output file here 
    OneFile *ofOut = oneFileOpenWriteNew (fnameTag(outPrefix,"1ref"), schema, "ref", true, nThread) ;
    if (!ofOut) die ("failed to open output file %s", argv[2]) ;
    oneAddProvenance (ofOut, "syngref", SYNG_VERSION, getCommandLine()) ;
    oneAddReference (ofOut, argv[0], 1) ;
    oneAddReference (ofOut, argv[1], 2) ;

    // read the khash  
    Params     params ;
    readParams (ofK, &params) ;                    // read the syncmer hash parameters
    fprintf(stdout, "params: k=%d, w=%d, seed=%d\n", params.k, params.w, params.seed);

    Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly
    KmerHash  *kh = kmerHashReadOneFile (ofK) ;    // read in the kmerhash
    
    if (kh->len != params.w + params.k)
        die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;

    printf("khash length: %d\n", kh->len);
    printf("khash dimension: %d\n", kh->dim);
    printf("khash max: %llu\n", kh->max);
    printf("khash count: %llu\n", kh->count);
    printf ("read %llu kmers from %s\n", kh->max, oneFileName (ofK)) ;
    
    oneFileClose (ofK) ;
    timeUpdate (stdout) ;

    // process the sequences serially
    I64 totSeq = 0;
    Array aSync = arrayCreate(1<<20, SyncPos); 

    while (seqIOread(sio)) {

        I64 seqLen = sio->seqLen;
        I64 sync;

        printf("read sequence %s length %" PRIu64 "\n", sqioId(sio), seqLen);
        char *s = sqioSeq(sio);
        printf("sequence: %s\n", s);

        char *id = sqioId(sio);
        totSeq += seqLen; 
        int pos = 0;
       
        SeqhashIterator *sit = syncmerIterator (sh, s, seqLen) ;
        while (syncmerNext (sit, 0, &pos, 0)) 
        {
            printf("found syncmer at position %d", pos);
            if (kmerHashFind (kh, s+pos, &sync)) {
                printf("found hit: %lld", sync);
                bool forward = true;
                if (sync < 0) {
                    sync = -sync;
                    forward = false;
                    printf("reverse complement\n"); 
                }
                SyncPos *sp = arrayp (aSync, sync, SyncPos) ;
              
                sp->pos = pos;
                sp->forward = forward;
                sp->sync = sync;
            }
             else {
                printf("no hit\n");
            }
        }
        seqhashIteratorDestroy (sit) ;
    }

    printf("SyncPos max: %llu\n", arrayMax(aSync));
    for (U64 i = 0; i < arrayMax(aSync); i++) {
      SyncPos *sp = arrp(aSync, i, SyncPos);
      printf("SyncPos %llu: sync=%d, pos=%d, forward=%s\n", i, sp->sync, sp->pos, sp->forward ? "true" : "false");
    }


    // fill out -- write to .1ref file (using ofOut)
    /*
    "P 3 ref                   REFERENCE INFORMATION\n"
    "O N 2 6 STRING 3 INT      name of sequence, e.g. chr1, and its length\n"
    "O I 1 8 INT_LIST          list of indexes into name table for each sync: 0 if missing, -ve if RC\n"
    "O P 1 8 INT_LIST          positions in sequence\n"

    for (j = 1 ; j < dictMax (nameDict) ; ++i)
    { 
        oneInt(of,1) = arr(chrSize,j,I64) ;
        char *name = dictName(nameDict,j) ;
        oneWriteLine (ofOut, 'N', strlen(name), name) ;
    }

    oneWriteLine (ofOut, 'I', kh->max, index) ;
    oneWriteLine (ofOut, 'P', kh->max, position) ;
    */

    // free allocated memory
    arrayDestroy(aSync); 
    seqIOclose (sio) ;
    oneFileClose (ofOut) ;
    kmerHashDestroy (kh) ;
    oneSchemaDestroy (schema) ;				   
    timeUpdate (stderr) ;
    timeTotal (stderr) ;    
}

