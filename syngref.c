/*  File: syngref.c
 *  Author: Anant Maheshwari (am3320@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: determine locations of syncmers from a <syng>.1khash within a reference (.fa) genome
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 12 01:08 2025 (am3320)
 * Created: Mon Jan 6 14:51  2025 (am3320)
 *-------------------------------------------------------------------
 */

#include "syng.h"
#include "dict.h"
#include "seqhash.h"
#include "seqio.h"
#include "ONElib.h"

/******************* handle syncmer parameters, data structures to store reference positions based on syncmer matches *************/

typedef struct {
  int w, k, seed ;
} Params ;

typedef struct {
  I32       pos ;       // position in reference
  I32      occurences;     // # of hits from reference
  int      chrom;          // chromosome number
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
    printf("arguments: %d %s %s \n", argc, argv[0], argv[1]);

    char *outPrefix = "syngref" ;
    if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }

    int nThread = 1 ;
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
    I64 nSync = kmerHashMax(kh);
    printf("nSync: %llu\n", nSync);

    if (kh->len != params.w + params.k)
        die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;
    
    oneFileClose (ofK) ;

    // process the sequences serially
    I64 totSeq = 0;
    int currentChrom = 1;
    Array aSync = arrayCreate(nSync, SyncPos); 
    DICT *nameDict = dictCreate(1<<12); // 4096 chromosomes should be a decent bound estimate

    while (seqIOread(sio)) {

        I64 seqLen = sio->seqLen;
        I64 sync;

        char *s = sqioSeq(sio);
        char *id = sqioId(sio);

        U64 index;
        dictAdd (nameDict, id, &index) ;
        totSeq += seqLen; 
        int pos = 0;
       
        SeqhashIterator *sit = syncmerIterator (sh, s, seqLen) ;
        while (syncmerNext (sit, 0, &pos, 0)) 
        {
            if (kmerHashFind (kh, s+pos, &sync)) {
                bool forward = true;
                if (sync < 0) {
                    sync = -sync;
                    forward = false;
                }
                SyncPos *sp = arrayp (aSync, sync, SyncPos) ;
                sp->pos = pos;
                sp->chrom = forward ? currentChrom : -currentChrom;
                sp->occurences += 1;
            }
        }
        seqhashIteratorDestroy (sit) ;
        currentChrom += 1;
    }

    for (int j = 0 ; j < dictMax (nameDict) ; ++j)
    { 
        char *name = dictName(nameDict,j) ;
        oneWriteLine (ofOut, 'N', strlen(name), name) ;
    }
    
    I64 *x = new (nSync, I64);
    I64 *y = new (nSync, I64);
    I64 *z = new (nSync, I64);

    for (U64 i = 0; i < nSync; i++) {
      SyncPos *sp = arrp(aSync, i, SyncPos);
      if (sp->occurences < 1) {
        x[i] = 0; // indicate missing
        y[i] = -1; // indicate missing (negative since 0 is a valid sequence position)
      } else {
        x[i] = sp->chrom;
        y[i] = sp->pos;
      }
      z[i] = sp->occurences;
    }

    oneWriteLine (ofOut, 'I', nSync, x);
    oneWriteLine (ofOut, 'P', nSync, y);
    oneWriteLine (ofOut, 'O', nSync, z);

    newFree (x, nSync, I64);
    newFree (y, nSync, I64);
    newFree (z, nSync, I64);
    
    // free allocated memory
    arrayDestroy(aSync); 
    dictDestroy(nameDict);
    seqIOclose (sio) ;
    oneFileClose (ofOut) ;
    kmerHashDestroy (kh) ;
    oneSchemaDestroy (schema) ;				   
    timeTotal (stderr) ;    
}

