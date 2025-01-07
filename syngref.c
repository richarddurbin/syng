/*  File: syngref.c
 *  Author: Anant Maheshwari (am3320@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: determine locations of syncmers in a <syng>.1khash within a reference (.fa) genome
 * Exported functions:
 * HISTORY:
 * Last edited: Jan  7 11:16 2025 (am3320)
 * Created: Mon Jan  6 14:51  2024 (am3320)
 *-------------------------------------------------------------------
 */

#include "syng.h"
#include "seqhash.h"
#include "seqio.h"
#include "ONElib.h"


/******************* package to handle syncmer parameters *************/

typedef struct {
  int w, k, seed ;
} Params ;

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
  "  -T <threads>           : [8] number of threads\n"
  "  -o <outfile prefix>    : default is syngref, file type is .1ref\n";

int main (int argc, char *argv[])
{
    fprintf(stdout, "Hello, world!\n");

    storeCommandLine (argc--, argv++) ;
    timeUpdate (0) ;

    if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }
      char *outPrefix = "syngref" ;

    if (argc != 2) die ("missing the three required arguments\n%s", usage) ;

    OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
    if (!schema) die ("problem creating schema - needs debugging/recompilation") ;
    
      // open all the files here, so we can die quickly if any file opens fail
    OneFile *ofK = oneFileOpenRead (argv[0], schema, "khash", 1) ;
    if (!ofK) die ("failed to open .1khash file %s", argv[0]) ;

    SeqIO   *sio = seqIOopenRead (argv[1], dna2indexConv, false) ;
    if (!sio) die ("failed to open sequence file %s", argv[1]) ;
    
    // set up output file here 
    // TODO - add output file

    // read the khash  
    Params     params ;
    readParams (ofK, &params) ;                    // read the syncmer hash parameters
    fprintf(stdout, "params: k=%d, w=%d, seed=%d\n", params.k, params.w, params.seed);

    Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here, awkwardly
    KmerHash  *kh = kmerHashReadOneFile (ofK) ;    // read in the kmerhash
    
    if (kh->len != params.w + params.k)
        die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;
    printf ("read %llu kmers from %s\n", kh->max, oneFileName (ofK)) ;
    oneFileClose (ofK) ;
    timeUpdate (stdout) ;

    int i, j ; // general index variables

    seqIOclose (sio) ;
    kmerHashDestroy (kh) ;
    oneSchemaDestroy (schema) ;				   
                    
    timeUpdate (stderr) ;
    timeTotal (stderr) ;
    
    // clean up memory
}

