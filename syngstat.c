/*  File: syngstat.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: statistics on syng files
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 15 16:24 2026 (rd109)
 * Created: Sun Mar 16 17:11:08 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include "seqhash.h"
#include "syng.h"

/**************** main program ********************/

static char usage[] =
  "Usage: syngstat <syng 1file>*\n" ;

int main (int argc, char *argv[])
{
  timeUpdate (0) ;

  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;

  while (argc--)
    { OneFile *of = oneFileOpenRead (*argv++, 0, 0, 1) ;
      if (!of) die ("failed to open ONEfile %s", *--argv) ;

      if (of->subType && !strcmp (of->subType, "gbwt"))
	{ I64 nZ = 0, nV = 0, nE = 0, ne = 0 ;
	  oneStats (of, 'V', &nV, 0, 0) ;
	  oneStats (of, 'E', &nE, 0, 0) ;
	  oneStats (of, 'e', &ne, 0, 0) ;
	  oneStats (of, 'Z', &nZ, 0, 0) ;
	  printf ("GBWT file %s with %'lld vertices, %'lld edges and %'lld sequences\n",
		  oneFileName(of), nV, (nE+ne)/2, nZ) ;
	  if (nV)
	    { SyngBWT *sb = syngBWTread (of) ;
	      syngBWTstat (sb) ;
	      syngBWTdestroy (sb) ;
	    }
	  timeUpdate (stdout) ;
	}

      oneFileClose (of) ;
    }
}
