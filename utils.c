/*  File: utils.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) R Durbin, 1996
 *-------------------------------------------------------------------
 * Description: core utility functions
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 29 18:31 2024 (rd109)
 * * Feb 22 14:52 2019 (rd109): added fzopen()
 * Created: Thu Aug 15 18:32:26 1996 (rd)
 *-------------------------------------------------------------------
 */

#ifdef __linux__
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include "utils.h"

void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  exit (-1) ;
}

void warn (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "WARNING: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
}

static char* commandLine = 0 ;

void storeCommandLine (int argc, char **argv)
{
  int i, totLen = 1 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  if (commandLine) free (commandLine) ;
  commandLine = new (totLen, char) ;
  if (argc) strcpy (commandLine, argv[0]) ; else *commandLine = 0 ;
  for (i = 1 ; i < argc ; ++i) { strcat (commandLine, " ") ; strcat (commandLine, argv[i]) ; }
}

char *getCommandLine (void) { return commandLine ; }

static unsigned long totalAllocated = 0 ;
static unsigned long maxAllocated = 0 ;

void *myalloc (size_t size)
{
  void *p = (void*) malloc (size) ;
  if (!p) die ("myalloc failure requesting %d bytes - totalAllocated %lu", size, totalAllocated) ;
  totalAllocated += size ;
  if (totalAllocated > maxAllocated) maxAllocated = totalAllocated ;
  return p ;
}

void *mycalloc (size_t number, size_t size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p)
    die ("mycalloc failure requesting %ld objects of size %ld - totalAllocated %lu",
	 number, size, totalAllocated) ;
  totalAllocated += size*number ;
  if (totalAllocated > maxAllocated) maxAllocated = totalAllocated ;
  return p ;
}

void  myfree   (void* x, size_t size)
{
  totalAllocated -= size ;
  if (x) free (x) ; // allows to reduce size
}

void *myresize (void* x, size_t nOld, size_t nNew, size_t size)
{
  void *z = mycalloc (nNew, size) ;
  if (nOld < nNew) memcpy (z, x, nOld*size) ; else memcpy (z, x, nNew*size) ;
  myfree (x, nOld*size) ;
  return z ;
}

char *fgetword (FILE *f)
{
  int n = 0 ;
  static char *buf = 0 ;
  int bufSize = 64 ;
  char *cp ;
  if (!buf) buf = (char*) myalloc (bufSize) ;
  cp = buf ;
  while (!feof (f) && (*cp = getc (f)))
    if (isgraph(*cp) && !isspace(*cp))
      { if (++n >= bufSize)
	  { bufSize *= 2 ;
	    if (!(buf = (char*) realloc (buf, bufSize)))
	      die ("fgetword realloc failure requesting %d bytes", bufSize) ;
	    cp = &buf[n] ;
	  }
	else
	  ++cp ;
      }
    else
      { while (*cp != '\n' && !feof(f) && (isspace(*cp) || !isgraph(*cp))) *cp = getc (f) ;
	ungetc (*cp, f) ;
	break ;
      }
  *cp = 0 ;
  return buf ;
}

#define WITH_ZLIB
#ifdef WITH_ZLIB
#include <zlib.h>
#endif

FILE *fzopen(const char *path, const char *mode)
{  /* very cool from https://stackoverflow.com/users/3306211/fernando-mut */
#if defined(WITH_ZLIB)
  gzFile zfp = 0 ;			/* fernando said *zfp - makes me worry.... */

  if (strlen(path) > 3 && !strcmp(&path[strlen(path)-3], ".gz")) // only gzopen on .gz files
    zfp = gzopen(path,mode) ; // try gzopen

  if (!zfp) return fopen(path,mode);

  /* open file pointer */
#ifdef __linux__
  { cookie_io_functions_t io_funcs ;
    io_funcs.read  = (void*) gzread ; // the void* are for cross-compiler happiness
    io_funcs.write = (void*) gzwrite ;
    io_funcs.seek  = (void*) gzseek ;
    io_funcs.close = (void*) gzclose ;
    return fopencookie (zfp, mode, io_funcs) ;
  }
#else // assume a Mac
  return funopen(zfp,
                 (int(*)(void*,char*,int))gzread,
                 (int(*)(void*,const char*,int))gzwrite,
                 (fpos_t(*)(void*,fpos_t,int))gzseek,
                 (int(*)(void*))gzclose) ;
#endif

#else
  return fopen(path,mode);
#endif
}

char* fnameTag (char* root, char* tag)
{
  char *fileName = new (strlen (root) + strlen (tag) + 2, char) ;
  strcpy (fileName, root) ;
  strcat (fileName, ".") ;
  strcat (fileName, tag) ;
  return fileName ;
}

FILE *fopenTag (char* root, char* tag, char* mode)
{
  char *fileName = fnameTag (root, tag) ;
  FILE *f = fzopen (fileName, mode) ;
  free (fileName) ;
  return f ;
}

/***************** rusage for timing information ******************/

#include <sys/resource.h>
#include <sys/time.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

#ifdef RUSAGE_STRUCTURE_DEFINITIONS
struct rusage {
  struct timeval ru_utime; /* user time used */
  struct timeval ru_stime; /* system time used */
  long ru_maxrss;          /* integral max resident set size */
  long ru_ixrss;           /* integral shared text memory size */
  long ru_idrss;           /* integral unshared data size */
  long ru_isrss;           /* integral unshared stack size */
  long ru_minflt;          /* page reclaims */
  long ru_majflt;          /* page faults */
  long ru_nswap;           /* swaps */
  long ru_inblock;         /* block input operations */
  long ru_oublock;         /* block output operations */
  long ru_msgsnd;          /* messages sent */
  long ru_msgrcv;          /* messages received */
  long ru_nsignals;        /* signals received */
  long ru_nvcsw;           /* voluntary context switches */
  long ru_nivcsw;          /* involuntary context switches */
};

struct timeval {
  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
  suseconds_t  tv_usec;  /* and microseconds */
} ;
#endif /* RUSAGE STRUCTURE_DEFINITIONS */

static struct rusage rOld, rFirst ;
static struct timeval tOld, tFirst ;

void timeUpdate (FILE *f)
{
  static bool isFirst = true ;
  struct rusage rNew ;
  struct timeval tNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  gettimeofday(&tNew, 0) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      secs = tNew.tv_sec - tOld.tv_sec ;
      usecs =  tNew.tv_usec - tOld.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\telapsed %d.%06d", secs, usecs) ;
      fprintf (f, "\talloc_max %lu", maxAllocated/1000000) ;
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fputc ('\n', f) ;
    }
  else
    { rFirst = rNew ;
      tFirst = tNew ;
      isFirst = false ;
    }

  rOld = rNew ;
  tOld = tNew ;
}

void timeTotal (FILE *f) { rOld = rFirst ; tOld = tFirst ; timeUpdate (f) ; }

/********************* end of file ***********************/
