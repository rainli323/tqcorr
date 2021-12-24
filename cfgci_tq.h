/* cfgci.h */
#if !defined (CFGHDR)
#define CFGHDR
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <float.h>
#include <libundmol.h>
#include <libundmolmpi.h>
#include <mpi.h>
#define _FILE_OFFSET_BITS 64
#define KEY          "#CITQ"
#define CINAME       "citq.dat"
#define SCRNAME      "scr.dat"
#define INFONAME     "undmol.dat"
#define MOINTSNAME   "moints.dat"
#define MCFGNAME     "mcfg.dat"
#define CIFTNAME     "cift.dat"
#define LENNAME      128
#define LENLINE      128
#define NTRIALMX     32
#define MAXCFGS      2048          /* arbitrary -- change if needed */
#define NSTMAX       8
#define MAXNBF       384
#define NINTMX       128          /* arbitrary -- change if needed */
#define NEXTMX       ((MAXNBF)-(NINTMX))
#define MAXVRTX      400
#define MXKMAX       19
#define MXSTATE      20
#define MXSPVRTX    ((MXKMAX/2+1)*(MXKMAX/2+2))
#define GRPMUL(i,j) ((i) ^ (j))
#define IOFF(i)     (((i)*((i)-1))/2)
#define MAX(i,j)    ((i)>(j)?(i):(j))
#define MIN(i,j)    ((i)<(j)?(i):(j))
#undef  DEBUG
#endif
