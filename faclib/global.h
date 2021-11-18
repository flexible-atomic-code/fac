/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _GLOBAL_H_
#define _GLOBAL_H_ 1


/*************************************************************
  Header file defining some global constants, macros.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <pthread.h>

#define PMALLOC_CHECK 2

#if PMALLOC_CHECK == 1
#include "pmalloc.h"
#elif PMALLOC_CHECK == 2
#include "mmalloc.h"
#else
#include "omalloc.h"
#endif

#include "sysdef.h"

#if USE_MPI == 1
#include <mpi.h>
#elif USE_MPI == 2
#include <omp.h>
#endif

/* define constants */
#include "consts.h"

//#define OMP_STAT 1

#define LOCK pthread_mutex_t
#define InitLock(x) pthread_mutex_init((x), NULL)
#define SetLockNT(x) pthread_mutex_lock((x))

#ifdef OMP_STAT
#define SetLock(x) SetLockWT((x))
#else
#define SetLock(x) pthread_mutex_lock((x))
#endif
#define TryLock(x) pthread_mutex_trylock((x))
#define ReleaseLock(x) pthread_mutex_unlock((x))
#define DestroyLock(x) pthread_mutex_destroy((x))

#include "mpiutil.h"

#define USEBF 1
#ifdef USEBF
#define TFILE BFILE
#define FOPEN(fn,m) BFileOpen((fn),(m),(-1))
#define FGETS BFileGetLine
#define FCLOSE BFileClose
#define FREWIND BFileRewind
#define FREAD BFileRead
#define FWRITE BFileWrite
#define FSEEK BFileSeek
#define FTELL BFileTell
#define FFLUSH BFileFlush
#else
#define TFILE FILE
#define FOPEN(fn, m) fopen((fn),(m))
#define FGETS fgets
#define FCLOSE fclose
#define FREWID rewind
#define FREAD fread
#define FWRITE fwrite
#define FSEEK fseek
#define FTELL ftell
#define FFLUSH fflush
#endif
      
/*
** VARIABLE:    DEBUG_RECOUPLE, DEBUG_STRUCTURE, FAC_DEBUG
** TYPE:        macro constants
** PURPOSE:     debugging flags.
** NOTE:        
*/
#define DEBUG_RECOUPLE  10
#define DEBUG_STRUCTURE 20
#define FAC_DEBUG 0
#if FAC_DEBUG
/*
** VARIABLE:    debug_log
** TYPE:        global.
** PURPOSE:     file handler for the output of debug information.
** NOTE:        it is only defined when FAC_DEBUG is true.
*/
extern FILE *debug_log;
#endif

/*
** VARIABLE:    PERFORM_STATISTICS
** TYPE:        macro constant.
** PURPOSE:     indicate whether the profiling 
**              information should be compiled in.
** NOTE:        normally, should be commented out.
*/

/*#define PERFORM_STATISTICS 1*/
#ifdef PERFORM_STATISTICS
extern FILE *perform_log;
#endif

#endif

