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

#ifndef _MPIUTIL_H_
#define _MPIUTIL_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#include "global.h"

#ifndef RBUFL
#define RBUFL 32000000
#endif

#define BUFLN 1024

typedef struct _RANDIDX_ {
  int i;
  double r;
} RANDIDX;

typedef struct _PTRIDX_ {
  int i;
  char *r;
} PTRIDX;

typedef struct _MPID_ {
  int myrank;
  int nproc;
  long long wid;
} MPID;

typedef struct _BFILE_ {
  char *fn;
  FILE *f;
  char *buf;
  int *w, p, n, nbuf;
  int nr, mr, eof;
  LOCK lock;
} BFILE;

int SkipWMPI(int w);
int SkipMPI();
void MPISeqBeg();
void MPISeqEnd();
void MPrintf(int ir, char *format, ...);
int MPIRank(int *np);
int MyRankMPI();
int NProcMPI();
long long WidMPI();
long long CWidMPI();
void SetWidMPI(long long w);
void ResetWidMPI(void);
double WallTime();
MPID *DataMPI();
double TimeSkip();
double TimeLock();
long long NumLock();
void SetLockWT(LOCK *x);
void SetLockMPI(void);
void ReleaseLockMPI(void);
int MPIReady(void);
void Abort(int r);
BFILE *BFileOpen(char *fn, char *md, int nb);
size_t BFileRead(void *ptr, size_t size, size_t nmemb, BFILE *f);
int BFileClose(BFILE *f);
void BFileRewind(BFILE *f);
size_t BFileWrite(void *ptr, size_t size, size_t nmemb, BFILE *bf);
char *BFileGetLine(char *s, int size, BFILE *f);
int BFileSeek(BFILE *bf, long offset, int whence);
long BFileTell(BFILE *bf);
int BFileFlush(BFILE *bf);
void InitializeMPI(int n, int m);
void FinalizeMPI(void);
RANDIDX *RandList(int n);
void RandIntList(int n, int *k);
void ArgSort(int n, double *r, int *k);
int CompareRandIdx(const void *p1, const void *p2);
int ComparePtrIdx(const void *p1, const void *p2);
#endif
