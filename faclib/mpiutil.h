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
#define RBUFL 1280000
#endif

#define BUFLN 1024

typedef struct _BFILE_ {
  char *fn;
  FILE *f;
  void *buf;
  int p, n, nbuf;
  int nr, mr, eof;
} BFILE;

void MPISeqBeg();
void MPISeqEnd();
void MPrintf(int ir, char *format, ...);
int MPIRank(int *np);
int MPIReady();
BFILE *BFileOpen(char *fn, char *md, int nb);
size_t BFileRead(void *ptr, size_t size, size_t nmemb, BFILE *f);
int BFileClose(BFILE *f);
void BFileRewind(BFILE *f);
char *BFileGetLine(char *s, int size, BFILE *f);

#endif
