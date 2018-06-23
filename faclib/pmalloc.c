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

#define _PMALLOC_H_ 1

#include <stdio.h>
#include <stdlib.h>

#define MAXNALLOC 1000000
typedef struct _MEM_INFO_ {
  long int base;
  char *f;
  int nline;
  size_t size;
} MEM_INFO;
  
static int n_alloc = 0;
static int n_free = 0;
static size_t _tsize = 0;
static MEM_INFO mem_alloc[MAXNALLOC];
static MEM_INFO mem_free[MAXNALLOC];
static FILE *pmalloc_log;

size_t pmsize(void) {
  return _tsize;
}

void *pmalloc(size_t size, char *f, int nline) {
  void *p;

  p = malloc(size);
  mem_alloc[n_alloc].base = (long int) p;
  mem_alloc[n_alloc].f = f;
  mem_alloc[n_alloc].nline = nline;
  mem_alloc[n_alloc].size = size;
  fprintf(pmalloc_log, "%8s %lx %16s %5d %zu\n", "MALLOC:", (long)p, f, nline, size);
  n_alloc++;
  if (n_alloc == MAXNALLOC) {
    fprintf(pmalloc_log, "MAXNALLOC reached\n");
    exit(1);
  }
  return p;
}

void *pcalloc(size_t n, size_t size, char *f, int nline) {
  void *p;

  p = calloc(n, size);
  mem_alloc[n_alloc].base = (long int) p;
  mem_alloc[n_alloc].f = f;
  mem_alloc[n_alloc].nline = nline;
  mem_alloc[n_alloc].size = size*n;
  fprintf(pmalloc_log, "%8s %lx %16s %5d %zu\n", "CALLOC:", (long)p, f, nline, size*n);
  
  n_alloc++;
  if (n_alloc == MAXNALLOC) {
    fprintf(pmalloc_log, "MAXNALLOC reached\n");
    exit(1);
  }
  return p;
}

void *prealloc(void *p, size_t size, char *f, int nline) {
  void *q;

  mem_free[n_free].base = (long int) p;
  mem_free[n_free].f = f;
  mem_free[n_free].nline = nline;
  fprintf(pmalloc_log, "%8s %lx %16s %5d\n", "REALLOC:", (long)p, f, nline);
  n_free++;

  q = realloc(p, size);
  mem_alloc[n_alloc].base = (long int) q;
  mem_alloc[n_alloc].f = f;
  mem_alloc[n_alloc].nline = nline;
  mem_alloc[n_alloc].size = size;
  n_alloc++;
  if (n_alloc == MAXNALLOC) {
    fprintf(pmalloc_log, "MAXNALLOC reached\n");
    exit(1);
  }

  return q;
}

void pfree(void *p, char *f, int nline) {  
  mem_free[n_free].base = (long int) p;
  mem_free[n_free].f = f;
  mem_free[n_free].nline = nline;
  fprintf(pmalloc_log, "%8s %lx %16s %5d\n", "FREE:", (long)p, f, nline);
  fflush(pmalloc_log);
  n_free++;
  
  free(p);
}

int CompareMemory(const void *p1, const void *p2) {
  long int b1, b2;

  b1 = ((MEM_INFO *) p1)->base;
  b2 = ((MEM_INFO *) p2)->base;

  if (b1 < b2) return -1;
  else if (b1 > b2) return 1;
  else return 0;
}

void pmalloc_open(void) {
  pmalloc_log = fopen("pmalloc.log", "w");
}

void pmalloc_check(void) {
  int i, j;
  int n_leaks;
  int n_leakm;
  int n_ilegal;
  int tmem;

  if (n_alloc == n_free) {
    fprintf(pmalloc_log, "no mem leak\n");
    return;
  }

  qsort(mem_alloc, n_alloc, sizeof(MEM_INFO), CompareMemory);
  qsort(mem_free, n_free, sizeof(MEM_INFO), CompareMemory);

  i = 0; 
  j = 0;
  n_leaks = 0;
  n_leakm = 0;
  n_ilegal = 0;
  tmem = 0;
  while (i < n_alloc && j < n_free) {
    if (mem_alloc[i].base < mem_free[j].base) {
      fprintf(pmalloc_log, "%6d: Leak = %lx, %5zu, %30s (%d)\n",  
	     n_leaks, mem_alloc[i].base, mem_alloc[i].size, 
	     mem_alloc[i].f, mem_alloc[i].nline);
      n_leaks++;
      n_leakm += mem_alloc[i].size;
      i++;
    } else if (mem_alloc[i].base == mem_free[j].base) {
      tmem += mem_alloc[i].size;
      i++;
      j++;
    } else {
      fprintf(pmalloc_log, "%6d: Illegal Free = %lx, in File %s (%d)\n", 
	     n_ilegal, mem_free[j].base, mem_free[j].f, mem_free[j].nline);
      fprintf(pmalloc_log, "%d %d %d %d\n", i, n_alloc, j, n_free);
      n_ilegal++;
      j++;
    }
  }

  tmem += n_leakm;

  fprintf(pmalloc_log, "Total:  %d\n", tmem);
  fprintf(pmalloc_log, "Leaked: %d\n", n_leakm);

  return;
}
