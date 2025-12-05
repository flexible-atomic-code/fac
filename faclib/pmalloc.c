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

#include <stdio.h>
#include <stdlib.h>

#define MAXNALLOC 16000000
typedef struct _MEM_INFO_ {
  size_t base;
  char *fa, *ff;
  int nla, nlf;
  size_t size;
} MEM_INFO;
  
static int n_alloc = 0;
static size_t _tsize = 0;
static double _dmsize = 0;
static MEM_INFO mem_alloc[MAXNALLOC];
static FILE *pmalloc_log;

size_t pmsize(void) {
  return _tsize;
}

void *pmalloc(size_t size, char *f, int nline) {
  size_t *p;

  p = (size_t *) malloc(size + sizeof(size_t));
  *p = n_alloc;
  mem_alloc[n_alloc].base = (size_t) &p[1];
  mem_alloc[n_alloc].fa = f;
  mem_alloc[n_alloc].nla = nline;
  mem_alloc[n_alloc].ff = NULL;
  mem_alloc[n_alloc].nlf = 0;
  mem_alloc[n_alloc].size = size;
#if PMALLOC == 11
  fprintf(pmalloc_log, "%8s %8d %lx %16s %5d %zu\n", "MALLOC:", n_alloc, (long)(&p[1]), f, nline, size);
  fflush(pmalloc_log);
#endif
#pragma omp atomic
  n_alloc++;
#pragma omp atomic
  _tsize += size;
#pragma omp atomic
  _dmsize += size;
  if (n_alloc == MAXNALLOC) {
    fprintf(pmalloc_log, "MAXNALLOC reached\n");
    exit(1);
  }
  return &p[1];
}

void *pcalloc(size_t n, size_t size, char *f, int nline) {
  size_t *p;
  size_t ns = n*size;
  
  p = (size_t *) calloc(ns+sizeof(size_t), 1);
  *p = n_alloc;
  mem_alloc[n_alloc].base = (size_t) &p[1];
  mem_alloc[n_alloc].fa = f;
  mem_alloc[n_alloc].nla = nline;
  mem_alloc[n_alloc].ff = NULL;
  mem_alloc[n_alloc].nlf = 0;
  mem_alloc[n_alloc].size = ns;
#if PMALLOC == 11
  fprintf(pmalloc_log, "%8s %8d %lx %16s %5d %zu\n", "CALLOC:", n_alloc, (long)(&p[1]), f, nline, ns);
  fflush(pmalloc_log);
#endif
#pragma omp atomic
  n_alloc++;
#pragma omp atomic
  _tsize += ns;
#pragma omp atomic
  _dmsize += ns;
  if (n_alloc == MAXNALLOC) {
    fprintf(pmalloc_log, "MAXNALLOC reached\n");
    exit(1);
  }
  return &p[1];
}

void *prealloc(void *p, size_t size, char *f, int nline) {
  size_t *ps, s0;

  if (p) {
    ps = (size_t *) p;
    ps--;    
    s0 = mem_alloc[*ps].size;
    if (mem_alloc[*ps].ff) {
      fprintf(pmalloc_log, "double free: %d %lx %s(%d), free in %s(%d), alloc in %s(%d)\n",
	      (int)(*ps), (long) &ps[1], f, nline,
	      mem_alloc[*ps].ff, mem_alloc[*ps].nlf,
	      mem_alloc[*ps].fa, mem_alloc[*ps].nla);
      fflush(pmalloc_log);
      mem_alloc[*ps].ff = f;
      mem_alloc[*ps].nlf = -nline;
    } else {
      mem_alloc[*ps].ff = f;
      mem_alloc[*ps].nlf = nline;
    }
#if PMALLOC == 11
    fprintf(pmalloc_log, "%8s %8ld %lx %16s %5d\n", "REALLOC:", *ps, (long)(&ps[1]), f, nline);
    fflush(pmalloc_log);
#endif
#pragma omp atomic
    _tsize -= s0;
  }
  
  ps = (size_t *) realloc(ps, size+sizeof(size_t));
  *ps = n_alloc;
  mem_alloc[n_alloc].base = (size_t) &ps[1];
  mem_alloc[n_alloc].fa = f;
  mem_alloc[n_alloc].nla = nline;
  mem_alloc[n_alloc].ff = NULL;
  mem_alloc[n_alloc].nlf = 0;
  mem_alloc[n_alloc].size = size;
#pragma omp atomic
  n_alloc++;
#pragma omp atomic
  _tsize += size;
#pragma omp atomic
  _dmsize += size;
  if (n_alloc == MAXNALLOC) {
    fprintf(pmalloc_log, "MAXNALLOC reached\n");
    exit(1);
  }

  return &ps[1];
}

void pfree(void *p, char *f, int nline) {
  size_t *ps;
  ssize_t s0;
  
  if (!p) return;
  ps = (size_t *) p;
  ps--;
  
  s0 = mem_alloc[*ps].size;
  if (mem_alloc[*ps].ff) {
    fprintf(pmalloc_log, "double free: %d %lx %s(%d), free in %s(%d), alloc in %s(%d)\n",
	    (int)(*ps), (long) &ps[1], f, nline,
	    mem_alloc[*ps].ff, mem_alloc[*ps].nlf,
	    mem_alloc[*ps].fa, mem_alloc[*ps].nla);
    fflush(pmalloc_log);
    mem_alloc[*ps].ff = f;
    mem_alloc[*ps].nlf = -nline;
  } else {
    mem_alloc[*ps].ff = f;
    mem_alloc[*ps].nlf = nline;
  }

#if PMALLOC == 11
  fprintf(pmalloc_log, "%8s %8d %lx %16s %5d\n", "FREE:", *ps, (long) (&ps[1]), f, nline);
  fflush(pmalloc_log);
#endif

#pragma omp atomic
  _tsize -= s0;
  
  free(ps);
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
  int i;
  int n_leaks;
  size_t n_leakm;
  int n_ilegal;
  size_t tmem;

  n_leaks = 0;
  n_ilegal = 0;
  n_leakm = 0;
  tmem = 0;
  for (i = 0; i < n_alloc; i++) {
    if (mem_alloc[i].ff == NULL) {
      fprintf(pmalloc_log, "%8d: Leak = %d %lx, %5zu, %s (%d)\n",  
	      n_leaks, i, mem_alloc[i].base, mem_alloc[i].size, 
	      mem_alloc[i].fa, mem_alloc[i].nla);
      n_leaks++;
      n_leakm += mem_alloc[i].size;
    } else if (mem_alloc[i].nlf < 0) {
      fprintf(pmalloc_log, "%d: Illegal = %d %lx, %5zu, %s(%d), %s(%d)\n",
	      n_ilegal, i, mem_alloc[i].base, mem_alloc[i].size,
	      mem_alloc[i].fa, mem_alloc[i].nla,
	      mem_alloc[i].ff, mem_alloc[i].nlf);
      n_ilegal++;
    }
    tmem += mem_alloc[i].size;
  }

  fprintf(pmalloc_log, "Total:  %ld %ld %ld\n", tmem, (size_t)_dmsize,  _tsize);
  fprintf(pmalloc_log, "Leaked: %d %ld\n", n_leaks, n_leakm);
  fprintf(pmalloc_log, "Illegal: %d\n", n_ilegal);

  fclose(pmalloc_log);
  return;
}
