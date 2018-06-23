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

#define _MMALLOC_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t _tsize = 0;

void *mmalloc(size_t size) {
  size_t *p = NULL;

  p = (size_t *) malloc(size+sizeof(size_t));
  if (p == NULL) {
    printf("malloc error: %zu %zu\n", size, _tsize);
    exit(1);
  }
  *p = size;
#pragma omp atomic
  _tsize += size;
  return &p[1];
}

void *mcalloc(size_t n, size_t size) {
  size_t *p;
  size_t ns = n*size;

  p = (size_t *) calloc(ns+sizeof(size_t), 1);
  if (p == NULL) {
    printf("calloc error: %zu %zu %zu\n", n, size, _tsize);
    exit(1);
  }
  *p = ns;
#pragma omp atomic
  _tsize += ns;
  return &p[1];
}

void *mrealloc(void *p, size_t size) {
  size_t *ps;

  if (p) {
    ps = (size_t *) p;
    ps--;
#pragma omp atomic
    _tsize -= ps[0];
  }
  ps = (size_t *) realloc(ps, size+sizeof(size_t));
  if (ps == NULL) {
    printf("realloc error: %zu %zu\n", size, _tsize);
    exit(1);
  }
  *ps = size;
#pragma omp atomic
  _tsize += size;
  return &ps[1];
}

void mfree(void *p) {
  size_t *ps;

  if (!p) return;
  ps = (size_t *) p;
  ps--;
#pragma omp atomic
  _tsize -= ps[0];
  free(ps);
}

size_t mmsize(void) {
  return _tsize;
}
