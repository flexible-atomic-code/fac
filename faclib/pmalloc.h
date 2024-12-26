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

#ifndef _PMALLOC_H_
#define _PMALLOC_H_ 1

#define malloc(x)      pmalloc((x), __FILE__, __LINE__)
#define calloc(n, x)   pcalloc((n), (x), __FILE__, __LINE__)
#define realloc(p, n)  prealloc((p), (n), __FILE__, __LINE__)
#define free(p)        pfree((p), __FILE__, __LINE__)
#define msize()        pmsize()

void *pmalloc(size_t size, char *f, int nline);
void *pcalloc(size_t n, size_t size, char *f, int nline);
void *prealloc(void *p, size_t size, char *f, int nline);
void pfree(void *p, char *f, int nline);
size_t pmsize(void);
void pmalloc_check(void);

#endif
