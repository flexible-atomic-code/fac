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

#ifndef _MMALLOC_H_
#define _MMALLOC_H_ 1

#define malloc(x)      mmalloc((x))
#define calloc(n, x)   mcalloc((n),(x))
#define realloc(p, n)  mrealloc((p),(n))
#define free(p)        mfree((p))
#define msize()        mmsize()

void *mmalloc(size_t size);
void *mcalloc(size_t n, size_t x);
void *mrealloc(void *p, size_t n);
void mfree(void *p);
size_t mmsize(void);

#endif
