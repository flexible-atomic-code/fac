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
#include "omalloc.h"

size_t omsize(void) {
#if PMALLOC == 2
  size_t ts;
  int r;
  
  MallocExtension_ReleaseFreeMemory();
  r = MallocExtension_GetNumericProperty("generic.current_allocated_bytes",&ts);
  if (r) {
    return ts;
  } else {
    printf("tcmalloc allocated error: %d\n", r);
    return 0;
  }
#elif PMALLOC == 3
  size_t ts, ss, ep;
  int r;
  
  ss = sizeof(size_t);
  r = mallctl("epoch", NULL, NULL, &ep, ss);
  if (r != 0) {
    printf("jemalloc epoch error: %d\n", r);
  }
  r = mallctl("stats.allocated", &ts, &ss, NULL, 0);
  if (r == 0) {
    return ts;
  } else {
    printf("jemalloc allocated error: %d\n", r);
  }
#else
  return 0;
#endif
}
