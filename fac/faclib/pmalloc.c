#define _PMALLOC_H_ 1

#include <stdio.h>
#include <stdlib.h>

#define MAXNALLOC 1000000
typedef struct _MEM_INFO_ {
  off_t base;
  char *f;
  int nline;
  size_t size;
} MEM_INFO;
  
static int n_alloc = 0;
static int n_free = 0;
static MEM_INFO mem_alloc[MAXNALLOC];
static MEM_INFO mem_free[MAXNALLOC];

void *pmalloc(size_t size, char *f, int nline) {
  void *p;

  p = malloc(size);
  mem_alloc[n_alloc].base = (off_t) p;
  mem_alloc[n_alloc].f = f;
  mem_alloc[n_alloc].nline = nline;
  mem_alloc[n_alloc].size = size;
  n_alloc++;
  if (n_alloc == MAXNALLOC) {
    printf("MAXNALLOC reached\n");
    exit(1);
  }
  return p;
}

void *pcalloc(size_t n, size_t size, char *f, int nline) {
  void *p;

  p = calloc(n, size);
  mem_alloc[n_alloc].base = (off_t) p;
  mem_alloc[n_alloc].f = f;
  mem_alloc[n_alloc].nline = nline;
  mem_alloc[n_alloc].size = size*n;
  
  n_alloc++;
  if (n_alloc == MAXNALLOC) {
    printf("MAXNALLOC reached\n");
    exit(1);
  }
  return p;
}

void *prealloc(void *p, size_t size, char *f, int nline) {
  void *q;

  mem_free[n_free].base = (off_t) p;
  mem_free[n_free].f = f;
  mem_free[n_free].nline = nline;
  n_free++;

  q = realloc(p, size);
  mem_alloc[n_alloc].base = (off_t) q;
  mem_alloc[n_alloc].f = f;
  mem_alloc[n_alloc].nline = nline;
  mem_alloc[n_alloc].size = size;
  n_alloc++;
  if (n_alloc == MAXNALLOC) {
    printf("MAXNALLOC reached\n");
    exit(1);
  }

  return q;
}

void pfree(void *p, char *f, int nline) {  
  mem_free[n_free].base = (off_t) p;
  mem_free[n_free].f = f;
  mem_free[n_free].nline = nline;
  n_free++;
  
  free(p);
}

int _CompareMemory(const void *p1, const void *p2) {
  off_t b1, b2;

  b1 = ((MEM_INFO *) p1)->base;
  b2 = ((MEM_INFO *) p2)->base;

  if (b1 < b2) return -1;
  else if (b1 > b2) return 1;
  else return 0;
}

void pmalloc_check(void) {
  int i, j;
  int n_leaks;
  int n_leakm;
  int n_ilegal;
  int tmem;

  if (n_alloc == n_free) {
    printf("no mem leak\n");
    return;
  }

  qsort(mem_alloc, n_alloc, sizeof(MEM_INFO), _CompareMemory);
  qsort(mem_free, n_free, sizeof(MEM_INFO), _CompareMemory);

  i = 0; 
  j = 0;
  n_leaks = 0;
  n_leakm = 0;
  n_ilegal = 0;
  tmem = 0;
  while (i < n_alloc && j < n_free) {
    if (mem_alloc[i].base < mem_free[j].base) {
      printf("%6d: Leak = %lX, %5d, %30s (%d)\n",  
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
      printf("%6d: Illegal Free = %lx, in File %s (%d)\n", 
	     n_ilegal, mem_free[j].base, mem_free[j].f, mem_free[j].nline);
      printf("%d %d %d %d\n", i, n_alloc, j, n_free);
      n_ilegal++;
      j++;
    }
  }

  tmem += n_leakm;

  printf("Total:  %d\n", tmem);
  printf("Leaked: %d\n", n_leakm);

  return;
}
