#ifndef _ARRAY_H_
#define _ARRAY_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "global.h"

typedef struct _DATA_ {
  void *dptr;
  struct _DATA_ *next;
} DATA;

typedef struct _ARRAY_ {
  short esize;
  short block;
  int   dim;
  DATA  *data;
} ARRAY;

typedef struct _MULTI_ {
  short ndim;
  short esize;
  short *block;
  ARRAY *array;
} MULTI;


#ifdef PERFORM_STATISTICS
typedef struct _ARRAY_TIMING_ {
  clock_t array;
  clock_t multi;
} ARRAY_TIMING;

int GetArrayTiming(ARRAY_TIMING *t);
#endif

int ArrayInit(ARRAY *a, int esize, int block);
void *ArrayGet(ARRAY *a, int i);
void *ArraySet(ARRAY *a, int i, void *d);
void *ArrayAppend(ARRAY *a, void *d);
int ArrayFree(ARRAY *a, void (*FreeElem)(void *));
int ArrayFreeData(DATA *p, int esize, int block, void (*FreeElem)(void *));

int MultiInit(MULTI *ma, int esize, int ndim, int *block);
void *MultiGet(MULTI *ma, int *k);
void *MultiSet(MULTI *ma, int *k, void *d);
int MultiFree(MULTI *ma, void (*FreeElem)(void *));
int MultiFreeData(ARRAY *a, int d, void (*FreeElem)(void *));

#endif
