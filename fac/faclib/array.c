#include "array.h"

static char *rcsid="$Id: array.c,v 1.6 2001/10/19 22:45:38 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#ifdef PERFORM_STATISTICS  
static ARRAY_TIMING timing = {0, 0};

int GetArrayTiming(ARRAY_TIMING *t) {
  memcpy(t, &timing, sizeof(timing)); 
  return 0;
}
#endif

/******************************************************************/
/* implements a variable length one- and multi- dimensional array */
/******************************************************************/

int ArrayInit(ARRAY *a, int esize, int block) {
  a->esize = esize;
  a->block = block;
  a->dim = 0;
  a->data = NULL;
  return 0;
}

void *ArrayGet(ARRAY *a, int i) {
  DATA *p;
  if (i < 0 || i >= a->dim) return NULL;
  p = a->data;
  while (i >= a->block) {
    p = p->next;
    i -= a->block;
  }
  if (p->dptr) {
    return ((char *) p->dptr) + i*(a->esize);
  } else {
    return NULL;
  }
}

void *ArraySet(ARRAY *a, int i, void *d) {
  void *pt;
  DATA *p;
  /*
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif
  */

  if (a->dim == 0) {
    a->data = (DATA *) malloc(sizeof(DATA));
    a->data->dptr = (void *) calloc(a->block, a->esize);
    a->data->next = NULL;
  }
  p = a->data;
  if (a->dim <= i) a->dim = i+1;
  while (i >= a->block) {
    if (!(p->next)) {
      p->next = (DATA *) malloc(sizeof(DATA));
      p->next->dptr = NULL;
      p->next->next = NULL;
    }
    p = p->next;
    i -= a->block;
  }

  if (!(p->dptr)) {
    p->dptr = (void *) calloc(a->block, a->esize);
  }
  pt = ((char *) p->dptr) + i*a->esize;
  
  if (d) memcpy(pt, d, a->esize);

  /*
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.array += stop-start;
#endif 
  */
  return pt;
}

void *ArrayAppend(ARRAY *a, void *d) {
  int i;  
  i = a->dim;
  return ArraySet(a, i, d);
}

int ArrayFreeData(DATA *p, int esize, int block, 
		  void (*FreeElem)(void *)) {
  void *pt;
  int i;

  if (p->next) {
    ArrayFreeData(p->next, esize, block, FreeElem);
  }
  
  if (FreeElem && p->dptr) {
    pt = p->dptr;
    for (i = 0; i < block; i++) {
      FreeElem(pt);
      pt = ((char *) pt) + esize;
    }
  }
  if (p->dptr) free(p->dptr);
  if (p) free(p);
  p = NULL;
  return 0;
}
    
int ArrayFree(ARRAY *a, void (*FreeElem)(void *)) {
  if (!a) return 0;
  if (a->dim == 0) return 0;
  ArrayFreeData(a->data, a->esize, a->block, FreeElem);
  a->dim = 0;
  a->data = NULL;
  return 0;
}

int MultiInit(MULTI *ma, int esize, int ndim, int *block) {
  int i;
  ma->ndim = ndim;
  ma->esize = esize;
  ma->block = (short *) malloc(sizeof(short)*ndim);
  for (i = 0; i < ndim; i++) ma->block[i] = block[i];
  ma->array = NULL;
  return 0;
}

void *MultiGet(MULTI *ma, int *k) {
  ARRAY *a;
  int i;
  a = ma->array;
  if (a == NULL) return NULL;
  for (i = 0; i < ma->ndim; i++) {
    a = (ARRAY *) ArrayGet(a, k[i]);
    if (a == NULL) return NULL;
  }
  
  return (void *) a;
}

void *MultiSet(MULTI *ma, int *k, void *d) {
  ARRAY *a;
  void *pt;
  int i, ndim1, ndim2;
  /*
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif
  */

  if (ma->array == NULL) {
    ma->array = (ARRAY *) malloc(sizeof(ARRAY));
    if (ma->ndim > 1) {
      ArrayInit(ma->array, sizeof(ARRAY), ma->block[0]);
    } else {
      ArrayInit(ma->array, ma->esize, ma->block[0]);
    }
  }
  a = ma->array;
  ndim1 = ma->ndim-1;
  ndim2 = ma->ndim-2;
  for (i = 0; i < ndim1; i++) {
    a = (ARRAY *) ArraySet(a, k[i], NULL);
    if (a->esize == 0) {
      if (i < ndim2) {
	ArrayInit(a, sizeof(ARRAY), ma->block[i+1]);
      } else {
	ArrayInit(a, ma->esize, ma->block[i+1]);
      }
    }
  }
    
  pt = ArraySet(a, k[i], d);
  /*
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.multi += stop-start;
#endif
  */
  return pt;
}

int MultiFree(MULTI *ma, void (*FreeElem)(void *)) {
  if (ma->ndim <= 0) return 0;
  MultiFreeData(ma->array, ma->ndim, FreeElem);
  free(ma->array);
  ma->array = NULL;
  free(ma->block);
  ma->block = NULL;
  ma->ndim = 0;
  return 0;
}

int MultiFreeData(ARRAY *a, int d, void (*FreeElem)(void *)) {
  int i, d1;
  ARRAY *b;
  if (a == NULL) return 0;
  if (d > 1) {
    d1 = d-1;
    for (i = 0; i < a->dim; i++) {
      b = (ARRAY *) ArrayGet(a, i);
      if (b) {
	MultiFreeData(b, d1, FreeElem);
      }
    }
    ArrayFree(a, NULL);
  } else {
    ArrayFree(a, FreeElem);
  }
  return 0;
}

