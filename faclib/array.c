#include "array.h"

static char *rcsid="$Id: array.c,v 1.7 2001/11/24 21:12:28 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/*************************************************************
  Implementation of module "array"
  
  This module implements a variable length one- and 
  multi-dimensional array.

  Author: M. F. Gu, mfgu@space.mit.edu
**************************************************************/

#ifdef PERFORM_STATISTICS  
static ARRAY_TIMING timing = {0, 0};
/* 
** FUNCTION:    GetArrayTiming
** PURPOSE:     retreive the timing information.
** INPUT:       {ARRAY_TIMING *t},
**              pointer to ARRAY_TIMING struct, holding the result.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/
int GetArrayTiming(ARRAY_TIMING *t) {
  memcpy(t, &timing, sizeof(timing)); 
  return 0;
}
#endif

/* 
** FUNCTION:    ArrayInit
** PURPOSE:     initialize the one-dimensional array.
** INPUT:       {ARRAY *a},
**              pointer to the array to be initialized.
**              {int esize},
**              size of the elements in bytes.
**              {int block},
**              number of elements in one block.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/
int ArrayInit(ARRAY *a, int esize, int block) {
  a->esize = esize;
  a->block = block;
  a->dim = 0;
  a->data = NULL;
  return 0;
}

/* 
** FUNCTION:    ArrayGet
** PURPOSE:     retrieve the i-th element of the array.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {int i},
**              index of the element.
** RETURN:      {void *},
**              pointer to the element. 
**              NULL, if does not exist.
** SIDE EFFECT: 
** NOTE:        
*/
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

/* 
** FUNCTION:    ArraySet
** PURPOSE:     set the i-th element.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {int i},
**              index of the element.
**              {void *d},
**              pointer to the data to be copied.
** RETURN:      {void *},
**              pointer to the element.
** SIDE EFFECT: 
** NOTE:        if d == NULL, this function simply retrieve the
**              i-th element. if the element does not exist,
**              an empty one is created.
*/
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

/* 
** FUNCTION:    ArrayAppend
** PURPOSE:     append an element to the array
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {void *d},
**              data to be appened.
** RETURN:      {void *},
**              pointer to the appended element.
** SIDE EFFECT: 
** NOTE:        
*/
void *ArrayAppend(ARRAY *a, void *d) {
  int i;  
  i = a->dim;
  return ArraySet(a, i, d);
}

/* 
** FUNCTION:    ArrayFreeData
** PURPOSE:     free the data stored in the array.
** INPUT:       {DATA *p},
**              pointer to the data to be freed
**              {int esize},
**              size of the element in bytes.
**              {int block},
**              number of elements in one block.
**              {void (*FreeElem)(void *)},
**              a function called before freeing the data.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        this function calls itself recursively.
*/
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

/* 
** FUNCTION:    ArrayFree
** PURPOSE:     deinitialize the array.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {void (*FreeElem)(void *)},
**              a function called before freeing each element.
** RETURN:      {int 0},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/    
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

