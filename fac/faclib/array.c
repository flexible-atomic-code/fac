#include "array.h"

static char *rcsid="$Id: array.c,v 1.10 2003/01/22 21:58:02 mfgu Exp $";
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
  char *ct;
  DATA *p;
  /*
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif
  */
 
  if (a->dim == 0) {
    a->data = (DATA *) malloc(sizeof(DATA));
    a->data->dptr = calloc(a->block, a->esize);
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
    p->dptr = calloc(a->block, a->esize);
  }
  
  ct = (char *) p->dptr;
  for (; i > 0; i--) {
    ct += a->esize;
  }
  pt = (void *) ct;

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
** FUNCTION:    ArrayContiguous
** PURPOSE:     Return a 1-d standard c-array contiguous in memory.
** INPUT:       {ARRAY *a},
**              pointer to the array.
** RETURN:      {void *},
**              pointer to the resulting array.
** SIDE EFFECT: 
*/
void *ArrayContiguous(ARRAY *a) {
  void *r, *rp;
  DATA *p;
  int i, m;

  if (a->dim == 0) return NULL;
  m = a->esize*a->block;
  r = malloc(a->esize*a->dim);
  p = a->data;
  i = a->dim;
  rp = r;
  while (1) {
    if (i <= a->block) {
      memcpy(rp, p->dptr, i*a->esize);
      break;
    } else {
      memcpy(rp, p->dptr, m);
      rp = ((char *)rp) + m;
      i -= a->block;
    }
    p = p->next;
  }
  
  return r;
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
  if (p->dptr) {
    free(p->dptr);
  }
  if (p) {
    free(p);
  }
    
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
** RETURN:      {int},
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

/* 
** FUNCTION:    ArrayTrim
** PURPOSE:     Trim the tail of an array to a given length.
** INPUT:       {ARRAY *a},
**              pointer to the array.
**              {int n},
**              length of the final array.
**              {void (*FreeElem)(void *)},
**              a function called before freeing each element.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        if the length of array is <= n, nothing happens.
*/    
int ArrayTrim(ARRAY *a, int n, void (*FreeElem)(void *)) {
  DATA *p;
  void *pt;
  int i;

  if (!a) return 0;
  if (a->dim <= n) return 0;
  
  if (n == 0) {
    ArrayFree(a, FreeElem);
  }

  i = n;
  p = a->data;
  while (i >= a->block) {
    p = p->next;
    i -= a->block;
  }

  if (i == 0) {
    ArrayFreeData(p, a->esize, a->block, FreeElem);
    p = NULL;
  } else {
    if (p->next) {
      ArrayFreeData(p->next, a->esize, a->block, FreeElem);
      p->next = NULL;
    }
    if (p->dptr && FreeElem) {
      pt = ((char *) p->dptr) + i*(a->esize);
      for (; i < a->block; i++) {
	FreeElem(pt);
	pt = ((char *) pt) + a->esize;
      }
    }
  }

  a->dim = n;

  return 0;
}         

/* 
** FUNCTION:    MultiInit
** PURPOSE:     initialize a multi-dimensional array.
** INPUT:       {MULTI *ma},
**              pointer to the multi-dimensional array.
**              {int esize},
**              size of each element in bytes.
**              {int ndim},
**              number of dimensions of the array.
**              {int *block},
**              integer array of length ndim,
**              giving the block size in each dimension.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/    
int MultiInit(MULTI *ma, int esize, int ndim, int *block) {
  int i;
  ma->ndim = ndim;
  ma->esize = esize;
  ma->block = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  for (i = 0; i < ndim; i++) ma->block[i] = block[i];
  ma->array = NULL;
  return 0;
}

/* 
** FUNCTION:    MultiGet
** PURPOSE:     get an element in a multi-dimensional array.
** INPUT:       {MULTI *ma},
**              pointer to the multi-dimensional array.
**              {int *k},
**              integer array of length ndim, 
**              giving the indexes in each dimension.
** RETURN:      {void *},
**              pointer to the element
** SIDE EFFECT: 
** NOTE:        
*/    
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

/* 
** FUNCTION:    MultiSet
** PURPOSE:     Set an element in a multi-dimensional array.
** INPUT:       {MULTI *ma},
**              pointer to the multi-dimensional array.
**              {int *k},
**              integer array of length ndim, 
**              giving the indexes in each dimension.
**              {void *d},
**              pointer to a piece of data to be copied to the array.
** RETURN:      {void *},
**              pointer to the element just set.
** SIDE EFFECT: 
** NOTE:        if d == NULL, returns an uninitialized element.
*/    
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

/* 
** FUNCTION:    MultiFree
** PURPOSE:     Free multi-dimensional array.
** INPUT:       {MULTI *ma},
**              pointer to the multi-dimensional array.
**              {void (*FreeElem)(void *)},
**              a function called before freeing each element.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/    
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

/* 
** FUNCTION:    MultiFreeData
** PURPOSE:     Free the data of a multi-dimensional array.
** INPUT:       {ARRAY *a},
**              pointer to an array, which is the data of MULTI.
**              {int d},
**              the number of dimensions the array contains.
**              {void (*FreeElem)(void *)},
**              a function called before freeing each element.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT: 
** NOTE:        
*/    
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

/*
** the following set of funcitons are a different implementation
** for the MULTI array, which turns out to be much worse than
** than the original one. It is not being used.
*/
int NMultiInit(MULTI *ma, int esize, int ndim, int *block) {
  int i, n;

  ma->ndim = ndim;
  ma->esize = esize;
  ma->block = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  n = 1;
  for (i = 0; i < ndim; i++) {
    ma->block[i] = block[i];
    n *= block[i];
  }
  ma->array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ma->array, sizeof(MDATA), n);
  return 0;
}

static int IndexMatch(unsigned short *i, int *k, int n) {
  int q;
  
  for (q = 0; q < n; q++) {
    if (i[q] != k[q]) return 0;
  }
  return 1;
}

static void IndexCopy(unsigned short *i, int *k, int n) {
  int q;

  for (q = 0; q < n; q++) {
    i[q] = k[q];
  }
}

void *NMultiGet(MULTI *ma, int *k) {
  ARRAY *a;
  MDATA *pt;
  DATA *p;
  int i, j, m;

  a = ma->array;
  p = a->data;
  i = a->dim;
  j = 0;
  while (p) {
    pt = (MDATA *) p->dptr;
    for (m = 0; m < a->block && j < i; j++, m++) {
      if (IndexMatch(pt->index, k, ma->ndim)) {
	return pt->data;
      }
      pt++;
    }
    p = p->next;
  }

  return NULL;
}

void *NMultiSet(MULTI *ma, int *k, void *d) {
  int i, j, m;
  MDATA *pt;
  ARRAY *a;
  DATA *p, *p0;
  
  a = ma->array;
  if (a->dim == 0) {
    a->data = (DATA *) malloc(sizeof(DATA));
    a->data->dptr = calloc(a->block, a->esize);
    a->data->next = NULL;
    pt = (MDATA *) a->data->dptr;
  } else {
    p = a->data;
    i = a->dim;
    j = 0;
    while (p) {
      pt = (MDATA *) p->dptr;
      for (m = 0; m < a->block && j < i; j++, m++) {
	if (IndexMatch(pt->index, k, ma->ndim)) {
	  if (d) {
	    memcpy(pt->data, d, ma->esize);
	  }
	  return pt->data;
	}
	pt++;
      }
      p0 = p;
      p = p->next;
    }  
    if (m == a->block) {
      p0->next = (DATA *) malloc(sizeof(DATA));
      p = p0->next;
      p->dptr = calloc(a->block, a->esize);
      p->next = NULL;
      pt = (MDATA *) p->dptr;
    }
  }

  pt->index = (unsigned short *) malloc(sizeof(unsigned short)*ma->ndim);
  IndexCopy(pt->index, k, ma->ndim);
  pt->data = calloc(1, ma->esize);
  if (d) memcpy(pt->data, d, ma->esize);
  (a->dim)++;

  return pt->data;
}

static int NMultiArrayFreeData(DATA *p, int esize, int block, 
			void (*FreeElem)(void *)) { 
  MDATA *pt;
  int i;
  
  if (p->next) {
    NMultiArrayFreeData(p->next, esize, block, FreeElem);
  }

  if (p->dptr) {
    pt = p->dptr;
    for (i = 0; i < block; i++) {
      free(pt->index);
      if (FreeElem) FreeElem(pt->data);
      free(pt->data);
      pt++;
    }
    free(p->dptr);
  }
  if (p) {
    free(p);
  }
  p = NULL;
  return 0;
}
    
int NMultiFreeData(ARRAY *a, int d, void (*FreeElem)(void *)) {
  if (!a) return 0;
  if (a->dim == 0) return 0;
  NMultiArrayFreeData(a->data, a->esize, a->block, FreeElem);
  a->dim = 0;
  a->data = NULL;
  return 0;
}

int NMultiFree(MULTI *ma, void (*FreeElem)(void *)) {
  if (!ma) return 0;
  if (ma->ndim <= 0) return 0;
  NMultiFreeData(ma->array, ma->ndim, FreeElem);
  free(ma->array);
  ma->array = NULL;
  free(ma->block);
  ma->block = NULL;
  ma->ndim = 0;
  return 0;
}
