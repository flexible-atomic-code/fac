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

#include "array.h"
#include "mpiutil.h"

static char *rcsid="$Id: array.c,v 1.17 2005/07/20 19:43:19 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/*************************************************************
  Implementation of module "array"
  
  This module implements a variable length one- and 
  multi-dimensional array.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

static double _maxsize = -1;
static double _totalsize = 0;
static double _overheadsize = 0;
static ARRAY *_multistats = NULL;

//#undef SetLock
//#define SetLock(x) SetLockWT((x))
void InitMultiStats(void) {
  if (_multistats == NULL) {
    _multistats = malloc(sizeof(ARRAY));
    ArrayInit(_multistats, sizeof(MULTI *), 256);
  }
}

void ReportMultiStats(void) {
  int i;
  if (_multistats == NULL) return;
  //if (MyRankMPI() != 0) return;
  for (i = 0; i < _multistats->dim; i++) {
    MULTI **pma = (MULTI **) ArrayGet(_multistats, i);
    if (pma == NULL) continue;
    MULTI *ma = *pma;
    if (ma == NULL) continue;
    if (ma->numelem > 0) {
      MPrintf(-1, "idx=%d, id=%s, nd=%d, hs=%d, ne=%d, me=%d, ts=%g, os=%g, ms=%g, c=%d, ch=%g, isize=%d, esize=%d, lock=%x\n", i, ma->id, ma->ndim, ma->hsize, ma->numelem, ma->maxelem, ma->totalsize, ma->overheadsize, ma->maxsize, ma->clean_flag, ma->cth, ma->isize, ma->esize, ma->lock);
    }
  }
}

void RemoveMultiLocks(void) {
  int i, j;
  if (_multistats == NULL) return;
  if (MyRankMPI() != 0) return;
  for (i = 0; i < _multistats->dim; i++) {
    MULTI **pma = (MULTI **) ArrayGet(_multistats, i);
    if (pma == NULL) continue;
    MULTI *ma = *pma;
    if (ma == NULL) continue;
    if (ma->lock) {
      DestroyLock(ma->lock);
      free(ma->lock);
      ma->lock = NULL;
    }
    for (j = 0; j < ma->hsize; j++) {
      ARRAY *a = &ma->array[j];
      if (a->lock) {
	DestroyLock(a->lock);
	free(a->lock);
	a->lock = NULL;
      }
    }
  }
}
      
void InitIntData(void *p, int n) {
  int *d;
  int i;
  
  d = (int *) p;
  for (i = 0; i < n; i++) {
    d[i] = 0;
  }
}

void InitDoubleData(void *p, int n) {
  double *d;
  int i;
  
  d = (double *) p;
  for (i = 0; i < n; i++) {
    d[i] = 0;
  }
}

void InitPointerData(void *p, int n) {
  void **d;
  int i;

  d = (void **) p;
  for (i = 0; i < n; i++) {
    d[i] = NULL;
  }
}

void InitArrayData(void *p, int n) {
  ARRAY *d;
  int i;

  d = (ARRAY *) p;
  for (i = 0; i < n; i++) {
    d[i].dim = 0;
    d[i].esize = 0;
  }
}

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
  a->bsize = ((int)esize)*((int)block);
  a->dim = 0;
  a->data = NULL;
#if USE_MPI == 2
  a->lock = (LOCK *) malloc(sizeof(LOCK));
  if (0 != InitLock(a->lock)) {
    printf("cannot InitLock in ArrayInit\n");
    free(a->lock);
    a->lock = NULL;
    Abort(1);
  }
#else  
  a->lock = NULL;
#endif
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
**              {void (*InitData)(void *, int)},
**              a function to be called to initialize the data
**              when first created.
** RETURN:      {void *},
**              pointer to the element.
** SIDE EFFECT: 
** NOTE:        if d == NULL, this function simply retrieve the
**              i-th element. if the element does not exist,
**              an empty one is created.
*/
void *ArraySet(ARRAY *a, int i, void *d, 
	       void (*InitData)(void *, int)) {
  void *pt;
  char *ct;
  DATA *p;
 
  if (a->dim == 0) {
    a->data = (DATA *) malloc(sizeof(DATA));
    a->data->dptr = malloc(a->bsize);
    if (InitData) InitData(a->data->dptr, a->block);
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
    p->dptr = malloc(a->bsize);
    if (InitData) InitData(p->dptr, a->block);
  }
  
  ct = (char *) p->dptr;
  for (; i > 0; i--) {
    ct += a->esize;
  }
  pt = (void *) ct;
  
  if (d) memcpy(pt, d, a->esize);
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
void *ArrayAppend(ARRAY *a, void *d, 
		  void (*InitData)(void *, int)) {
  int i;  
  i = a->dim;
  return ArraySet(a, i, d, InitData);
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

  if (p == NULL) return 0;

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
  free(p);
  p = NULL;

  return 0;
}

int ArrayFreeLock(ARRAY *a, void (*FreeElem)(void *)) {
  if (a == NULL) return 0;
  ArrayFree(a, FreeElem);
  if (a->lock) {
    DestroyLock(a->lock);
    free(a->lock);
    a->lock = NULL;
  }
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
    return 0;
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
** FUNCTION:    SMultiInit
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
int SMultiInit(MULTI *ma, int esize, int ndim, int *block, char *id) {
  int i;
  if (id != NULL) {
    strncpy(ma->id, id, MULTI_IDLEN-1);
  } else {
    ma->id[0] = '\0';
  }
  ma->maxsize = -1;
  ma->totalsize = 0;
  ma->cth = 0;
  ma->clean_mode = -1;
  ma->clean_flag = 0;
  ma->ndim = ndim;
  ma->esize = esize;
  ma->block = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  for (i = 0; i < ndim; i++) ma->block[i] = block[i];
  ma->array = NULL;
  return 0;
}

/* 
** FUNCTION:    SMultiGet
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
void *SMultiGet(MULTI *ma, int *k, LOCK **lock) {
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
** FUNCTION:    SMultiSet
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
void *SMultiSet(MULTI *ma, int *k, void *d, LOCK **lock,
		void (*InitData)(void *, int),
		void (*FreeElem)(void *)) {
  ARRAY *a;
  void *pt;
  int i, ndim1, ndim2;

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
    a = (ARRAY *) ArraySet(a, k[i], NULL, InitArrayData);
    if (a->esize == 0) {
      if (i < ndim2) {
	ArrayInit(a, sizeof(ARRAY), ma->block[i+1]);
      } else {
	ArrayInit(a, ma->esize, ma->block[i+1]);
      }
    }
  }
  
  pt = ArraySet(a, k[i], d, InitData);
  return pt;
}

/* 
** FUNCTION:    SMultiFreeDataOnly
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
int SMultiFreeDataOnly(ARRAY *a, int d, void (*FreeElem)(void *)) {
  int i, d1;
  ARRAY *b;
  if (a == NULL) return 0;
  if (d > 1) {
    d1 = d-1;
    for (i = 0; i < a->dim; i++) {
      b = (ARRAY *) ArrayGet(a, i);
      if (b) {
	SMultiFreeDataOnly(b, d1, FreeElem);
      }
    }
    ArrayFree(a, NULL);
  } else {
    ArrayFree(a, FreeElem);
  }
  return 0;
}

int SMultiFreeData(MULTI *ma, void (*FreeElem)(void *)) {
  return SMultiFreeDataOnly(ma->array, ma->ndim, FreeElem);
}

/* 
** FUNCTION:    SMultiFree
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
int SMultiFree(MULTI *ma, void (*FreeElem)(void *)) {
  if (ma->ndim <= 0) return 0;
  SMultiFreeData(ma, FreeElem);
  free(ma->array);
  ma->array = NULL;
  free(ma->block);
  ma->block = NULL;
  ma->ndim = 0;
  return 0;
}


/* hash table size 2^NHASH */
typedef unsigned long int ub4;
typedef struct _MDATA_ {
  int *index;
  LOCK *lock;
  void *data;
} MDATA;

void InitMDataData(void *p, int n) {
  MDATA *d;
  int i;
  
  d = (MDATA *) p;
  for (i = 0; i < n; i++) {
    d[i].index = NULL;
    d[i].data = NULL;
    d[i].lock = NULL;
  }
}

#define HashSize(n) ((ub4)1<<(((n)/2)+16))
#define HashMask(n) (HashSize(n)-1)

#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

#define mix(a,b,c) \
{ \
  a -= c;  a ^= rot(c, 4);  c += b; \
  b -= a;  b ^= rot(a, 6);  a += c; \
  c -= b;  c ^= rot(b, 8);  b += a; \
  a -= c;  a ^= rot(c,16);  c += b; \
  b -= a;  b ^= rot(a,19);  a += c; \
  c -= b;  c ^= rot(b, 4);  b += a; \
}

#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}

#define Mix(a, b, c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}
  
static inline int Hash2a(int *id, ub4 length, ub4 initval, int n, int m) {
  ub4 a, b, c, len, i, *k, kd[32];

  k = kd;
  for (i = 0; i < length; i++) k[i] = id[i];
  
  len = length;
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;           /* the previous hash value */
  
  /*---------------------------------------- handle most of the key */
  while (len >= 3) {
    a += k[0];
    b += k[1];
    c += k[2];
    Mix(a,b,c);
    k += 3; len -= 3;
  }

  /*-------------------------------------- handle the last 2 ub4's */
  c += length;
  switch(len) {
    /* c is reserved for the length */
  case 2 : b+=k[1];
  case 1 : a+=k[0];
    /* case 0: nothing left to add */
  }
  Mix(a,b,c);
  /*-------------------------------------------- report the result */
  return (int) (c & m);
}

static inline int Hash2(int *k, ub4 length, ub4 initval, int n, int m) {
  ub4 a,b,c;
  
  /* Set up the internal state */
  a = b = c = 0xdeadbeef + (((ub4)length)<<2) + initval;
  
  /*------------------------------------------------- handle most of the key */
  while (length > 3) {
    a += k[0];
    b += k[1];
    c += k[2];
    mix(a,b,c);
    length -= 3;
    k += 3;
  }
  
  /*------------------------------------------- handle the last 3 uint32_t's */
  switch(length)                     /* all the case statements fall through */
    { 
    case 3 : c+=k[2];
    case 2 : b+=k[1];
    case 1 : a+=k[0];
      final(a,b,c);
    case 0:     /* case 0: nothing left to add */
      break;
    }
  /*------------------------------------------------------ report the result */
  return c&m;
}

void AddMultiSize(MULTI *ma, int size) {
  ma->totalsize += size;
  _totalsize += size;
}

void LimitMultiSize(MULTI *ma, double r) {
  if (ma == NULL) {
    _maxsize = r;
  } else {
    ma->maxsize = r;
  }
}

double TotalSize() {
#if PMALLOC_CHECK == 2
  return (double) msize();
#else
  return _totalsize;
#endif
}

double TotalArraySize() {
  return _totalsize;
}

inline int IdxCmp(int *i0, int *i1, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (i0[i] != i1[i]) return 1;
  }
  return 0;
}

int NMultiInit(MULTI *ma, int esize, int ndim, int *block, char *id) {
  int i, n, s;
  if (id != NULL) {
    strncpy(ma->id, id, MULTI_IDLEN-1);
  } else {
    ma->id[0] = '\0';
  }
  ma->maxsize = -1;
  ma->totalsize = 0;
  ma->cth = 0;
  ma->clean_mode = -1;
  ma->clean_flag = 0;
  ma->ndim = ndim;
  ma->isize = sizeof(int)*ndim;
  ma->esize = esize;
  s = sizeof(unsigned short)*ndim;
  ma->block = (unsigned short *) malloc(s);
  ma->overheadsize += s;
  _overheadsize += s;
  n = HashSize(ma->ndim);
  ma->hsize = n;
  ma->hmask = ma->hsize-1;
  s = sizeof(ARRAY)*n;
  ma->array = (ARRAY *) malloc(s);
  ma->overheadsize += s;
  _overheadsize += s;
  for (i = 0; i < n; i++) {
    ArrayInit(&(ma->array[i]), sizeof(MDATA), 8);
  }
#if USE_MPI == 2
  ma->lock = (LOCK *) malloc(sizeof(LOCK));
  if (0 != InitLock(ma->lock)) {
    printf("cannot InitLock0 in NMultiInit: %s\n", ma->id);
    free(ma->lock);
    ma->lock = NULL;
    Abort(1);
  }

  ma->clean_lock = (LOCK *) malloc(sizeof(LOCK));
  if (0 != InitLock(ma->clean_lock)) {
    printf("cannot InitLock1 in NMultiInit: %s\n", ma->id);
    free(ma->clean_lock);
    ma->clean_lock = NULL;
    Abort(1);
  }
#else
  ma->lock = NULL;
  ma->clean_lock = NULL;
#endif
  if (_multistats != NULL && ma->id[0]) {
    ArrayAppend(_multistats, &ma, InitPointerData);    
  }
  ma->iset = 0;
  return 0;
}

void *NMultiGet(MULTI *ma, int *k, LOCK **lock) {
  ARRAY *a;
  MDATA *pt;
  DATA *p;
  int i, j, m, h;

  h = Hash2(k, ma->ndim, 0, ma->ndim, ma->hmask);
  a = &(ma->array[h]);
  p = a->data;
  i = a->dim;
  j = 0;
  while (p) {
    pt = (MDATA *) p->dptr;
    for (m = 0; m < a->block && j < i; j++, m++) {
      if (IdxCmp(pt->index, k, ma->ndim) == 0) {
	if (lock) *lock = pt->lock;
	return pt->data;
      }
      pt++;
    }
    p = p->next;
  }
  
  return NULL;
}

void *NMultiSet(MULTI *ma, int *k, void *d, LOCK **lock,
		void (*InitData)(void *, int),
		void (*FreeElem)(void *)) {
  int i, j, m, h, size, locked = 0, clocked = 0;
  MDATA *pt;
  ARRAY *a;
  DATA *p, *p0;
  int myrank = MyRankMPI()+1;
  int cleanmode;
  cleanmode = ma->clean_mode;
  if (cleanmode < 0) {
    if (_maxsize > 0 && ma->cth > 0) {
      double ts = TotalSize();
      double ats = TotalArraySize();
      if (ts >= _maxsize && ma->totalsize > ma->cth*ats) {
	cleanmode = 1;
      }
    } else if (ma->maxsize > 0 && ma->totalsize >= ma->maxsize) {
      cleanmode = 0;
    } else if (ma->clean_flag > 0) {
      cleanmode = 0;
    }
  }
  if (cleanmode >= 0) {
    if (ma->clean_lock) {
      SetLock(ma->clean_lock);
      clocked = 1;
    }
    if (ma->iset == 0) {
      ma->clean_mode = cleanmode;
      NMultiFreeData(ma, FreeElem);
    }
  }
#pragma omp atomic
  ma->iset += myrank;  
  if (clocked) {
    ReleaseLock(ma->clean_lock);
  }

  h = Hash2(k, ma->ndim, 0, ma->ndim, ma->hmask);
  a = &(ma->array[h]);
  locked = 0;
  if (a->dim == 0) {
    if (a->lock) {
      SetLock(a->lock);
      locked = 1;
    }
    if (a->dim == 0) {
      size = sizeof(DATA);
      a->data = (DATA *) malloc(size);
      a->data->dptr = malloc(a->bsize);
      size += a->bsize;
      ma->totalsize += size;
      _totalsize += size;
      InitMDataData(a->data->dptr, a->block);
      a->data->next = NULL;
      pt = (MDATA *) a->data->dptr;
    }
  } 
  if (a->dim > 0) {
    p = a->data;
    i = a->dim;
    j = 0;
    p0 = p;
    pt = (MDATA *) p->dptr;
    while (p && j < i) {
      for (m = 0; m < a->block && j < i; j++, m++) {
	if (IdxCmp(pt->index, k, ma->ndim) == 0) {
	  if (d) {
	    memcpy(pt->data, d, ma->esize);
	  }
	  if (lock) *lock = pt->lock;
	  if (locked) {
	    ReleaseLock(a->lock);
	  }	  
	  return pt->data;
	}
	pt++;
      }
      if (m == a->block) {
	p0 = p;
	p = p->next;
	if (p) {
	  pt = (MDATA *) p->dptr;
	} else {
	  break;
	}
      }
    }  
    if (!locked && a->lock) {
      SetLock(a->lock);
      locked = 1;
    }
    if (a->dim > i) {
      if (m == a->block) {
	p = p0->next;
	pt = (MDATA *) p->dptr;
	m = 0;
      }
      while (p && j < a->dim) {
	for (; m < a->block && j < a->dim; j++, m++) {
	  if (IdxCmp(pt->index, k, ma->ndim) == 0) {
	    if (d) {
	      memcpy(pt->data, d, ma->esize);
	    }
	    if (lock) *lock = pt->lock;
	    if (locked) {
	      ReleaseLock(a->lock);
	    }
	    return pt->data;
	  }
	  pt++;
	}
	if (m == a->block) {
	  p0 = p;
	  p = p->next;
	  if (p) {
	    pt = (MDATA *) p->dptr;
	    m = 0;
	  } else {
	    break;
	  }
	}
      }  
    }
    if (m == a->block) {	
      size = sizeof(DATA);
      p0->next = (DATA *) malloc(size);
      p = p0->next;
      p->dptr = malloc(a->bsize);
      size += a->bsize;
      ma->totalsize += size;
      _totalsize += size;
      InitMDataData(p->dptr, a->block);
      p->next = NULL;
      pt = (MDATA *) p->dptr;
    }
  }

  size = sizeof(LOCK);
#if USE_MPI == 2
  pt->lock = (LOCK *) malloc(sizeof(LOCK));
  if (0 != InitLock(pt->lock)) {
    free(pt->lock);
    pt->lock = NULL;
  }
#else
  pt->lock = NULL;
#endif
  pt->data = malloc(ma->esize);
  size += ma->esize + ma->isize;
  ma->totalsize += size;
  ma->numelem++;
  _totalsize += size;
  if (InitData) InitData(pt->data, 1);
  if (d) memcpy(pt->data, d, ma->esize);
  if (lock) *lock = pt->lock;  
  int *idx = (int *) malloc(ma->isize);
  memcpy(idx, k, ma->isize);
  pt->index = idx;
  (a->dim)++;  
  if (locked) {
    ReleaseLock(a->lock);
  }
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
      pt->index = NULL;
      if (pt->lock) {
	DestroyLock(pt->lock);
	free(pt->lock);
	pt->lock = NULL;
      }
      if (FreeElem && pt->data) FreeElem(pt->data);
      free(pt->data);
      pt->data = NULL;
      pt++;
    }
    free(p->dptr);
    p->dptr = NULL;
  }
  if (p) {
    free(p);
    p = NULL;
  }
  return 0;
}

void SetMultiCleanFlag(MULTI *ma) {
  if (ma->lock) {
    SetLock(ma->lock);
  }
  //MPrintf(-1, "SMC: %s %d\n", ma->id, ma->clean_flag);
  if (ma->clean_flag <= 0 &&
      ma->cth > 0 &&
      ma->totalsize > ma->cth*TotalArraySize()) {
    ma->clean_flag = 1;
  }
  if (ma->lock) {
    ReleaseLock(ma->lock);
  }
}

int NMultiFreeDataOnly(ARRAY *a, void (*FreeElem)(void *)) {
  if (!a) return 0;
  if (a->dim == 0) return 0;
  NMultiArrayFreeData(a->data, a->esize, a->block, FreeElem);
  a->dim = 0;
  a->data = NULL;
  return 0;
}

int NMultiFreeData(MULTI *ma, void (*FreeElem)(void *)) {
  ARRAY *a;
  int i;
#pragma omp flush
  if (ma->lock) SetLock(ma->lock);
  int clean = 1;
  if (ma->clean_mode == 0) {
    if (ma->totalsize < ma->maxsize && ma->clean_flag <= 0) clean = 0;
    if (clean) {
      MPrintf(-1,
	      "clean0 %s t=%g o=%g m=%g tt=%g to=%g tm=%g c=%d ch=%g\n",
	      ma->id, ma->totalsize, ma->overheadsize, ma->maxsize,	    
	      _totalsize, _overheadsize, _maxsize, ma->clean_flag, ma->cth);
    }
  } else if (ma->clean_mode == 1) {
    double ts = TotalSize();
    double ats = TotalArraySize();
    if (ts < _maxsize || ma->totalsize <= ma->cth*ats) clean = 0;
    if (clean) {
      MPrintf(-1,
	      "clean1: %s t=%g o=%g m=%g tt=%g to=%g tm=%g c=%d ch=%g\n",
	      ma->id, ma->totalsize, ma->overheadsize, ma->maxsize,
	      _totalsize, _overheadsize, _maxsize, ma->clean_flag, ma->cth);
    }
  } else {
    if (ma->totalsize <= 0 && ma->clean_flag <= 0) clean = 0;
  }
  if (clean) {
    ma->clean_thread = MyRankMPI();
    if (ma->iset > 0 && ma->clean_mode >= 0) {
      printf("invalid clean with iset: %s: %d %lud\n",
	     ma->id, ma->clean_thread, ma->iset);
      Abort(1);
    }
    for (i = 0; i < ma->hsize; i++) {
      a = &(ma->array[i]);
      if (a->lock) SetLock(a->lock);
      NMultiFreeDataOnly(a, FreeElem);
      if (a->lock) ReleaseLock(a->lock);
    }
    _totalsize -= ma->totalsize;
    ma->totalsize = 0;
    ma->numelem = 0;
    ma->clean_flag = 0;
  }
  ma->clean_mode = -1;
#pragma omp flush
  if (ma->lock) ReleaseLock(ma->lock);
  return 0;
}

int NMultiFree(MULTI *ma, void (*FreeElem)(void *)) {
  if (!ma) return 0;
  if (ma->ndim <= 0) return 0;
  NMultiFreeData(ma, FreeElem);
  free(ma->array);
  ma->array = NULL;
  free(ma->block);
  ma->block = NULL;
  ma->ndim = 0;
  ma->iset = 0;
  return 0;
}

int MultiFreeLock(MULTI *ma, void (*FreeElem)(void *)) {
  MultiFree(ma, FreeElem);
  if (ma->lock) {
    DestroyLock(ma->lock);
    free(ma->lock);
    ma->lock = NULL;
  }
  if (ma->clean_lock) {
    DestroyLock(ma->clean_lock);
    free(ma->clean_lock);
    ma->clean_lock = NULL;
  }
  return 0;
}

int MMultiInit(MULTI *ma, int esize, int ndim, int *block, char *id) {
  int i, n;
  if (id == NULL) {
    strncpy(ma->id, id, MULTI_IDLEN-1);
  } else {
    ma->id[0] = '\0';
  }
  ma->maxsize = -1;
  ma->totalsize = 0;
  ma->cth = 0;
  ma->clean_mode = -1;
  ma->clean_flag = 0;
  ma->ndim = ndim;
  ma->ndim1 = ndim-1;
  ma->isize = sizeof(unsigned short)*ndim;
  ma->esize = esize;
  ma->block = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  ma->iidx = (int *) malloc(sizeof(int)*ndim);
  ma->sidx = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  ma->ridx = (unsigned short *) malloc(sizeof(unsigned short)*ndim);
  ma->iblock = (int *) malloc(sizeof(int)*ndim);
  for (i = 0; i < ndim; i++) {
    ma->block[i] = block[i];
    if (i > 0) {
      ma->iblock[i] = block[i]*ma->iblock[i-1];
    } else {
      ma->iblock[i] = block[i];
    }
  }
  
  n = HashSize(ma->ndim);
  ma->hsize = n;
  ma->ia = (ARRAY *) malloc(sizeof(ARRAY)*n);
  ma->da = (ARRAY *) malloc(sizeof(ARRAY)*n);
  for (i = 0; i < n; i++) {
    ArrayInit(&(ma->ia[i]), ma->isize, 16);
    ArrayInit(&(ma->da[i]), ma->esize, ma->iblock[ma->ndim1]);
  }
  ma->array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(ma->array, ma->esize, ma->iblock[ma->ndim1]);
  return 0;
}

void MMultiIndex(MULTI *ma, int *k) {
  int i;
  ma->isf = 1;
  for (i = 0; i < ma->ndim; i++) {
    if (k[i] >= ma->block[i]) {
      ma->iidx[i] = k[i]/ma->block[i];
      ma->sidx[i] = (unsigned short) (ma->iidx[i]);
      ma->ridx[i] = (unsigned short) (k[i]%ma->block[i]);
      ma->isf = 0;
    } else {
      ma->iidx[i] = 0;
      ma->sidx[i] = 0;
      ma->ridx[i] = (unsigned short) k[i];
    }
  }  
  ma->aidx = ma->ridx[0];
  for (i = 1; i < ma->ndim; i++) {
    ma->aidx += ma->ridx[i]*ma->iblock[i-1];
  }
}

void *MMultiGet(MULTI *ma, int *k, LOCK **lock) {
  int i, j, m, h;
  ARRAY *a, *da;
  char *pt;
  DATA *p;
  
  MMultiIndex(ma, k);
  if (ma->isf) {
    return ArrayGet(ma->array, ma->aidx);
  }

  h = Hash2(ma->iidx, ma->ndim, 0, ma->ndim, ma->hmask);
  a = &(ma->ia[h]);
  da = &(ma->da[h]);
  p = a->data;
  i = a->dim;
  j = 0;
  while (p) {
    pt = (char *) p->dptr;
    for (m = 0; m < a->block && j < i; j++, m++) {
      if (bcmp(pt, ma->sidx, ma->isize) == 0) {
	h = ma->aidx + j*ma->iblock[ma->ndim1];
	return ArrayGet(da, h);
      }
      pt += ma->isize;
    }
    p = p->next;
  }
  return NULL;
}

void *MMultiSet(MULTI *ma, int *k, void *d, LOCK **lock,
		void (*InitData)(void *, int),
		void (*FreeElem)(void *)) {
  int i, j, m, h;
  ARRAY *a, *da;
  char *pt;
  DATA *p, *p0;

  MMultiIndex(ma, k);
  if (ma->isf) {
    return ArraySet(ma->array, ma->aidx, d, InitData);
  }
  
  h = Hash2(ma->iidx, ma->ndim, 0, ma->ndim, ma->hmask);
  a = &(ma->ia[h]);
  da = &(ma->da[h]);
  if (a->dim == 0) {
    a->data = (DATA *) malloc(sizeof(DATA));
    a->data->dptr = malloc(a->bsize);
    pt = (char *) a->data->dptr;
    for (i = 0; i < a->bsize; i++, pt++) {
      *pt = 0xFF;
    }
    a->data->next = NULL;
    pt = (char *) a->data->dptr;
  } else {
    p = a->data;
    i = a->dim;
    j = 0;
    while (p) {
      pt = (char *) p->dptr;
      for (m = 0; m < a->block && j < i; j++, m++) {
	if (bcmp(pt, ma->sidx, ma->isize) == 0) {
	  h = ma->aidx + j*ma->iblock[ma->ndim1];
	  return ArraySet(da, h, d, InitData);
	}      
	pt += ma->isize;
      }
      p0 = p;
      p = p->next;
    }

    if (m == a->block) {
      p0->next = (DATA *) malloc(sizeof(DATA));
      p = p0->next;
      p->dptr = malloc(a->bsize);
      pt = (char *) p->dptr;
      for (i = 0; i < a->bsize; i++, pt++) {
	*pt = 0xFF;
      }
      p->next = NULL;
      pt = (char *) p->dptr;
    }
  }

  memcpy(pt, ma->sidx, ma->isize);
  j = ma->aidx + a->dim*ma->iblock[ma->ndim1];
  a->dim++;
  return ArraySet(da, j, d, InitData);
}

int MMultiFreeData(MULTI *ma, void (*FreeElem)(void *)) {
  ARRAY *a;
  int i;
  if (!ma) return 0;
  if (ma->ndim == 0) return 0;
  
  a = ma->array;
  ArrayFreeData(a->data, a->esize, a->block, FreeElem);
  a->dim = 0;
  a->data = NULL;
  for (i = 0; i < ma->hsize; i++) {
    a = &(ma->ia[i]);
    ArrayFreeData(a->data, a->esize, a->block, NULL);
    a->dim = 0;
    a->data = NULL;
    a = &(ma->da[i]);
    ArrayFreeData(a->data, a->esize, a->block, FreeElem);
    a->dim = 0;
    a->data = NULL;
  }
  return 0;
}

int MMultiFree(MULTI *ma, void (*FreeElem)(void *)) {
  if (!ma) return 0;
  if (ma->ndim <= 0) return 0;
  MMultiFreeData(ma, FreeElem);
  free(ma->array);  
  ma->array = NULL;
  free(ma->ia);
  ma->ia = NULL;
  free(ma->da);
  ma->da = NULL;
  free(ma->block);
  free(ma->iblock);
  free(ma->iidx);
  free(ma->sidx);
  free(ma->ridx);
  ma->block = NULL;
  ma->ndim = 0;
  return 0;
}

void InitIdxAry(IDXARY *ia, int n, int *d) {
  int k;
  ia->n = n;
  if (n == 0) {
    ia->d = NULL;
    ia->m0 = 0;
    ia->m1 = 0;
    ia->m = 0;
    ia->i = NULL;
    return;
  }
  ia->d = d;
  ia->m0 = ia->d[0];
  ia->m1 = ia->d[0];
  for (k = 1; k < ia->n; k++) {
    if (ia->m0 > ia->d[k]) ia->m0 = ia->d[k];
    if (ia->m1 < ia->d[k]) ia->m1 = ia->d[k];
  }

  ia->m = 1+ia->m1-ia->m0;
  ia->i = malloc(sizeof(int)*ia->m);
  for (k = 0; k < ia->m; k++) {
    ia->i[k] = -1;
  }
  for (k = 0; k < ia->n; k++) {
    ia->i[ia->d[k]-ia->m0] = k;
  }
}

int IdxGet(IDXARY *ia, int d) {
  if (d < ia->m0) return -1;
  if (d > ia->m1) return -2;
  return ia->i[d - ia->m0];
}
  
void FreeIdxAry(IDXARY *ia, int md) {
  if (md == 0) {
    if (ia->n > 0) free(ia->d);
    if (ia->m > 0) free(ia->i);
    ia->n = 0;
    ia->m = 0;
    return;
  }
  if (md == 1) {
    if (ia->n > 0) free(ia->d);
    ia->n = 0;
    return;
  }
  if (md == 2) {
    if (ia->m > 0) free(ia->i);
    ia->m = 0;
    return;
  }  
}
int Bisect(void *p0, int n, int m, void *p,
	   int (*comp)(const void *, const void *)) {
  int i, i0, i1, k;

  if (n == 0) return -1;
  i0 = 0;
  i1 = n-1;
  while (i1-i0 > 1) {
    i = (i0+i1)/2;
    k = comp(p0, p+m*i);
    if (k == 0) return i;
    else if (k < 0) i1 = i;
    else i0 = i;
  }
  
  k = comp(p0, p);
  if (k == 0) return i0;
  k = comp(p0, p+(n-1)*m);
  if (k == 0) return i1;
  
  return -1;
}
  
int IBisect(int b, int n, int *a) {
  int i, i0, i1;

  if (n == 0) return -1;
  i0 = 0;
  i1 = n - 1;
  while (i1 - i0 > 1) {
    i = (i0 + i1)/2;
    if (b == a[i]) return i;
    else if (b < a[i]) i1 = i;
    else i0 = i;
  }
  
  if (b == a[i0]) return i0;
  else if (b == a[i1]) return i1;
  else return -1;
}

