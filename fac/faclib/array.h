#ifndef _ARRAY_H_
#define _ARRAY_H_ 1

/*************************************************************
  Header of module "array"
  
  This module implements a variable length one- and 
  multi-dimensional array.

  Author: M. F. Gu, mfgu@space.mit.edu
**************************************************************/

/* 
<** The following format is used for documenting the source **>
*/

/* documenting a struct */
/*
** STRUCT:      
** PURPOSE:     
** FIELDS:      
** NOTE:        
*/

/* documenting a function */
/* 
** FUNCTION:    
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a macro function */
/* 
** MACRO:       
** PURPOSE:     
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/

/* documenting a global, static varialbe or a macro constant */
/*
** VARIABLE:    
** TYPE:        
** PURPOSE:     
** NOTE:        
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "global.h"

/*
** STRUCT:      DATA
** PURPOSE:     data structure of the array elements
** FIELDS:      {void *dptr},
**              pointer to a block of data.
**              {DATA *next},
**              pointer to the next block.
** NOTE:        
*/
typedef struct _DATA_ {
  void *dptr;
  struct _DATA_ *next;
} DATA;

/*
** STRUCT:      ARRAY
** PURPOSE:     a one-dimensional array.
** FIELDS:      {short esize},
**              the size of each element in bytes.
**              {short block},
**              number of elements in each block.
**              {int dim},
**              the size of the array.
** NOTE:        
*/
typedef struct _ARRAY_ {
  short esize;
  short block;
  int   dim;
  DATA  *data;
} ARRAY;

/*
** STRUCT:      MULTI
** PURPOSE:     a multi-dimensional array.
** FIELDS:      {short ndim},
**              the dimension of the array.
**              {short esize},
**              size of each array element in bytes.
**              {short *block},
**              number of elements in each block for each dimension.
**              {ARRAY *array},
**              the multi-dimensional array is implemented as array 
**              of arrays. 
** NOTE:        
*/
typedef struct _MULTI_ {
  short ndim;
  short esize;
  short *block;
  ARRAY *array;
} MULTI;


#ifdef PERFORM_STATISTICS
/*
** STRUCT:      ARRAY_TIMING
** PURPOSE:     timing informaiton for the modue "array".
** FIELDS:      {clock_t array, multi},
**              time spent in ARRAY and MULTI operations.
** NOTE:        
*/
typedef struct _ARRAY_TIMING_ {
  clock_t array;
  clock_t multi;
} ARRAY_TIMING;

int   GetArrayTiming(ARRAY_TIMING *t);
#endif

int   ArrayInit(ARRAY *a, int esize, int block);
void *ArrayGet(ARRAY *a, int i);
void *ArraySet(ARRAY *a, int i, void *d);
void *ArrayAppend(ARRAY *a, void *d);
int   ArrayTrim(ARRAY *a, int n, void(*FreeElem)(void *));
int   ArrayFree(ARRAY *a, void (*FreeElem)(void *));
int   ArrayFreeData(DATA *p, int esize, int block, void (*FreeElem)(void *));

int   MultiInit(MULTI *ma, int esize, int ndim, int *block);
void *MultiGet(MULTI *ma, int *k);
void *MultiSet(MULTI *ma, int *k, void *d);
int   MultiFree(MULTI *ma, void (*FreeElem)(void *));
int   MultiFreeData(ARRAY *a, int d, void (*FreeElem)(void *));

#endif
