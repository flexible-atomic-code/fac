#ifndef _ARRAY_H_
#define _ARRAY_H_ 1

/*************************************************************
  Header of module "array"
  
  This module implements a variable length one- and 
  multi-dimensional array.

  Author: M. F. Gu, mfgu@stanford.edu
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
#include <stddef.h>
#include <time.h>

#define USE_NMULTI 1
/* choose MULTI implementation */
#ifdef USE_NMULTI
#define MultiInit NMultiInit
#define MultiGet NMultiGet
#define MultiSet NMultiSet
#define MultiFreeData NMultiFreeData
#define MultiFree NMultiFree
#else
#define MultiInit SMultiInit
#define MultiGet SMultiGet
#define MultiSet SMultiSet
#define MultiFreeData SMultiFreeData
#define MultiFree SMultiFree
#endif /*USE_NMULTI*/


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
  unsigned short esize;
  unsigned short block;
  int bsize;
  int dim;
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
  int numelem, maxelem;
  unsigned short ndim;
  unsigned short isize;
  unsigned short esize;
  unsigned short *block;
  ARRAY *array;
} MULTI;

int   ArrayInit(ARRAY *a, int esize, int block);
void *ArrayGet(ARRAY *a, int i);
void *ArraySet(ARRAY *a, int i, void *d, 
	       void (*InitData)(void *, int));
void *ArrayContiguous(ARRAY *a);
void *ArrayAppend(ARRAY *a, void *d, 
		  void (*InitData)(void *, int));
int   ArrayTrim(ARRAY *a, int n, 
		void(*FreeElem)(void *));
int   ArrayFree(ARRAY *a, void (*FreeElem)(void *));
int   ArrayFreeData(DATA *p, int esize, int block, 
		    void (*FreeElem)(void *));

int   SMultiInit(MULTI *ma, int esize, int ndim, int *block);
void *SMultiGet(MULTI *ma, int *k);
void *SMultiSet(MULTI *ma, int *k, void *d, 
		void (*InitData)(void *, int),
		void (*FreeElem)(void *));
int   SMultiFree(MULTI *ma, void (*FreeElem)(void *));
int   SMultiFreeDataOnly(ARRAY *a, int d, void (*FreeElem)(void *));
int   SMultiFreeData(MULTI *ma, void (*FreeElem)(void *));

/*
** the following set of funcitons are a different implementation
** for the MULTI array,
*/
int   NMultiInit(MULTI *ma, int esize, int ndim, int *block);
void *NMultiGet(MULTI *ma, int *k);
void *NMultiSet(MULTI *ma, int *k, void *d, 
		void (*InitData)(void *, int),
		void (*FreeElem)(void *));
int   NMultiFree(MULTI *ma, 
		 void (*FreeElem)(void *));
int   NMultiFreeDataOnly(ARRAY *a, void (*FreeElem)(void *));
int   NMultiFreeData(MULTI *ma, void (*FreeElem)(void *));

void  InitIntData(void *p, int n);
void  InitDoubleData(void *p, int n);
void  InitPointerData(void *p, int n);
void  InitArrayData(void *p, int n);
void  InitMDataData(void *p, int n);

#endif
