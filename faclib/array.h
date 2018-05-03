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

#include <unistd.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <time.h>

#include "global.h"

#define USE_NMULTI 1

/* choose MULTI implementation */
#if USE_NMULTI == 1
#define MultiInit NMultiInit
#define MultiGet NMultiGet
#define MultiSet NMultiSet
#define MultiFreeData NMultiFreeData
#define MultiFree NMultiFree
#elif USE_NMULTI == 0
#define MultiInit SMultiInit
#define MultiGet SMultiGet
#define MultiSet SMultiSet
#define MultiFreeData SMultiFreeData
#define MultiFree SMultiFree
#else
#define MultiInit MMultiInit
#define MultiGet MMultiGet
#define MultiSet MMultiSet
#define MultiFreeData MMultiFreeData
#define MultiFree MMultiFree
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
  char id[MULTI_IDLEN];
  int esize;
  int block;
  int bsize;
  volatile int dim;
  DATA  *data;
  LOCK *lock;
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
  char id[MULTI_IDLEN];
  int numelem, maxelem;
  double totalsize, overheadsize, maxsize, cth;
  int clean_mode, clean_thread, clean_flag;
  unsigned long iset;
  unsigned short ndim, ndim1;
  unsigned short isize;
  unsigned short esize;
  unsigned short *block, *sidx, *ridx;
  int *iidx, *iblock;
  int isf, hsize, hmask, aidx;
  ARRAY *array;
  ARRAY *ia, *da;
  LOCK *lock, *clean_lock;
} MULTI;

typedef struct _IDXARY_ {
  int n, m, m0, m1;
  int *d, *i;
} IDXARY;

int IdxCmp(int *i0, int *i1, int n);
void InitMultiStats(void);
void ReportMultiStats(void);
void RemoveMultiLocks(void);
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
int   ArrayFreeLock(ARRAY *a, void (*FreeElem)(void *));
int   ArrayFreeData(DATA *p, int esize, int block, 
		    void (*FreeElem)(void *));

int   SMultiInit(MULTI *ma, int esize, int ndim, int *block, char *id);
void *SMultiGet(MULTI *ma, int *k, LOCK **lock);
void *SMultiSet(MULTI *ma, int *k, void *d, LOCK **lock,
		void (*InitData)(void *, int),
		void (*FreeElem)(void *));
int   SMultiFree(MULTI *ma, void (*FreeElem)(void *));
int   SMultiFreeDataOnly(ARRAY *a, int d, void (*FreeElem)(void *));
int   SMultiFreeData(MULTI *ma, void (*FreeElem)(void *));

/*
** the following set of funcitons are a different implementation
** for the MULTI array,
*/
int   NMultiInit(MULTI *ma, int esize, int ndim, int *block, char *id);
void *NMultiGet(MULTI *ma, int *k, LOCK **lock);
void *NMultiSet(MULTI *ma, int *k, void *d, LOCK **lock,
		void (*InitData)(void *, int),
		void (*FreeElem)(void *));
int   NMultiFree(MULTI *ma, 
		 void (*FreeElem)(void *));
int   MultiFreeLock(MULTI *ma, 
		    void (*FreeElem)(void *));
int   NMultiFreeDataOnly(ARRAY *a, void (*FreeElem)(void *));
int   NMultiFreeData(MULTI *ma, void (*FreeElem)(void *));

/*
** yet another implementation of MULTI array
*/
int   MMultiInit(MULTI *ma, int esize, int ndim, int *block, char *id);
void *MMultiGet(MULTI *ma, int *k, LOCK **lock);
void *MMultiSet(MULTI *ma, int *k, void *d, LOCK **lock,
		void (*InitData)(void *, int),
		void (*FreeElem)(void *));
int   MMultiFree(MULTI *ma, 
		 void (*FreeElem)(void *));
int   MMultiFreeData(MULTI *ma, void (*FreeElem)(void *));
void AddMultiSize(MULTI *ma, int size);
void LimitMultiSize(MULTI *ma, double d);

void  InitIntData(void *p, int n);
void  InitDoubleData(void *p, int n);
void  InitPointerData(void *p, int n);
void  InitArrayData(void *p, int n);
void  InitMDataData(void *p, int n);

void InitIdxAry(IDXARY *ia, int n, int *d);
void FreeIdxAry(IDXARY *ia, int md);
int IdxGet(IDXARY *ia, int d);
int          IBisect(int k, int n, int *a);
int          Bisect(void *p0, int n, int m, void *p,
		    int (*comp)(const void *, const void *));
void SetMultiCleanFlag(MULTI *ma);
double TotalSize(void);
double TotalArraySize(void);
#endif
