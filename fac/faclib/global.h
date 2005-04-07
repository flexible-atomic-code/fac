#ifndef _GLOBAL_H_
#define _GLOBAL_H_ 1


/*************************************************************
  Header file defining some global constants, macros.

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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#ifdef PMALLOC_CHECK
#include "pmalloc.h"
#endif

#include "sysdef.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "consts.h"

/* choose multi-array implementation */
/*
#define MultiInit     NMultiInit 
#define MultiGet      NMultiGet 
#define MultiSet      NMultiSet 
#define MultiFree     NMultiFree 
#define MultiFreeData NMultiFreeData
*/

/*
** VARIABLE:    DEBUG_RECOUPLE, DEBUG_STRUCTURE, FAC_DEBUG
** TYPE:        macro constants
** PURPOSE:     debugging flags.
** NOTE:        
*/
#define DEBUG_RECOUPLE  10
#define DEBUG_STRUCTURE 20
#define FAC_DEBUG 0
#if FAC_DEBUG
/*
** VARIABLE:    debug_log
** TYPE:        global.
** PURPOSE:     file handler for the output of debug information.
** NOTE:        it is only defined when FAC_DEBUG is true.
*/
extern FILE *debug_log;
#endif

/*
** VARIABLE:    PERFORM_STATISTICS
** TYPE:        macro constant.
** PURPOSE:     indicate whether the profiling 
**              information should be compiled in.
** NOTE:        normally, should be commented out.
*/

/*#define PERFORM_STATISTICS 1*/
#ifdef PERFORM_STATISTICS
extern FILE *perform_log;
#endif

#endif

