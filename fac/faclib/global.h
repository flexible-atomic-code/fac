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

/*
** VARIABLE:    VERSION, SUBVERSION, SUBSUBVERSION.
** TYPE:        macro constants.
** PURPOSE:     tracking the version.
** NOTE:        
*/
#define VERSION        1
#define SUBVERSION     0
#define SUBSUBVERSION  5
#define VersionGE(h, a, b, c)    (((h)->version >= (a)) &&\
                                  ((h)->sversion >= (b)) &&\
                                  ((h)->ssversion >= (c)))
#define VersionLE(h, a, b, c)    (((h)->version <= (a)) &&\
                                  ((h)->sversion <= (b)) &&\
                                  ((h)->ssversion <= (c)))

/* 
** MACRO:       IsOdd, IsEven
** PURPOSE:     determin if an integer is Odd or Even.
** INPUT:       {int x},
**              integer.
** RETURN:      {int},
**              0: false.
**              1: true.
** SIDE EFFECT: 
** NOTE:        
*/
#define IsOdd(x)  ((abs((x))&0x01)?1:0)
#define IsEven(x) ((abs((x))&0x01)?0:1)

/* 
** MACRO:       Max, Min
** PURPOSE:     the larger and lesser of two numbers.
** INPUT:       {generic a},
**              number participating the comparison.
**              {generic b},
**              number participating the comparison.
** RETURN:      {generic},
**              the larger or lesser of a and b.
** SIDE EFFECT: 
** NOTE:        
*/
#define Max(a, b) (((a)>(b))?(a):(b))
#define Min(a, b) (((a)<(b))?(a):(b))

#define IsNan(a)  (!((a)>0) && !((a)<=0)) 

/*
** VARIABLE:    EPS1, ..., EPS16.
** TYPE:        macro constants.
** PURPOSE:     some small numbers.
** NOTE:        
*/
#define EPS16 1E-16
#define EPS12 1E-12
#define EPS10 1E-10
#define EPS8  1E-08
#define EPS6  1E-06
#define EPS5  1E-05
#define EPS4  1E-04
#define EPS3  1E-03
#define EPS2  1E-02
#define EPS1  1E-01

/*
** VARIABLE:    PI, TWO_PI
** TYPE:        macro constants, double
** PURPOSE:     PI, and 2*PI
** NOTE:        
*/
#define PI         3.14159265359
#define TWO_PI     6.28318530718
/*
** VARIABLE:    SQRT2
** TYPE:        macro constant, double
** PURPOSE:     sqrt(2)
** NOTE:        
*/
#define SQRT2      1.41421356237

/*
** VARIABLE:    HARTREE_EV
** TYPE:        macro constant
** PURPOSE:     1 Hartree in eV.
** NOTE:        
*/
#define HARTREE_EV 27.2113962 

/*
** VARIABLE:    RATE_AU, RATE_AU10, RATE_AU12
** TYPE:        macro constants.
** PURPOSE:     atomic units of rate in 1/s, 10^10 1/s, 10^12 1/s
** NOTE:        
*/
#define RATE_AU    4.13413733E16
#define RATE_AU10  4.13413733E06
#define RATE_AU12  4.13413733E04

/*
** VARIABLE:    AREA_AU20
** TYPE:        macro constants.
** PURPOSE:     atomic units of area in 10^-20 cm2
** NOTE:        
*/
#define AREA_AU20  2.80028560859E3

/*
** VARIABLE:    RBOHR
** TYPE:        macro constant.
** PURPOSE:     Bohr radius in Angstrom.
** NOTE:        
*/
#define RBOHR      0.529177249

/*
** VARIABLE:    FINE_STRUCTURE_CONST, FINE_STRUCTURE_CONST2
** TYPE:        macro constants.
** PURPOSE:     fine structure constant and its square.
** NOTE:        
*/
#define FINE_STRUCTURE_CONST  7.29735308E-3
#define FINE_STRUCTURE_CONST2 5.32513620E-5
#define AMU  1836.153

/* nucleus */
#define N_ELEMENTS 109

/* radial QK modes */
#define QK_DEFAULT    -1
#define QK_EXACT       0
#define QK_INTERPOLATE 1
#define QK_FIT         2
#define QK_CB          3
#define QK_DW          4
#define QK_BED         5

/* blocks for multi arrays */
#define MULTI_BLOCK2   128
#define MULTI_BLOCK3   64
#define MULTI_BLOCK4   25
#define MULTI_BLOCK5   15
#define MULTI_BLOCK6   10

/* orbital */
#define MAXRP      3000 /* maximum radial mesh */
#define DMAXRP     1200 /* default radial mesh points */
#define GRIDASYMP  36   /* no. points in one wavelength near infinity */
#define GRIDRATIO  1.1  /* ratio of successive mesh near origin */
#define ENERELERR  1E-4 /* relative energy error */
#define ENEABSERR  1E-5 /* absolute energy error */

/* config */
#define NSPARTITION        100
#define MAX_SPEC_SYMBOLS   21
#define LEVEL_NAME_LEN     128
#define GROUP_NAME_LEN     64
#define MAX_GROUPS         4096
#define MAX_SYMMETRIES     256
#define CONFIGS_BLOCK      256
#define STATES_BLOCK       2048

/* radial */
#define ORBITALS_BLOCK     1024
#define OPTSTABLE          0.5
#define OPTTOL             1E-6
#define OPTNITER           128
#define OPTPRINT           0
#define QEDSE              5
#define QEDVP              2
#define QEDNMS             1
#define QEDSMS             1
#define QEDBREIT           5

/* structure */
#define LEVELS_BLOCK       1024
#define ANGZ_BLOCK         32
#define ANGZCUT            1E-6
#define MIXCUT             1E-5
#define NPRINCIPLE         2
#define MAXDN              3

/* transition */
#define G_COULOMB          1
#define G_BABUSHKIN        2
#define M_FR               0
#define M_NR               1
#define DGAUGE             G_BABUSHKIN
#define DMODE              M_NR
#define ERANK              4
#define MRANK              4
#define TRCUT              1E-4

/* recouple */
#define MAXRANK            16

/* coulomb */
#define NHYDROGEN          10
#define LHYDROGEN          7
#define NHYDROGENMAX       512
#define LHYDROGENMAX       20
#define DIPOLE_BLOCK       64

/* grid */
#define MAXNKL             50
#define MAXKL              512
#define MAXNUSR            30
#define MAXNE              20
#define MAXNTE             6
#define TE_MIN_MAX         (1.0/5.0)

/* excitation */
#define NGOSK              64
#define EXCLQR             0
#define EXCLMAX            36
#define EXCLCB             20
#define EXCTOL             5E-2
#define XBORN              (-0.35)

/* ionization */
#define IONMAXK            8
#define IONLQR             0
#define IONLMAX            36
#define IONLEJEC           8
#define IONLCB             20
#define IONTOL             5E-2

/* recombination */
#define RECNMAX            512
#define RECNSPEC           8
#define RECNFROZEN         8
#define RECLMAX            12
#define AICUT              0.0

/* polarization */
#define MAXPOL             4 /* maximum multipol for polarization */

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

