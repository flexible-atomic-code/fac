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

#ifndef _CONSTS_H_
#define _CONSTS_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "global.h"

/*
** VARIABLE:    VERSION, SUBVERSION, SUBSUBVERSION.
** TYPE:        macro constants.
** PURPOSE:     tracking the version.
** NOTE:        
*/
#define VERSION        1
#define SUBVERSION     1
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
** VARIABLE:    EPS1, ..., EPS30.
** TYPE:        macro constants.
** PURPOSE:     some small numbers.
** NOTE:        
*/
#define EPS30 1E-30
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
#define FOUR_PI     12.56637061436
#define ONETHIRD   0.33333333333
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
#define HARTREE_EV 27.211386018
#define RYDBERG_EV 13.605693009

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
#define RBOHR      0.52917721067

/*
** VARIABLE:    MBOHR
** TYPE:        macro constant.
** PURPOSE:     Bohr magenton in eV/Gauss.
** NOTE:        
*/
#define MBOHR      5.7883818012E-9

/*
** VARIABLE:    FINE_STRUCTURE_CONST, FINE_STRUCTURE_CONST2
** TYPE:        macro constants.
** PURPOSE:     fine structure constant and its square.
** NOTE:        
*/
#define FINE_STRUCTURE_CONST  7.2973525664E-3
#define FINE_STRUCTURE_CONST2 5.325135447834E-5
#define AMU  1822.888486

/* nucleus */
#define N_ELEMENTS 120
#define NISO 58
#define NFERMI 5001

#define NORBMAP0 500
#define NORBMAP1 1000
#define NORBMAP2 5000
#define KORBMAP 500

/* radial QK modes */
#define QK_DEFAULT    -1
#define QK_EXACT       0
#define QK_INTERPOLATE 1
#define QK_FIT         2
#define QK_CB          3
#define QK_DW          4
#define QK_BED         5

/* blocks for multi arrays */
#define MULTI_BLOCK2   512
#define MULTI_BLOCK3   64
#define MULTI_BLOCK4   32
#define MULTI_BLOCK5   4
#define MULTI_BLOCK6   10
#define MULTI_IDLEN 64
#define ARYCTH 0.05

/* orbital */
#define MAXRP      3000  /* maximum radial mesh */
#define DMAXRP     1200  /* default radial mesh points */
#define GRIDASYMP  36    /* no. points in one wavelength near infinity */
#define GRIDRATIO  1.1   /* ratio of successive mesh near origin */
#define GRIDRMIN   1E-6  /* starting point of the mesh is GRIDRMIN/Z */
#define GRIDRMINN0  1E-4  /* starting point relative to nucleus radius */
#define GRIDRMINN1  1E-2  /* starting point relative to nucleus radius */
#define GRIDQR 0.5 /* grid transform non-log term index */
#define ENERELERR  1E-6  /* relative energy error */
#define ENEABSERR  1E-4  /* absolute energy error */
#define ENERELERR1 1E-8

/* config */
#define MCHSHELL           2048
#define MAX_SPEC_SYMBOLS   21
#define LEVEL_NAME_LEN     1024
#define GROUP_NAME_LEN     64
#define MAX_GROUPS         1024
#define MAX_SYMMETRIES     256
#define CONFIGS_BLOCK      1024
#define STATES_BLOCK       2048

/* radial */
#define ORBITALS_BLOCK     1024
#define OPTSTABLE          0.5
#define OPTTOL             3.0
#define OPTNITER           512
#define OPTPRINT           0
#define POTMODE            0
#define POTHXS             1.0
#define POTIHX             -2.0
#define POTHX0             0.427
#define POTHX1             0.075
#define NKSEP              5
#define QEDSE              5
#define QEDMSE             41
#define QEDVP              3
#define QEDNMS             3
#define QEDSMS             3
#define QEDBREIT           10
#define QEDMBREIT          0
#define QEDNBREIT          5

/* structure */
#define MAX_HAMS           2000
#define LEVELS_BLOCK       1024
#define ANGZ_BLOCK         1024
#define ANGZxZ_BLOCK       8192
#define ANGZCUT            1E-5
#define MIXCUT             1E-5
#define MIXCUT2            1.0
#define NPRINCIPLE         2
#define MAXDN              3
#define MBCLOSE            8        
#define MAXLEVEB           1000000
#define DIAGMAXITER        0
#define DIAGMAXTOL         0.1
#define DGEEVMODE          0
#define PERTURBMAXITER     50
#define PERTURBEXPDIM      0.1
#define PERTURBEXPDIMZ      1e-4

/* transition */
#define G_COULOMB          1
#define G_BABUSHKIN        2
#define M_FR               0
#define M_NR               1
#define DGAUGE             G_BABUSHKIN
#define DMODE              M_NR
#define ERANK              4
#define MRANK              4
#define TRCUT0             1E-4
#define TRCUT              1E-4

/* recouple */
#define MAXRANK            20

/* coulomb */
#define NHYDROGEN          20
#define LHYDROGEN          10
#define NHYDROGENMAX       512
#define LHYDROGENMAX       20
#define DIPOLE_BLOCK       64
#define CBMULT             2
#define MAXNCB             ((CBMULT*(CBMULT+3))/2)
#define CBLMIN             15
#define CBLMAX             150

/* grid */
#define MAXNKL             50
#define MAXKL              512
#define MAXNUSR            30
#define MAXNE              20
#define MAXNTE             6
#define MAXNTHETA          30
#define MAXNPHI            60
#define TE_MIN_MAX         0.2

/* excitation */
#define NGOSK              256
#define EXCLQR             0
#define EXCLMAX            36
#define EXCLCB             36
#define EXCTOL             5E-2
#define XBORN              (-0.5)
#define XBORN1             (-1.0)
#define XBORN0             0.25
#define EBORN              100.0
#define MAXCECACHE         1000000

/* ionization */
#define IONMAXK            6
#define IONLQR             0
#define IONLMAX            36
#define IONLEJEC           4
#define IONLCB             36
#define IONTOL             5E-2

/* recombination */
#define RECNMAX            512
#define RECNSPEC           8
#define RECNFROZEN         8
#define RECLMAX            12
#define AICUT              0.0
#define MAXAICACHE         1000000

/* polarization */
#define MAXPOL             4 /* maximum multipol for polarization */

#endif
