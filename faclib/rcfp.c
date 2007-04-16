#include "angular.h"
#include "rcfp.h"

static char *rcsid="$Id: rcfp.c,v 1.10 2005/08/20 21:26:19 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/*************************************************************
  Implementation of "rcfp".
  This module calculates the reduced coefficients of fractional
  parentage, the reduced matrix elements of creation, and
  annihilation operators between the single subshell states.

  It is mainly translated from the F90 package of 
  Gaigalas et al. CPC 139 (2001) 263.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/


/*
** VARIABLE:    terms_jj
** TYPE:        static array.
** PURPOSE:     tabulates the RCFP_TERMs for j <= 9/2
** NOTE:        the macro subshell_term() is to strip off the
**              identifier directly copied from CPC 139, 263.
*/
#define subshell_term(a, b, c, d, e) {a, b, c, d, e}
static RCFP_TERM terms_jj[63] = {      
  /* j = 1/2; these terms have indices terms_jj(1:2) */
  subshell_term(1, 0, 1, 1, 0), subshell_term(1, 1, 0, 0, 0),  
     
  /* j = 3/2; these terms have indices terms_jj(3:5) */
  subshell_term(3, 1, 1, 3, 0), subshell_term(3, 2, 0, 0, 0),  
  subshell_term(3, 0, 2, 4, 0),                                                
  /* j = 5/2; these terms have indices terms_jj(6:11) */
  subshell_term(5, 2, 1, 5, 0), subshell_term(5, 0, 3, 3, 0),   
  subshell_term(5, 0, 3, 9, 0), subshell_term(5, 3, 0, 0, 0),  
  subshell_term(5, 1, 2, 4, 0), subshell_term(5, 1, 2, 8, 0),  
      
  /* j = 7/2; these terms have indices terms_jj(12:25) */
  subshell_term(7, 1, 3, 3, 0), subshell_term(7, 1, 3, 5, 0),  
  subshell_term(7, 3, 1, 7, 0), subshell_term(7, 1, 3, 9, 0),  
  subshell_term(7, 1, 3,11, 0), subshell_term(7, 1, 3,15, 0),  
  subshell_term(7, 4, 0, 0, 0), subshell_term(7, 2, 2, 4, 0),  
  subshell_term(7, 0, 4, 4, 0), subshell_term(7, 2, 2, 8, 0),  
  subshell_term(7, 0, 4, 8, 0), subshell_term(7, 0, 4,10, 0),  
  subshell_term(7, 2, 2,12, 0), subshell_term(7, 0, 4,16, 0),  
      
  /* j = 9/2; these terms have indices terms_jj(26:63) */
  subshell_term(9, 0, 5, 1, 0), subshell_term(9, 2, 3, 3, 0),  
  subshell_term(9, 2, 3, 5, 0), subshell_term(9, 0, 5, 5, 0),  
  subshell_term(9, 2, 3, 7, 0), subshell_term(9, 0, 5, 7, 0),  
  subshell_term(9, 4, 1, 9, 0), subshell_term(9, 2, 3, 9, 0),  
  subshell_term(9, 0, 5, 9, 0), subshell_term(9, 2, 3,11, 0),  
  subshell_term(9, 0, 5,11, 0), subshell_term(9, 2, 3,13, 0),  
  subshell_term(9, 0, 5,13, 0), subshell_term(9, 2, 3,15, 0),  
  subshell_term(9, 0, 5,15, 0), subshell_term(9, 2, 3,17, 0),  
  subshell_term(9, 0, 5,17, 0), subshell_term(9, 0, 5,19, 0),  
  subshell_term(9, 2, 3,21, 0), subshell_term(9, 0, 5,25, 0),  
  subshell_term(9, 5, 0, 0, 0), subshell_term(9, 1, 4, 0, 0),  
  subshell_term(9, 3, 2, 4, 0), subshell_term(9, 1, 4, 4, 0),  
  subshell_term(9, 1, 4, 6, 0), subshell_term(9, 3, 2, 8, 0),  
  subshell_term(9, 1, 4, 8, 0), subshell_term(9, 1, 4, 8, 1),  
  subshell_term(9, 1, 4,10, 0), subshell_term(9, 3, 2,12, 0),  
  subshell_term(9, 1, 4,12, 0), subshell_term(9, 1, 4,12, 1),  
  subshell_term(9, 1, 4,14, 0), subshell_term(9, 3, 2,16, 0),  
  subshell_term(9, 1, 4,16, 0), subshell_term(9, 1, 4,18, 0),  
  subshell_term(9, 1, 4,20, 0), subshell_term(9, 1, 4,24, 0) };

/*
** VARIABLE:    rcfp_***, W_***
** TYPE:        static arrays
** PURPOSE:     tabulates the reduced matrix elements of 
**              A, AxW, WxA, WxW, for j <= 9/2.
** NOTE:        the macro reduced_coeff() is to strip off 
**              the identifier directly copied from CPC 139, 263.
*/
#define reduced_coeff(a, b, c) {a, b, c}
/* Set-up the reduced coefficients of fractional parentage (rcfp) */
/* j = 1/2  */
static REDUCED_COEFF rcfp_one_half[1][1] = {                            
  reduced_coeff(-1,      4,     1) };

/* j = 3/2 */            
static REDUCED_COEFF rcfp_three_half[2][1] = {                           
  reduced_coeff(-1,     12,     1), reduced_coeff(-1,     20,     1)};

/* j = 5/2 */         
static REDUCED_COEFF rcfp_five_half[3][3] = {
  {reduced_coeff(-1,     24,     1), reduced_coeff( 0,      1,     1), 
   reduced_coeff( 0,      1,     1)}, 
  {reduced_coeff(-1,     30,     1), reduced_coeff( 1,    120,     7), 
   reduced_coeff(-1,     90,     7)},
  {reduced_coeff(-1,     54,     1), reduced_coeff(-1,     48,     7), 
   reduced_coeff( 1,    330,     7)}};

/* j = 7/2	*/
static REDUCED_COEFF rcfp_seven_half[8][6] = {
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,     40,     1), reduced_coeff( 0,      1,     1), 
   reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  {reduced_coeff(-1,     54,     7), reduced_coeff(-1,     33,     1), 
  reduced_coeff(-1,     40,     1), reduced_coeff(-1,     65,     7), 
   reduced_coeff( 1,     30,     1), reduced_coeff( 0,      1,     1)},
  {reduced_coeff( 1,     88,     7), reduced_coeff(-1,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   1755,    77), 
   reduced_coeff(-1,     40,    11), reduced_coeff( 0,      1,     1)},
  {reduced_coeff( 1,    198,     7), reduced_coeff(-1,     36,    11), 
  reduced_coeff(-1,     72,     1), reduced_coeff( 1,   4500,    77), 
   reduced_coeff(-1,    234,    11), reduced_coeff(-1,    360,    11)}, 
  {reduced_coeff(-1,     78,    35), reduced_coeff( 1,    156,     5),
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,    108,    91), 
   reduced_coeff(-1,     30,     1), reduced_coeff( 1,     96,    13)},
  {reduced_coeff( 1,     66,     5), reduced_coeff( 1,     49,     5), 
   reduced_coeff( 0,      1,     1), reduced_coeff( 1,     27,     1),
   reduced_coeff( 1,    130,     7), reduced_coeff( 1,    136,     7)},
  {reduced_coeff( 0,      1,     1), reduced_coeff( 1,    195,    11), 
   reduced_coeff(-1,    104,     1), reduced_coeff(-1,    245,    11),
   reduced_coeff(-1,    624,    11), reduced_coeff( 1,   1224,    11)},
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   2720,   143), 
  reduced_coeff(-1,   2448,    77), reduced_coeff(-1,   7752,    91)} };

/* j = 9/2 */
static REDUCED_COEFF rcfp_nine_half_ket_46_48[3][20] = {
  /* rcfp for term_jj(46) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,     60,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(47) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,     12,     1), 
  reduced_coeff(-1,      8,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(48) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,     20,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,   1664,    33), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,     50,     1), reduced_coeff(-1,    130,    33), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,    272,    11), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,    560,    11), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}};

static REDUCED_COEFF rcfp_nine_half_ket_49_51[3][20] = {
  /* rcfp for term_jj(49) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,    130,     7), reduced_coeff(-1,     48,     7), 
  reduced_coeff( 1,     44,    21), reduced_coeff(-1,    340,    77), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,    140,    33), 
  reduced_coeff(-1,     70,    11), reduced_coeff( 1,   3808,   143), 
  reduced_coeff(-1,   2280,   143), reduced_coeff(-1,    110,    13), 
  reduced_coeff(-1,    918,   143), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(50) */
  {reduced_coeff( 0,      1,     1), reduced_coeff(-1,     63,     5), 
  reduced_coeff(-1,    312,    55), reduced_coeff(-1,      4,    11), 
  reduced_coeff( 1,      9,    11), reduced_coeff( 1,    255,    11), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   5040,   143), 
  reduced_coeff( 1,    840,   143), reduced_coeff(-1,   1428,   143), 
  reduced_coeff(-1,     95,   143), reduced_coeff( 1,    504,   715), 
  reduced_coeff(-1,   2856,   143), reduced_coeff( 1,  13566,   715), 
  reduced_coeff(-1,    850,   143), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(51) */
  {reduced_coeff( 0,      1,     1), reduced_coeff(-1,    384,    11), 
  reduced_coeff( 1,    156,    11), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,   1920,   143), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,     90,     1), reduced_coeff( 1,   7350,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   8160,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   1512,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   3648,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   9000,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}};


static REDUCED_COEFF rcfp_nine_half_ket_52_54[3][20] = {
  /* rcfp for term_jj(52) */
  {reduced_coeff(-1,   2184,   253), reduced_coeff(-1,     63,    23), 
  reduced_coeff(-1,  59904, 19481), reduced_coeff(-1, 302460, 19481), 
  reduced_coeff( 1,   5265,   161), reduced_coeff( 1, 848691,253253), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1, 145152, 36179), 
  reduced_coeff( 1, 217728, 36179), reduced_coeff( 1,1049580, 36179), 
  reduced_coeff( 1, 287337, 36179), reduced_coeff( 1,   5184,   299), 
  reduced_coeff( 1, 261120, 36179), reduced_coeff( 1, 691866, 36179), 
  reduced_coeff(-1,    750, 36179), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,  76608,  3289), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(53) */
  {reduced_coeff( 1,   1224,  1265), reduced_coeff( 1,   2652,   115), 
  reduced_coeff( 1, 188598, 13915), reduced_coeff(-1,  31824,  2783), 
  reduced_coeff( 1,    204,    23), reduced_coeff(-1,  12996, 13915), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,  25500,  2783), 
  reduced_coeff(-1,  38250,  2783), reduced_coeff(-1,    768,  2783), 
  reduced_coeff( 1,  81396, 13915), reduced_coeff( 1,   3213,   115), 
  reduced_coeff(-1,  60543,  2783), reduced_coeff(-1,3066144,236555), 
  reduced_coeff( 1, 727776, 47311), reduced_coeff( 1,    207,    17), 
  reduced_coeff(-1,  41553, 21505), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(54) */
  {reduced_coeff( 1,     52,     5), reduced_coeff(-1,     36,     5), 
  reduced_coeff( 1,     84,     5), reduced_coeff(-1,     56,    13), 
  reduced_coeff( 1,    126,    13), reduced_coeff(-1,   3570,   325), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,    360,    13), 
  reduced_coeff( 1,     60,    13), reduced_coeff( 1,    102,    13), 
  reduced_coeff(-1,   1064,    65), reduced_coeff(-1,   2142,    65), 
  reduced_coeff( 1,     42,    13), reduced_coeff( 1,    228,    65), 
  reduced_coeff( 1,    252,    13), reduced_coeff( 1,    342,    13), 
  reduced_coeff(-1,     66,    65), reduced_coeff( 1,    230,    13), 
   reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1)}};


static REDUCED_COEFF rcfp_nine_half_ket_55_57[3][20] = {
  /* rcfp for term_jj(55) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 1,    144,    11), 
  reduced_coeff( 1,    416,    11), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,     32,   165), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,    130,     1), reduced_coeff(-1,   1922,    33), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,    896,    55), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,    476,    11), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   1344,    11), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   2052,    55), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,    308,     5), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(56) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 1,    132,     5), 
  reduced_coeff(-1,  52728,  6655), reduced_coeff( 1,    196,  1331), 
  reduced_coeff(-1,     50,    11), reduced_coeff(-1,  24990,  1331), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,  21160,  1331), 
  reduced_coeff( 1,  31740,  1331), reduced_coeff( 1,  26250,  1331), 
  reduced_coeff( 1,  12920,  1331), reduced_coeff(-1,    357,    55), 
  reduced_coeff( 1,   9583,  1331), reduced_coeff(-1, 344988,  6655), 
  reduced_coeff( 1,   5700,  1331), reduced_coeff( 1,    171,    11), 
  reduced_coeff(-1,     15,   121), reduced_coeff(-1,   4830,   121), 
  reduced_coeff(-1,     84,    11), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(57) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1, 209950,  9317), reduced_coeff( 1,  77520,  9317), 
  reduced_coeff( 1,  12920,   231), reduced_coeff(-1,  25688,  9317), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   4522,  3993), 
  reduced_coeff( 1,   2261,  1331), reduced_coeff(-1,  48640, 22627), 
  reduced_coeff(-1, 285144,  9317), reduced_coeff( 1,    931,    44), 
  reduced_coeff( 1, 273885, 90508), reduced_coeff(-1, 112908,429913), 
  reduced_coeff(-1,2138580,158389), reduced_coeff( 1, 137781,  3740), 
  reduced_coeff( 1,6654375,156332), reduced_coeff(-1,  59616, 39083), 
  reduced_coeff( 1, 284089, 17765), reduced_coeff( 0,      1,     1)}};
 

static REDUCED_COEFF rcfp_nine_half_ket_58_60[3][20] = {
  /* rcfp for term_jj(58) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 1,   1530,    77), reduced_coeff( 1,  13056,  1001), 
  reduced_coeff( 1,  29376,  1001), reduced_coeff( 1,    720,  1001), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   1890,   143), 
  reduced_coeff(-1,    315,   143), reduced_coeff( 1, 101124,  2431), 
  reduced_coeff(-1,   4560,  1001), reduced_coeff(-1,  13965,   572), 
  reduced_coeff(-1,  35131,  9724), reduced_coeff( 1,  13500,  2431), 
  reduced_coeff(-1, 685900, 17017), reduced_coeff(-1,   1197,   884), 
  reduced_coeff(-1,  28875,   884), reduced_coeff(-1,   5060,   221), 
  reduced_coeff( 1,    759,    17), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(59) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 1,  22848,   715), reduced_coeff( 0,      1,     1), 
  reduced_coeff(-1,    170,     1), reduced_coeff( 1,    918,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,  32832,   715), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   9044,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,    576,    13), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   7524,    65), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 1,   1012,     5), reduced_coeff( 0,      1,     1)}, 
  /* rcfp for term_jj(60) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,   2128,   143), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   1938,   143), 
  reduced_coeff( 1,   2907,   143), reduced_coeff(-1,   4860,   143), 
  reduced_coeff(-1,  17136,  2717), reduced_coeff(-1,   1309,    52), 
  reduced_coeff(-1,   8505,   572), reduced_coeff( 1,  15876,   247), 
  reduced_coeff( 1,    420,    13), reduced_coeff( 1,   1287,    20), 
  reduced_coeff(-1,   6075,   988), reduced_coeff(-1,    132,    19), 
  reduced_coeff(-1,    253,    95), reduced_coeff(-1,    650,    19)}};
 
static REDUCED_COEFF rcfp_nine_half_ket_61_63[3][20] = {
  /* rcfp for term_jj(61) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 1,    570,    13), 
  reduced_coeff( 1,     95,    13), reduced_coeff( 1,   4104,   221), 
  reduced_coeff( 1,    504,    65), reduced_coeff(-1,   1463,   260), 
  reduced_coeff(-1,  39501,   884), reduced_coeff(-1,  60516,  1105), 
  reduced_coeff(-1,   1596,   221), reduced_coeff( 1,   3933,    68), 
  reduced_coeff( 1,    621,  3740), reduced_coeff(-1,    840,    17), 
  reduced_coeff(-1,    805,    17), reduced_coeff( 1,    390,    11)}, 
  /* rcfp for term_jj(62) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,   5796,   221), 
  reduced_coeff( 1,  17664,  1235), reduced_coeff(-1,   5313,    65), 
  reduced_coeff( 1,   1771,   221), reduced_coeff(-1,  16632,  1615), 
  reduced_coeff(-1,     88,    17), reduced_coeff( 1,    693,    17), 
  reduced_coeff(-1,  94269,  1615), reduced_coeff( 1, 192500,  7429), 
  reduced_coeff( 1,  30030,   323), reduced_coeff( 1,  24570,   437)}, 
  /* rcfp for term_jj(63) */
  {reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff( 0,      1,     1), 
  reduced_coeff( 0,      1,     1), reduced_coeff(-1,  15000,   323), 
  reduced_coeff( 1,    280,    17), reduced_coeff(-1,   1170,    17), 
  reduced_coeff(-1,  48750,  3553), reduced_coeff(-1,  15600,   437), 
  reduced_coeff( 1,   3510,    19), reduced_coeff(-1,  33930,   253)}};
				   

static REDUCED_COEFF W_00_one_half[2] = {
  reduced_coeff(-1,     2,     1), reduced_coeff(-1,       2,     1)};
 				   
   
static REDUCED_COEFF W_00_three_half[3] = {
  reduced_coeff(-1,    16,     1), reduced_coeff(-1,       6,     1), 
  reduced_coeff(-1,    10,     1)                                   };
    
   
static REDUCED_COEFF W_00_five_half[6] = {
  reduced_coeff(-1,    54,     1), reduced_coeff(-1,      12,     1), 
  reduced_coeff(-1,    30,     1), reduced_coeff(-1,      12,     1), 
  reduced_coeff(-1,    30,     1), reduced_coeff(-1,      54,     1)};
    
   
static REDUCED_COEFF W_00_seven_half[14] = {
  reduced_coeff(-1,    32,     1), reduced_coeff(-1,      48,     1), 
  reduced_coeff(-1,   128,     1), reduced_coeff(-1,      80,     1), 
  reduced_coeff(-1,    96,     1), reduced_coeff(-1,     128,     1), 
  reduced_coeff(-1,    20,     1), reduced_coeff(-1,      60,     1), 
  reduced_coeff(-1,    20,     1), reduced_coeff(-1,     108,     1), 
  reduced_coeff(-1,    36,     1), reduced_coeff(-1,      44,     1), 
  reduced_coeff(-1,   156,     1), reduced_coeff(-1,      68,     1)};
 
   
static REDUCED_COEFF W_10_one_half[2] = {
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,       6,     1)};
    
   
static REDUCED_COEFF W_10_three_half[3] = {
  reduced_coeff(-1,    12,     1), reduced_coeff(-1,      12,     1), 
  reduced_coeff( 0,     1,     1)                                   };
    
   
static REDUCED_COEFF W_10_five_half[6] = {
  reduced_coeff(-1,    48,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      20,     1), 
  reduced_coeff(-1,    10,     1), reduced_coeff(-1,      18,     1)};
    
   
static REDUCED_COEFF W_10_seven_half[14] = {
  reduced_coeff(-1,     6,     1), reduced_coeff(-1,       9,     1), 
  reduced_coeff(-1,   120,     1), reduced_coeff(-1,      15,     1), 
  reduced_coeff(-1,    18,     1), reduced_coeff(-1,      24,     1), 
  reduced_coeff(-1,    30,     1), reduced_coeff(-1,      30,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      54,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,    78,     1), reduced_coeff( 0,       1,     1)};
    
   
static REDUCED_COEFF W_01_one_half[2] = {
  reduced_coeff(-1,     6,     1), reduced_coeff( 0,       1,     1)};


static REDUCED_COEFF W_01_three_half[3] = {
  reduced_coeff(-1,    12,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,    12,     1)                                   };


static REDUCED_COEFF W_01_five_half[6] = {
  reduced_coeff(-1,   126,     7), reduced_coeff(-1,      12,     7), 
  reduced_coeff(-1,   198,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,    48,     7), reduced_coeff(-1,     288,     7)};
    
   
static REDUCED_COEFF W_01_seven_half[14] = {
  reduced_coeff(-1,    10,     7), reduced_coeff(-1,       5,     1), 
  reduced_coeff(-1,   168,     7), reduced_coeff(-1,     165,     7), 
  reduced_coeff(-1,   286,     7), reduced_coeff(-1,     680,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      30,     7), 
  reduced_coeff(-1,    10,     7), reduced_coeff(-1,     180,     7), 
  reduced_coeff(-1,    60,     7), reduced_coeff(-1,     110,     7), 
  reduced_coeff(-1,   546,     7), reduced_coeff(-1,     408,     7)};
    
  
static REDUCED_COEFF W_12_three_half_odd[1][1] = {
  reduced_coeff( 1,    60,     1)                   };
    
   
static REDUCED_COEFF W_12_three_half_even[2][2] = {
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,      30,     1)}, 
  {reduced_coeff(-1,    30,     1), reduced_coeff( 0,       1,     1)}};
    
   
static REDUCED_COEFF W_12_five_half_odd[3][3] = {
  {reduced_coeff( 1,   420,     7), reduced_coeff( 1,     360,     7), 
  reduced_coeff( 1,   270,     7)}, 
  {reduced_coeff( 1,     360,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1)}, 
  {reduced_coeff(-1,   270,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1)}                   };
    
   
static REDUCED_COEFF W_12_five_half_even[3][3] = {
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,    1960,    49), 
  reduced_coeff( 0,     1,     1)}, 
  {reduced_coeff(-1,    1960,    49), 
  reduced_coeff(-1,  1000,    49), reduced_coeff( 1,    2430,    49)}, 
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,    2430,    49), 
  reduced_coeff( 1,  1980,    49)}                   };
    
   
static REDUCED_COEFF W_12_seven_half_odd[21] = {
  reduced_coeff(-1,   252,    10), reduced_coeff( 1,    1056,    70), 
  reduced_coeff( 1,   144,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   507,    10), reduced_coeff(-1,      88,     1), 
  reduced_coeff(-1,   325,    14), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     200,     3), 
  reduced_coeff(-1,   520,    21), reduced_coeff( 1,      80,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,    6125,   462), 
  reduced_coeff(-1,   560,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   390,   539), reduced_coeff(-1,    3840,    49), 
  reduced_coeff( 1,  2040,    49)                                   };
    
   
static REDUCED_COEFF W_12_seven_half_even[36] = {
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      50,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,  1280,    49), reduced_coeff( 1,     990,    49), 
  reduced_coeff( 1,  2640,    49), reduced_coeff( 1,    1950,    49), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,   480,    49), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,     360,   539), 
  reduced_coeff( 1,  1872,    49), reduced_coeff( 1,      42,     1), 
  reduced_coeff( 1,   390,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,    48,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     234,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,    1040,    11), 
  reduced_coeff( 1,   340,     7), reduced_coeff( 0,       1,     1)};

 
static REDUCED_COEFF W_03_three_half_odd[1][1] = {
  reduced_coeff(-1,    28,     1)                   };
    
   
static REDUCED_COEFF W_03_three_half_even[2][2] = {
  {reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1)}, 
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,      28,     1)} };
    
   
static REDUCED_COEFF W_03_five_half_odd[3][3] = {
  {reduced_coeff(-1,   882,    21), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1)}, 
  {reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   384,    21), reduced_coeff(-1,     400,    21)}, 
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,     400,    21), 
  reduced_coeff( 1,   286,    21) }                  };
    
   
static REDUCED_COEFF W_03_five_half_even[3][3] = {
  {reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1)}, 
  {reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   162,     7), reduced_coeff(-1,     300,     7)}, 
  {reduced_coeff( 0,     1,     1), reduced_coeff(-1,     300,     7), 
  reduced_coeff( 1,    22,     7)}                   };
    
   
static REDUCED_COEFF W_03_seven_half_odd[21] = {
  reduced_coeff(-1,  1188,    70), reduced_coeff(-1,     196,    10), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,     234,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   189,   110), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,  1911,   242), reduced_coeff( 1,    1470,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      56,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,  394805, 22022), 
  reduced_coeff( 1,  5250,   121), reduced_coeff( 1,   53760,  1573), 
  reduced_coeff(-1,    78,   847), reduced_coeff( 1,   17408,   847), 
  reduced_coeff( 1, 12920,  1001)                                   };
    
   
static REDUCED_COEFF W_03_seven_half_even[36] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   110,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,   240,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,    1920,    77), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,      52,     7), 
  reduced_coeff(-1,   224,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,   32490,   847), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,  7644,   121), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,    12,    70), reduced_coeff( 1,     224,    10), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,    52,    70), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,  2040,    77), reduced_coeff(-1,     364,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,    1292,    77)};
    
  
static REDUCED_COEFF W_14_five_half_odd[3][3] = {
  {reduced_coeff( 1,  3780,    35), reduced_coeff(-1,     720,    35), 
  reduced_coeff(-1,  4950,    35)}, {reduced_coeff(-1,     720,    35), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1)}, 
  {reduced_coeff( 1,  4950,    35), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1)}                   };
    
   
static REDUCED_COEFF W_14_five_half_even[3][3] = {
  {reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  3528,    49)}, {reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  2430,    49), reduced_coeff( 1,    1980,    49)}, 
  {reduced_coeff(-1,  3528,    49), reduced_coeff( 1,    1980,    49), 
  reduced_coeff(-1,  7722,    49)}                   };
    
   
static REDUCED_COEFF W_14_seven_half_odd[21] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,      54,     7), 
  reduced_coeff(-1,   528,     7), reduced_coeff( 1,     546,    11), 
  reduced_coeff(-1,   378,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,  4335,   110), reduced_coeff(-1,      96,    11), 
  reduced_coeff( 1, 11271,  1694), reduced_coeff( 1,    3822,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     120,     1), 
  reduced_coeff( 1, 12000,    77), reduced_coeff(-1,     624,    11), 
  reduced_coeff(-1,   960,    11), reduced_coeff(-1,   30345,  3146), 
  reduced_coeff(-1,   210,   121), reduced_coeff(-1,  228480,  1573), 
  reduced_coeff( 1,580476,  5929), reduced_coeff( 1,  146880,  5929), 
  reduced_coeff(-1,627912,  7007)                                   };
       
   
static REDUCED_COEFF W_14_seven_half_even[36] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      90,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  2640,    49), reduced_coeff( 1,     480,    49), 
  reduced_coeff(-1,   360,   539), reduced_coeff( 1,   20592,   539), 
  reduced_coeff(-1,    42,     1), reduced_coeff( 1,     390,    11), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,468180,  5929), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,   21840,   847), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,  359424,  5929), 
  reduced_coeff(-1,  6750,   539), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,917280,  5929), reduced_coeff(-1,   10710,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,    36,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,     858,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,    5304,   121), 
  reduced_coeff(-1, 69768,   847), reduced_coeff( 0,       1,     1)};
      
   
static REDUCED_COEFF W_05_five_half_odd[3][3] = {
  {reduced_coeff(-1,  1386,    21), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1)}, {reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,     440,    21)}, 
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,     440,    21), 
  reduced_coeff(-1,  1430,    21) }                  };
       
   
static REDUCED_COEFF W_05_five_half_even[3][3] = {
  {reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1)}, {reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     330,     7)}, 
  {reduced_coeff( 0,     1,     1), reduced_coeff( 1,     330,     7), 
  reduced_coeff( 1,   572,     7) }                  };
       
   
static REDUCED_COEFF W_05_seven_half_odd[21] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     144,     7), 
  reduced_coeff( 1,    14,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   975,    22), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1, 50421,  2002), reduced_coeff( 1,     336,    11), 
  reduced_coeff(-1,  7000,   143), reduced_coeff(-1,      88,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,    6845, 26026), 
  reduced_coeff(-1,  3360,   143), reduced_coeff( 1,   14280,  1859), 
  reduced_coeff(-1,  1836,    77), reduced_coeff(-1,  103360,  1001), 
  reduced_coeff(-1, 28424,143143)                                   };
       
   
static REDUCED_COEFF W_05_seven_half_even[36] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   390,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,      70,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,      32,     7), 
  reduced_coeff( 1,    14,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,     576,    77), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,   210,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1, 63888,  1183), reduced_coeff(-1,     154,    13), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,    8568,   169), 
  reduced_coeff(-1,   176,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  1938,    91), reduced_coeff( 1,    1088,    11), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,   28424,  1183)};
  
   
static REDUCED_COEFF W_16_seven_half_odd[21] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     576,    11), 
  reduced_coeff( 1,   390,    77), reduced_coeff(-1,     144,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,     520,    11), 
  reduced_coeff(-1,    49,   121), reduced_coeff(-1,   12480,   121), 
  reduced_coeff(-1,   408,    11), reduced_coeff( 1,     520,     3), 
  reduced_coeff(-1,  1960,    33), reduced_coeff(-1,    1664,    11), 
  reduced_coeff( 1,  3264,    11), reduced_coeff( 1,  552250,  4719), 
  reduced_coeff(-1, 43520,   847), reduced_coeff(-1,   38760, 11011), 
  reduced_coeff(-1,  2652,   121), reduced_coeff( 1,   15504,   121), 
  reduced_coeff( 1, 38760,   143)                                   };
    
   
static REDUCED_COEFF W_16_seven_half_even[36] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1,   130,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   390,    11), reduced_coeff(-1,      48,     1), 
  reduced_coeff(-1,   234,     7), reduced_coeff( 1,    1040,    11), 
  reduced_coeff( 1,   340,     7), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  3120,   121), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,    1170,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,   18720,   121), 
  reduced_coeff(-1,    36,    11), reduced_coeff( 1,     858,     7), 
  reduced_coeff(-1,  5304,   121), reduced_coeff(-1,   69768,   847), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  1020,    11), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff(-1,   33592,   121), 
  reduced_coeff( 1, 31654,   121), reduced_coeff( 0,       1,     1)};
    
   
static REDUCED_COEFF W_07_seven_half_odd[21] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,   162,     7), reduced_coeff( 1,     272,     7), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1, 11025,  1573), reduced_coeff( 1,    4624,   121), 
  reduced_coeff(-1,  1632,   143), reduced_coeff(-1,     120,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1, 1224510, 20449), 
  reduced_coeff( 1,306000, 11011), reduced_coeff( 1,12558240,143143), 
  reduced_coeff(-1,  6460,   121), reduced_coeff( 1,   77520,  1573), 
  reduced_coeff(-1,297160,  1859)                                   };
    
   
static REDUCED_COEFF W_07_seven_half_even[36] = {
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,      60,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  1600,    77), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  2040,    77), reduced_coeff( 1,    4410,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1, 18360,   121), reduced_coeff( 0,       1,     1), 
  reduced_coeff(-1, 18816,   845), reduced_coeff(-1,   11016,   455), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,   11628,  1183), 
  reduced_coeff( 1,    34,     5), reduced_coeff( 0,       1,     1), 
  reduced_coeff( 1,  7752,   143), reduced_coeff( 1,    9690,   121), 
  reduced_coeff( 0,     1,     1), reduced_coeff( 1,  222870,  1859)};

/*
** VARIABLE:    rcfp_min_odd, rcfp_max_odd, *_even.
** TYPE:        static array
** PURPOSE:     the RCFP_TERM index limits for the 
**              dummy intermediate states inserted in between 
**              the operators.
** NOTE:        
*/
static int rcfp_min_odd[63] = {
  1, 0, 3, 2, 2, 8, 8, 8, 5, 5, 5,17,17,17,17,17,17,11,11,11,
  11,11,11,11,11,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,
  45,45,45,45,45,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,
  25,25,25 };

static int rcfp_max_odd[63] = {
  1, 0, 4, 2, 2,10,10,10, 7, 7, 7,24,24,24,24,24,24,16,16,16,
  16,16,16,16,16,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,
  62,62,62,62,62,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,
  44,44,44 };

static int rcfp_min_even[63] = {
  0, 1, 2, 3, 3, 5, 5, 5, 8, 8, 8,11,11,11,11,11,11,17,17,17,
  17,17,17,17,17,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,
  25,25,25,25,25,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,
  45,45,45 };

static int rcfp_max_even[63] = {
  0, 1, 2, 4, 4, 7, 7, 7,10,10,10,16,16,16,16,16,16,24,24,24,
  24,24,24,24,24,44,44,44,44,44,44,44,44,44,44,44,44,44,44,44,
  44,44,44,44,44,62,62,62,62,62,62,62,62,62,62,62,62,62,62,62,
  62,62,62 };


/* 
** FUNCTION:    ReducedCFP
** PURPOSE:     calculates the reduced coefficients 
**              of fractional parentage by looking up the table.
**              it's basically the reduced matrix elements of 
**              the creation or annihilation operator A in both
**              the angular and quasi-spin space.
** INPUT:       {int no_bra},
**              the RCFP_TEMR index of the bra state.
**              {int no_ket},
**              the RCFP_TERM index of the ket state.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        when the indexes are out of range, 
**              0.0 is returned, and warning is issued.
*/
double ReducedCFP(int no_bra, int no_ket) {
  double coeff;
  int no_a, no_b, phase, nom, denom;

  coeff = 0.0;

  if (rcfp_min_odd[no_bra] != rcfp_min_even[no_ket] ||
      !Triangle(terms_jj[no_ket].Q, 1, terms_jj[no_bra].Q) ||
      !Triangle(terms_jj[no_bra].subshellJ, 
		terms_jj[no_bra].j, 
		terms_jj[no_ket].subshellJ)) 
    return coeff;

  if (no_bra <= no_ket) {
    no_a = no_bra;  
    no_b = no_ket;
  } else {
    no_a = no_ket;   
    no_b = no_bra;
  }
  switch (terms_jj[no_bra].j) {
  case 1:
    phase = rcfp_one_half[no_b-1][no_a].phase;
    nom   = rcfp_one_half[no_b-1][no_a].nom;
    denom = rcfp_one_half[no_b-1][no_a].denom;
    break;
  case 3:
    phase = rcfp_three_half[no_b-3][no_a-2].phase;
    nom   = rcfp_three_half[no_b-3][no_a-2].nom;
    denom = rcfp_three_half[no_b-3][no_a-2].denom;
    break;
  case 5:
    phase = rcfp_five_half[no_b-8][no_a-5].phase;
    nom   = rcfp_five_half[no_b-8][no_a-5].nom;
    denom = rcfp_five_half[no_b-8][no_a-5].denom;
    break;
  case 7:
    phase = rcfp_seven_half[no_b-17][no_a-11].phase;
    nom   = rcfp_seven_half[no_b-17][no_a-11].nom;
    denom = rcfp_seven_half[no_b-17][no_a-11].denom;
    break;
  case 9:
    if (no_b > 44 && no_b < 48) {
      phase = rcfp_nine_half_ket_46_48[no_b-45][no_a-25].phase;
      nom   = rcfp_nine_half_ket_46_48[no_b-45][no_a-25].nom;
      denom = rcfp_nine_half_ket_46_48[no_b-45][no_a-25].denom;
    }
    if (no_b > 47 && no_b < 51) {
      phase = rcfp_nine_half_ket_49_51[no_b-48][no_a-25].phase;
      nom   = rcfp_nine_half_ket_49_51[no_b-48][no_a-25].nom;
      denom = rcfp_nine_half_ket_49_51[no_b-48][no_a-25].denom;
    }
    if (no_b > 50 && no_b < 54) {
      phase = rcfp_nine_half_ket_52_54[no_b-51][no_a-25].phase;
      nom   = rcfp_nine_half_ket_52_54[no_b-51][no_a-25].nom;
      denom = rcfp_nine_half_ket_52_54[no_b-51][no_a-25].denom;
    }
    if (no_b > 53 && no_b < 57) {
      phase = rcfp_nine_half_ket_55_57[no_b-54][no_a-25].phase;
      nom   = rcfp_nine_half_ket_55_57[no_b-54][no_a-25].nom;
      denom = rcfp_nine_half_ket_55_57[no_b-54][no_a-25].denom;
    }
    if (no_b > 56 && no_b < 60) {
      phase = rcfp_nine_half_ket_58_60[no_b-57][no_a-25].phase;
      nom   = rcfp_nine_half_ket_58_60[no_b-57][no_a-25].nom;
      denom = rcfp_nine_half_ket_58_60[no_b-57][no_a-25].denom;
    }
    if (no_b > 59 && no_b < 63) {
      phase = rcfp_nine_half_ket_61_63[no_b-60][no_a-25].phase;
      nom   = rcfp_nine_half_ket_61_63[no_b-60][no_a-25].nom;
      denom = rcfp_nine_half_ket_61_63[no_b-60][no_a-25].denom;
    }
    break;

  default:
    /*
    printf("improper input in ReducedCFP\n");
    */
    return 0.0;
  }

  if (phase != 0) {
    coeff = phase * sqrt(((double) nom)/denom);    
    if (no_bra >= no_ket  &&
	IsOdd((terms_jj[no_bra].Q - terms_jj[no_ket].Q +
	       terms_jj[no_bra].subshellJ - 
	       terms_jj[no_ket].subshellJ + 
	       terms_jj[no_bra].j - 1) / 2))
      coeff = - coeff;
  }
  return coeff;

}
       
/* 
** FUNCTION:    CompleteReducedW
** PURPOSE:     reduced matrix elements of W=AxA in both the
**              angular and quasi-spin space.
** INPUT:       {int no_bra},
**              bra state index.
**              {int no_ket},
**              ket state index.
**              {int k_q},
**              rank of the coupled operator in the quasi-spin space.
**              {int k_j},
**              rank of the coupled operator in the angular space.
** RETURN:      {double}
**              result
** SIDE EFFECT: 
** NOTE:        it checks if the angular momentum of the subshell 
**              is <= 9/2, in which case the coeff. are looked up 
**              in the table. otherwise, it is calculated using
**              decoupling formula of Racah.
*/
double CompleteReducedW(int no_bra, int no_ket, int k_q, int k_j) {
  double coeff;
  int jbra, jket, Jbra, Jket, Qbra, Qket, jrun, Jrun, Qrun;
  double w6j1, w6j2;
  int no_run;
  int phase;
  int kj2, kq2;

  coeff = 0.0;

  jbra = terms_jj[no_bra].j;
  if (jbra < 9) {
    coeff = CompleteReducedWFromTable(no_bra, no_ket, k_q, k_j);
    return coeff;
  }

  if (jbra == 9) {
    if (rcfp_min_even[no_bra] != rcfp_min_even[no_ket]) 
      return coeff;
    jket = terms_jj[no_ket].j;
    Jbra = terms_jj[no_bra].subshellJ;
    Jket = terms_jj[no_ket].subshellJ;
    Qbra = terms_jj[no_bra].Q;
    Qket = terms_jj[no_ket].Q;
    kq2 = k_q * 2;
    kj2 = k_j * 2;
    if (!Triangle(Qket, kq2, Qbra))
      return coeff;
    if (!Triangle(Jbra, kj2, Jket))
      return coeff;

    for (no_run = rcfp_min_odd[no_bra]; 
	 no_run <= rcfp_max_odd[no_bra];
	 no_run++) {
      Jrun = terms_jj[no_run].subshellJ;
      jrun = terms_jj[no_run].j;
      Qrun = terms_jj[no_run].Q;
      if ((w6j1 = W6j(jbra, jbra, kj2, Jket, Jbra, Jrun)) &&
	  (w6j2 = W6j(1, 1, kq2, Qket, Qbra, Qrun))) {
	coeff += w6j1 * w6j2 * (ReducedCFP(no_bra, no_run) *
				ReducedCFP(no_run, no_ket));
      }
    }

    if (fabs(coeff) > EPS30) {
      coeff *= sqrt((kq2 + 1.0) * (kj2 + 1.0));
      phase = Qbra + Jbra + Qket + Jket + kq2 + kj2;
      if (IsOdd(phase/2)) coeff = -coeff;
    }
    return coeff;
  }
  return coeff;
}


/* 
** FUNCTION:    CompleteReducedWFromTable
** PURPOSE:     reduced matrix elements of W=AxA in both the
**              angular and quasi-spin space by looking up the table
** INPUT:       {int no_bra},
**              bra state index.
**              {int no_ket},
**              ket state index.
**              {int k1},
**              rank of the coupled operator in the quasi-spin space.
**              {int k2},
**              rank of the coupled operator in the angular space.
** RETURN:      {double}
**              result     
** SIDE EFFECT: 
** NOTE:        the angular momentum of the subshell must <= 9/2
*/
double CompleteReducedWFromTable(int no_bra, int no_ket, int k1, int k2) {

  if (rcfp_min_even[no_bra] != rcfp_min_even[no_ket]) 
    return 0.0;
  if (!Triangle(terms_jj[no_ket].Q, 2*k1, terms_jj[no_bra].Q))
    return 0.0;
  if (!Triangle(terms_jj[no_bra].subshellJ, 2*k2, terms_jj[no_ket].subshellJ))
    return 0.0;
  
  if (k1 == 0 && k2 == 0) 
    return CompleteReducedWAll(W_00_one_half, W_00_three_half,
			       W_00_five_half, W_00_seven_half,
			       no_bra, no_ket);
  
  if (k1 == 1 && k2 == 0) 
    return CompleteReducedWAll(W_10_one_half, W_10_three_half, 
			       W_10_five_half,W_10_seven_half,
			       no_bra,no_ket);
  
  if (k1 == 0 && k2 == 1) 
    return CompleteReducedWAll(W_01_one_half, W_01_three_half, 
			       W_01_five_half,W_01_seven_half,
			       no_bra,no_ket);

  if (k1 == 1 && k2 == 2) {
    switch (terms_jj[no_bra].j) {
    case 3:
      return CompleteReducedW3(W_12_three_half_odd[0], 
			       W_12_three_half_even[0],
			       no_bra,no_ket);

    case 5:
      return CompleteReducedW5(W_12_five_half_odd[0],
			       W_12_five_half_even[0], 
			       no_bra,no_ket);

    case 7:
      return CompleteReducedW7(W_12_seven_half_odd,
			       W_12_seven_half_even,
			       no_bra,no_ket);
    default: 
      return 0.0;
    }
  }

  if (k1 == 0 && k2 == 3) {
    switch (terms_jj[no_bra].j) {
    case 3:
      return CompleteReducedW3(W_03_three_half_odd[0], 
			       W_03_three_half_even[0],
			       no_bra,no_ket);

    case 5:
      return CompleteReducedW5(W_03_five_half_odd[0],
			       W_03_five_half_even[0], 
			       no_bra,no_ket);

    case 7:
      return CompleteReducedW7(W_03_seven_half_odd,
			       W_03_seven_half_even,
			       no_bra,no_ket);
    default: 
      return 0.0;
    }
  }

  if (k1 == 1 && k2 == 4) {
    switch (terms_jj[no_bra].j) {
    case 5:
      return CompleteReducedW5(W_14_five_half_odd[0],
			       W_14_five_half_even[0], 
			       no_bra,no_ket);

    case 7:
      return CompleteReducedW7(W_14_seven_half_odd,
			       W_14_seven_half_even,
			       no_bra,no_ket);
    default: 
      return 0.0;
    }
  }

  if (k1 == 0 && k2 == 5) {
    switch (terms_jj[no_bra].j) {
    case 5:
      return CompleteReducedW5(W_05_five_half_odd[0],
			       W_05_five_half_even[0], 
			       no_bra,no_ket);

    case 7:
      return CompleteReducedW7(W_05_seven_half_odd,
			       W_05_seven_half_even,
			       no_bra,no_ket);
    default: 
      return 0.0;
    }
  }

  if (k1 == 1 && k2 == 6) {
    switch (terms_jj[no_bra].j) {
    case 7:
      return CompleteReducedW7(W_16_seven_half_odd,
			       W_16_seven_half_even,
			       no_bra,no_ket);
    default: 
      return 0.0;
    }
  }

  if (k1 == 0 && k2 == 7) {
    switch (terms_jj[no_bra].j) {
    case 7:
      return CompleteReducedW7(W_07_seven_half_odd,
			       W_07_seven_half_even,
			       no_bra,no_ket);
    default: 
      return 0.0;
    }
  }
  return 0.0;
}


/* 
** FUNCTION:    CompleteReducedW**
** PURPOSE:     perform actual table looking up for different cases.
** INPUT:       {REDUCED_COEFF *w}, 2 or 4,
**              pointer to the tables. 
**              {int no_bra},
**              bra state index.
**              {int no_ket},
**              ket state index.
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
double CompleteReducedWAll(REDUCED_COEFF *w1, REDUCED_COEFF *w3,
			   REDUCED_COEFF *w5, REDUCED_COEFF *w7,
			   int no_bra, int no_ket) {
  double coeff;
  int phase, nom, denom;
  
  coeff = 0.0;
  
  if (no_bra != no_ket) return coeff;

  switch (terms_jj[no_bra].j) {
  case 1:
    phase = w1[no_bra].phase;
    nom = w1[no_bra].nom;
    denom = w1[no_bra].denom;
    break;
    
  case 3:
    phase = w3[no_bra-2].phase;
    nom = w3[no_bra-2].nom;
    denom = w3[no_bra-2].denom;
    break;
    
  case 5:
    phase = w5[no_bra-5].phase;
    nom = w5[no_bra-5].nom;
    denom = w5[no_bra-5].denom;
    break;

  case 7:
    phase = w7[no_bra-11].phase;
    nom = w7[no_bra-11].nom;
    denom = w7[no_bra-11].denom;
    break;
  }
  
  if (phase) coeff = phase * sqrt(((double) nom)/denom);
  
  return coeff;
}

double CompleteReducedW3(REDUCED_COEFF *w3_o, REDUCED_COEFF *w3_e,
			 int no_bra, int no_ket) {
  double coeff;
  int phase, nom, denom;
  int i;

  coeff = 0.0;

  if (no_bra == 2) {
    phase = w3_o[0].phase;
    nom   = w3_o[0].nom;
    denom = w3_o[0].denom;
  } else if (no_bra > 2 && no_bra < 5) {    
    i = 2*(no_ket-3) + no_bra-3;
    phase = w3_e[i].phase;
    nom   = w3_e[i].nom;
    denom = w3_e[i].denom;
  } else {
    /*
    printf("improper input in CompleteReducedW3\n");
    */
    return 0.0;
  }

  if (phase) coeff = phase * sqrt(((double) nom)/denom);
  
  return coeff;
}


double CompleteReducedW5(REDUCED_COEFF *w5_o, REDUCED_COEFF *w5_e,
			 int no_bra, int no_ket) {
  double coeff;
  int phase, nom, denom;
  int i;

  coeff = 0.0;

  if (no_bra > 4 && no_bra < 8) {
    i = 3*(no_ket-5) + no_bra - 5;
    phase = w5_o[i].phase;
    nom   = w5_o[i].nom;
    denom = w5_o[i].denom;
  } else if (no_bra > 7 && no_bra < 11) {  
    i = 3*(no_ket-8) + no_bra - 8;  
    phase = w5_e[i].phase;
    nom   = w5_e[i].nom;
    denom = w5_e[i].denom;
  } else {
    /*
    printf("improper input in CompleteReducedW5\n");
    */
    return 0.0;
  }

  if (phase) coeff = phase * sqrt(((double) nom)/denom);
  
  return coeff;
}

double CompleteReducedW7(REDUCED_COEFF *w7_o, REDUCED_COEFF *w7_e,
			 int no_bra, int no_ket) {
  int phase, nom, denom, j;
  double coeff;

  int limits1[6] = {0, 5, 9, 12, 14, 15};
  int limits2[8] = {0, 7, 13, 18, 22, 25, 27, 28};
  int no_a, no_b;

  coeff = 0.0;
  if (no_bra < 11 && no_bra > 24) return coeff;
  if (no_ket < 11 && no_ket > 24) return coeff;

  if (no_bra > no_ket) {
    no_a = no_ket;
    no_b = no_bra;
  } else {
    no_a = no_bra;
    no_b = no_ket;
  }
  
  if (no_bra > 16) {
    no_a -= 17;
    no_b -= 17;
    j = limits2[no_a] + no_b;
    phase = w7_e[j].phase;
    nom = w7_e[j].nom;
    denom = w7_e[j].denom;
  } else {
    no_a -= 11;
    no_b -= 11;
    j = limits1[no_a] + no_b;
    phase = w7_o[j].phase;
    nom = w7_o[j].nom;
    denom = w7_o[j].denom;
  } 
  if (phase) {
    coeff = phase * sqrt(((double) nom)/denom); 
    if (no_bra > no_ket &&
	IsOdd((terms_jj[no_ket].subshellJ - terms_jj[no_bra].subshellJ +
	       terms_jj[no_ket].Q - terms_jj[no_bra].Q) / 2))
      coeff = -coeff;
  }
  return coeff;
}

/* 
** FUNCTION:    ReducedW
** PURPOSE:     reduced matrix elements of W=AxA only in
**              the angular space. 
** INPUT:       {RCFP_STATE *bra, *ket},
**              bra and ket states
**              {int k_j},
**              rank of the coupled operator in the angular space.
**              {int q_m1, q_m2},
**              the quasi-spin projection of the components.
**              i.e., whether they are creation or annihilaiton.
** RETURN:      {double}
**              result
** SIDE EFFECT: 
** NOTE:        if the input is inconsistent, 0.0 is returned,
**              and warning is issued.
*/
double ReducedW(RCFP_STATE *bra, RCFP_STATE *ket, 
		int k_j, int q_m1, int q_m2){
  double coeff, w6j1, a1, a2;
  RCFP_STATE run;
  int k_q, run_nu, Jrun, min_run, max_run;
  int run_i, Jbra, jbra, Jket, jket, Qbra, Qket, kj2, kq2;

  coeff = 0.0;
  kj2 = 2 * k_j;
  if (bra->state < 63 && ket->state < 63) {
    Jbra = terms_jj[bra->state].subshellJ;
    Jket = terms_jj[ket->state].subshellJ;
    Qbra = terms_jj[bra->state].Q;
    Qket = terms_jj[ket->state].Q;
    jbra = terms_jj[bra->state].j;
    jket = terms_jj[ket->state].j;
    if (rcfp_min_even[bra->state] != rcfp_min_even[ket->state] ||
	!Triangle(Jbra, kj2, Jket) ||
	QSpaceDelta(bra) == 0 ||
	QSpaceDelta(ket) == 0 ||
	bra->nq - ket->nq - q_m1 - q_m2)
      return coeff;
    
    if (q_m1 == q_m2) {
      if (IsOdd(k_j)) return coeff;
      if (!Triangle(Qket, 2, Qbra)) return coeff;
      coeff = ClebschGordanQuasispin(Qket, ket->subshellMQ, 2,
				     q_m1+q_m2, Qbra, bra->subshellMQ);
      if (fabs(coeff) < EPS30) return coeff;
      coeff = coeff * CompleteReducedW(bra->state, ket->state, 1, k_j);
      if (fabs(coeff) < EPS30) return coeff;
      coeff = coeff / sqrt(Qbra + 1.0);
    } else {
      if (k_j == 0) {
	if (bra->state != ket->state) return coeff;
	if (q_m1 == 1) coeff = -bra->nq;
	else coeff = jbra + 1.0 - bra->nq;
	coeff *= sqrt((Jbra + 1.0) / (jbra + 1.0));
      } else {
	if (IsEven(k_j)) k_q = 1;
	else k_q = 0;
	kq2 = 2*k_q;
	if (!Triangle(Qket, kq2, Qbra)) return coeff;
	coeff = ClebschGordanQuasispin(Qket, ket->subshellMQ, kq2,
				       q_m1 + q_m2, Qbra, bra->subshellMQ);
	if (fabs(coeff) < EPS30) return coeff;
	coeff *= CompleteReducedW(bra->state, ket->state, k_q, k_j);
	if (fabs(coeff) < EPS30) return coeff;
	coeff /= sqrt(2.0*(Qbra + 1.0));
	if (q_m1 == -1 && IsOdd(k_j)) coeff = -coeff;
      }
    }
  } else if (bra->state > 63 && ket->state > 63) {
    UnpackRCFPState(bra->state, &jbra, &run_nu, &Jbra);
    UnpackRCFPState(ket->state, &jket, &run_nu, &Jket);
    if (ket->nq > 2 || bra->nq > 2) {
      /*
      printf("1. improper input in ReducedW, \n");
      */
      return 0.0;
    }
    if (q_m1 + q_m2 + ket->nq != bra->nq) return 0.0;

    run.nq = ket->nq + q_m2;
    switch (run.nq) {
    case 0:
      min_run = 0;
      max_run = 0;
      break;
    case 1:
      min_run = jbra;
      max_run = jbra;
      break;
    case 2:
      min_run = 0;
      max_run = 2*jbra - 2;
      break;
    default:
      return 0.0;
    }
    for (run_i = min_run; run_i <= max_run; run_i += 4) {
      Jrun = run_i;
      if ((w6j1 = W6j(jbra, jbra, kj2, Jket, Jbra, Jrun))) {
	run_nu = run.nq;
	if (run.nq == 2 && Jrun == 0) run_nu = 0;
	run.subshellMQ = run.nq - (jbra+1)/2;
	run.state = PackRCFPState(jbra, run_nu, Jrun);
	a1 = ReducedA(bra, &run, q_m1);
	a2 = ReducedA(&run, ket, q_m2);
	coeff += (a1*a2) * w6j1;
      }
    }

    coeff *= sqrt(kj2 + 1.0);
    if (IsOdd((Jbra + Jket + kj2)/2)) coeff = -coeff;
  } else {
    /*
    printf("2. improper input in ReducedW\n");
    */
    return 0.0;
  }
  return coeff;
}

/* 
** FUNCTION:    ReducedWxW0
** PURPOSE:     reduced matrix element of WxW with final rank 0
** INPUT:       {RCFP_STATE *bra, *ket},
**              pointer to the bra and ket states.
**              {int k_j},
**              rank of the component W.
**              {int q_m1, q_m2, q_m3, q_m4},
**              quasi-spin projection of all components of W.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        if the input is inconsistent, 0.0 is returned,
**              and warning is issued.
*/
double ReducedWxW0(RCFP_STATE *bra, RCFP_STATE *ket,
		   int k_j, int q_m1, int q_m2, int q_m3, int q_m4) {
  int no_run;
  RCFP_STATE run;
  double coeff, coeff1;
  int run_nu, Jrun, jrun, Qrun, min_run, max_run;
  int run_i, Jbra, jbra, Jket, jket, Qbra, Qket, kj2;
  double w6j1;

  coeff = 0.0;
  
  kj2 = 2*k_j;
  if (bra->state < 63 && ket->state < 63) {
    Jbra = terms_jj[bra->state].subshellJ;
    Jket = terms_jj[ket->state].subshellJ;
    Qbra = terms_jj[bra->state].Q;
    Qket = terms_jj[ket->state].Q;
    jbra = terms_jj[bra->state].j;
    jket = terms_jj[ket->state].j;
    if (QSpaceDelta(bra) == 0) return coeff;
    if (QSpaceDelta(ket) == 0) return coeff;
    if (Jbra != Jket) return coeff;
    if (rcfp_min_even[bra->state] != rcfp_min_even[ket->state]) return coeff;
    if (rcfp_max_even[bra->state] != rcfp_max_even[ket->state]) return coeff;
    if (bra->nq - ket->nq - q_m1 - q_m2 - q_m3 - q_m4) return coeff;

    min_run = rcfp_min_even[bra->state];
    max_run = rcfp_max_even[bra->state];
    for (no_run = min_run; no_run <= max_run; no_run++) {
      jrun = terms_jj[no_run].j;
      Qrun = terms_jj[no_run].Q;
      Jrun = terms_jj[no_run].subshellJ;
      run.state = no_run;
      run.n = bra->n;
      run.nq = ket->nq + q_m3 + q_m4;
      run.subshellMQ = run.nq - (jrun +1)/2;
      if (Qrun >= abs(run.subshellMQ)) {
	if (W6jTriangle(kj2, kj2, 0, Jket, Jbra, Jrun)) {
	  coeff1 = (ReducedW(bra, &run, k_j, q_m1, q_m2) *
		    ReducedW(&run, ket, k_j, q_m3, q_m4));
	  if (IsOdd((kj2 - Jbra + Jrun)/2)) coeff1 = -coeff1;
	  coeff += coeff1;
	}
      }
    }

    coeff /= sqrt((kj2 + 1.0) * (Jbra + 1.0));
  } else if (bra->state > 63 && ket->state > 63) {
    if (bra->state != ket->state) return coeff;
    UnpackRCFPState(bra->state, &jbra, &run_nu, &Jbra);
    UnpackRCFPState(ket->state, &jket, &run_nu, &Jket);
    run.nq = ket->nq + q_m3 + q_m4;
    switch (run.nq) {
    case 0:
      min_run = 0;
      max_run = 0;
      break;
    case 1:
      min_run = jbra;
      max_run = jbra;
      break;
    case 2:
      min_run = 0;
      max_run = 2*jbra - 2;
      break;
    default:
      /*
      printf("improper input in ReducedWxW0\n");
      */
      return 0.0;
    }

    for (run_i = min_run; run_i <= max_run; run_i += 4) {
      Jrun = run_i;
      if ((w6j1 = W6j(kj2, kj2, 0, Jket, Jbra, Jrun))) {
	run_nu = run.nq;
	if (run.nq == 2 && Jrun == 0) run_nu = 0;
	run.subshellMQ = run.nq - (jbra + 1)/2;
	run.state = PackRCFPState(jbra, run_nu, Jrun);
	coeff += (ReducedW(bra, &run, k_j, q_m1, q_m2) *
		  ReducedW(&run, ket, k_j, q_m3, q_m4) * w6j1);
      }
    }

    if (IsOdd((Jbra + Jket)/2)) coeff = -coeff;
  }
  return coeff;
}
 
/* 
** FUNCTION:    ReducedAxW
** PURPOSE:     reduced matrix element of AxW.
** INPUT:       {RCFP_STATE *bra, *ket},
**              pointers to the bra and ket states.
**              {int k_j1},
**              rank of W.
**              {int kk_j2},
**              rank of the final coupled operator.
**              {int q_m1, q_m2, q_m3},
**              quasi-spin projection of the 3 compenents.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        if the input is inconsistent, 0.0 is returned,
**              and warning is issued.
*/    
double ReducedAxW(RCFP_STATE *bra, RCFP_STATE *ket,
		  int k_j1, int kk_j2, int q_m1, int q_m2, int q_m3) {
  double coeff, coeff1;
  int no_run;
  RCFP_STATE run;
  int run_nu, Jrun, jrun, Qrun, min_run, max_run;
  int run_i, Jbra, jbra, Jket, jket, Qbra, Qket, kj2;
  double w6j1;

  coeff = 0.0;
  
  kj2 = 2 * k_j1;
  if (bra->state < 63 && ket->state < 63) {    
    Jbra = terms_jj[bra->state].subshellJ;
    Jket = terms_jj[ket->state].subshellJ;
    Qbra = terms_jj[bra->state].Q;
    Qket = terms_jj[ket->state].Q;
    jbra = terms_jj[bra->state].j;
    jket = terms_jj[ket->state].j;
    if (!Triangle(Jbra, kk_j2, Jket)) return coeff;
    if (QSpaceDelta(bra) == 0) return coeff;
    if (QSpaceDelta(ket) == 0) return coeff;

    if (rcfp_min_odd[bra->state] != rcfp_min_even[ket->state]) return coeff;
    if (rcfp_max_odd[bra->state] != rcfp_max_even[ket->state]) return coeff;
    if (bra->nq - ket->nq - q_m1 - q_m2 - q_m3) return coeff;

    min_run = rcfp_min_odd[bra->state];
    max_run = rcfp_max_odd[bra->state];
    for (no_run = min_run; no_run <= max_run; no_run++) {
      jrun = terms_jj[no_run].j;
      Qrun = terms_jj[no_run].Q;
      Jrun = terms_jj[no_run].subshellJ;
      run.state = no_run;
      run.n = bra->n;
      run.nq = ket->nq + q_m2 + q_m3;
      run.subshellMQ = run.nq - (jrun +1)/2;

      if (Qrun >= abs(run.subshellMQ)) {
	if ((w6j1 = W6j(jbra, kj2, kk_j2, Jket, Jbra, Jrun))) {
	  coeff1 = ClebschGordanQuasispin(Qrun, run.subshellMQ, 1, 
					  q_m1, Qbra, bra->subshellMQ);
	  if (fabs(coeff1) > EPS30) {
	    coeff1 *= ReducedCFP(bra->state, run.state) *
	      ReducedW(&run, ket, k_j1, q_m2, q_m3) * w6j1;
	    coeff += coeff1/sqrt(Qbra + 1.0);
	  }
	}
      }
    }
    coeff *= sqrt(kk_j2 + 1.0);
    if (IsEven((kk_j2 + Jbra + Jket)/2)) coeff = -coeff;
  } else if (bra->state > 63 && ket->state > 63) {
    UnpackRCFPState(bra->state, &jbra, &run_nu, &Jbra);
    UnpackRCFPState(ket->state, &jket, &run_nu, &Jket);
    run.nq = ket->nq + q_m3 + q_m2;
    switch (run.nq) {
    case 0:
      min_run = 0;
      max_run = 0;
      break;
    case 1:
      min_run = jbra;
      max_run = jbra;
      break;
    case 2:
      min_run = 0;
      max_run = 2*jbra - 2;
      break;
    default:
      /*
      printf("improper input in ReducedAxW\n");
      */
      return 0.0;
    }

    for (run_i = min_run; run_i <= max_run; run_i += 4) {
      Jrun = run_i;
      if ((w6j1 = W6j(jbra, kj2, kk_j2, Jket, Jbra, Jrun))) {
	run_nu = run.nq;
	if (run.nq == 2 && Jrun == 0) run_nu = 0;
	run.subshellMQ = run.nq - (jbra + 1)/2;
	run.state = PackRCFPState(jbra, run_nu, Jrun);
	coeff += (ReducedA(bra, &run, q_m1) *
		  ReducedW(&run, ket, k_j1, q_m2, q_m3)) * w6j1;
      }
    }
    coeff *= sqrt(kk_j2 + 1.0);
    if (IsOdd((Jbra + kk_j2 + Jket)/2)) coeff = -coeff;
  }
  return coeff;
    
}

/* 
** FUNCTION:    ReducedWxA
** PURPOSE:     reduced matrix element of WxA.
** INPUT:       {RCFP_STATE *bra, *ket},
**              pointers to the bra and ket states.
**              {int k_j1},
**              rank of W.
**              {int kk_j2},
**              rank of the final coupled operator.
**              {int q_m1, q_m2, q_m3},
**              quasi-spin projection of the 3 compenents.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        if the input is inconsistent, 0.0 is returned,
**              and warning is issued.
*/    
double ReducedWxA(RCFP_STATE *bra, RCFP_STATE *ket,
		  int k_j1, int kk_j2, int q_m1, int q_m2, int q_m3) {
  double coeff, coeff1;
  int no_run;
  RCFP_STATE run;
  int run_nu, Jrun, jrun, Qrun, min_run, max_run;
  int run_i, Jbra, jbra, Jket, jket, Qbra, Qket, kj2;
  double w6j1;

  coeff = 0.0;
  
  kj2 = 2 * k_j1;
  if (bra->state < 63 && ket->state < 63) {    
    Jbra = terms_jj[bra->state].subshellJ;
    Jket = terms_jj[ket->state].subshellJ;
    Qbra = terms_jj[bra->state].Q;
    Qket = terms_jj[ket->state].Q;
    jbra = terms_jj[bra->state].j;
    jket = terms_jj[ket->state].j;
    if (!Triangle(Jbra, kk_j2, Jket)) return coeff;
    if (QSpaceDelta(bra) == 0) return coeff;
    if (QSpaceDelta(ket) == 0) return coeff;
    if (rcfp_min_even[bra->state] != rcfp_min_odd[ket->state]) return coeff;
    if (rcfp_max_even[bra->state] != rcfp_max_odd[ket->state]) return coeff;
    if (bra->nq - ket->nq - q_m1 - q_m2 - q_m3) return coeff;

    min_run = rcfp_min_even[bra->state];
    max_run = rcfp_max_even[bra->state];
    for (no_run = min_run; no_run <= max_run; no_run++) {
      jrun = terms_jj[no_run].j;
      Qrun = terms_jj[no_run].Q;
      Jrun = terms_jj[no_run].subshellJ;
      run.state = no_run;
      run.n = bra->n;
      run.nq = ket->nq + q_m3;
      run.subshellMQ = run.nq - (jrun +1)/2;
      if (Qrun >= abs(run.subshellMQ)) {
	if ((w6j1 = W6j(kj2, jbra, kk_j2, Jket, Jbra, Jrun))) {
	  coeff1 = ClebschGordanQuasispin(Qket, ket->subshellMQ, 1, 
					  q_m3, Qrun, run.subshellMQ);
	  if (fabs(coeff1) > EPS30) {
	    coeff1 *= (ReducedCFP(run.state, ket->state) *
		       ReducedW(bra, &run, k_j1, q_m1, q_m2)) * w6j1;
	    coeff += coeff1/sqrt(Qrun + 1.0);
	  }
	}
      }
    }
    coeff *= sqrt(kk_j2 + 1.0);
    if (IsEven((kk_j2 + Jbra + Jket)/2)) coeff = -coeff;
  } else if (bra->state > 63 && ket->state > 63) {
    UnpackRCFPState(bra->state, &jbra, &run_nu, &Jbra);
    UnpackRCFPState(ket->state, &jket, &run_nu, &Jket);
    run.nq = ket->nq + q_m3;
    switch (run.nq) {
    case 0:
      min_run = 0;
      max_run = 0;
      break;
    case 1:
      min_run = jbra;
      max_run = jbra;
      break;
    case 2:
      min_run = 0;
      max_run = 2*jbra - 2;
      break;
    default:
      /*
      printf("improper input in ReducedWxA\n");
      */
      return 0.0;
    }

    for (run_i = min_run; run_i <= max_run; run_i += 4) {
      Jrun = run_i;
      if ((w6j1 = W6j(kj2, jbra, kk_j2, Jket, Jbra, Jrun))) {
	run_nu = run.nq;
	if (run.nq == 2 && Jrun == 0) run_nu = 0;
	run.subshellMQ = run.nq - (jbra + 1)/2;
	run.state = PackRCFPState(jbra, run_nu, Jrun);
	coeff += (ReducedA(&run, ket, q_m1) *
		  ReducedW(bra, &run, k_j1, q_m2, q_m3)) * w6j1;
      }
    }
    coeff *= sqrt(kk_j2 + 1.0);
    if (IsOdd((Jbra + kk_j2 + Jket)/2)) coeff = -coeff;
  }
  return coeff;
    
}

/* 
** FUNCTION:    ReducedA
** PURPOSE:     reduced matrix element of A. 
**              i.e., the regular coeff. of fractional parentage.
** INPUT:       {RCFP_STATE *bra, *ket},
**              pointers to the bra and ket states.
**              {int q_m},
**              quasi-spin projection of the operator,
**              which indicates whether it is creation or annihilaiton.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        if the input is inconsistent, 0.0 is returned,
**              and warning is issued.
*/    
double ReducedA(RCFP_STATE *bra, RCFP_STATE *ket, int q_m) {
  double coeff;
  int bra_nu, ket_nu;  
  int Jbra, jbra, Jket, jket, Qbra, Qket;

  coeff = 0.0;

  if (bra->state < 63 && ket->state < 63) {
    if (rcfp_min_odd[bra->state] != rcfp_min_even[ket->state]) return coeff;
    if (bra->nq - ket->nq - q_m != 0) return coeff;
    if (QSpaceDelta(bra) == 0) return coeff;
    if (QSpaceDelta(ket) == 0) return coeff;

    Jbra = terms_jj[bra->state].subshellJ;
    Jket = terms_jj[ket->state].subshellJ;
    Qbra = terms_jj[bra->state].Q;
    Qket = terms_jj[ket->state].Q;
    jbra = terms_jj[bra->state].j;
    jket = terms_jj[ket->state].j;
    
    if (!Triangle(Qket, 1, Qbra)) return coeff;
    if (!Triangle(Jbra, jbra, Jket)) return coeff;

    coeff = -ClebschGordanQuasispin(Qket, ket->subshellMQ, 1,
				     q_m, Qbra, bra->subshellMQ); 
    if (fabs(coeff) < EPS30) return coeff;
    coeff *= ReducedCFP(bra->state, ket->state);
    coeff /= sqrt(Qbra + 1.0);
  } else if (bra->state > 63 && ket->state > 63) {
    UnpackRCFPState(bra->state, &jbra, &bra_nu, &Jbra);
    if (q_m == -1) {
      UnpackRCFPState(ket->state, &jket, &ket_nu, &Jket);
      coeff = sqrt(ket->nq * (Jket + 1.0));
      if (IsOdd((Jbra - Jket -jbra)/2 + ket->nq)) coeff = -coeff;
    } else if (q_m == 1) {
      coeff = sqrt(bra->nq * (Jbra + 1.0));
      if (IsOdd(bra->nq)) coeff = -coeff;
    }
  } else {
    /*
    printf("improper input in ReducedA\n");
    */
    return 0.0;
  }

  return coeff;
}

/* 
** FUNCTION:    QSpaceDelta
** PURPOSE:     determine if the quasi-spin values are consistent.
** INPUT:       {RCFP_STATE *s},
**              pointer to the state.
** RETURN:      {int},
**              0: inconsistent.
**              1: consistent.
** SIDE EFFECT: 
** NOTE:        
*/
int QSpaceDelta(RCFP_STATE *s) {
  if (terms_jj[s->state].Q < abs(s->subshellMQ)) return 0;
  if (IsOdd(terms_jj[s->state].Q + s->subshellMQ)) return 0;
  return 1;
}

/* 
** FUNCTION:    ClebschGordanQuasispin
** PURPOSE:     the C-G coeff. occuring in quasi-spin space are
**              usually of small angular momentum, this used more
**              specialized formulae whenever possible.
** INPUT:       {int ja, ma},
**              angular momentum and its projection.
**              {int jb, mb},
**              angular momentum and its projection.
**              {int jab, mab},
**              coupled angular momentum and its projection.
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/
double ClebschGordanQuasispin(int ja, int ma, int jb, int mb, 
			      int jab, int mab) {
  double CG;

  CG = 0.0;
  if (ma + mb != mab) return CG;
  if (!Triangle(ja, jb, jab)) return CG;
  if (abs(ma) > ja || abs(mb) > jb || abs(mab) > jab) return CG;
  if (IsOdd(abs(ma) + ja) || IsOdd(abs(mb) + jb) ||
      IsOdd(abs(mab) + jab))
    return CG;

  switch (jb) {
  case 0:
    if (ja != jab) return CG;
    if (ma != mab) return CG;
    CG = 1.0;
    break;
  case 1:
    if (ja + 1 == jab) 
      CG = sqrt(0.5 * (jab + mb*mab) / jab);
    else if (ja - 1 == jab) 
      CG = -mb * sqrt(0.5 * (jab - mb*mab + 2.0) / (jab + 2.0));
    break;
  case 2:
    if (mb == 0) {
      if (ja + 2 == jab) 
	CG = sqrt((0.5 * (jab + mab) * (jab - mab)) / ((jab - 1.0) * jab));
      else if (ja == jab) 
	CG = (0.5 * mab) / (0.5 * sqrt(jab * (jab + 2.0)));
      else if (ja - 2 == jab) 
	CG = -sqrt(0.5 * ((jab + mab + 2.0) * (jab - mab + 2.0))/
		   ((jab + 2.0) * (jab + 3.0)));
    } else if (mb == 2 || mb == -2) {
      if (ja + 2 == jab) 
	CG = 0.5*sqrt(((jab + mb*mab*0.5 - 2.0) * (jab + mb*0.5*mab))/       
		      ((jab - 1.0) * jab));
      else if (ja == jab) 
	CG = -mb * 0.5 * sqrt(0.5 * ((jab - mb*mab*0.5 + 2.0) * 
				     (jab+mb*mab*0.5)) /
			      ((jab + 2.0) * jab));
      else if (ja -2 == jab) 
	CG = 0.5 * sqrt(((jab - mb*mab*0.5 + 2.0) * 
			 (jab - mb*mab*0.5 + 4.0)) /
			((jab + 2.0) * (jab + 3.0)));
    }
    break;
  default:
    CG = ClebschGordan(ja, ma, jb, mb, jab, mab);
  }
  return CG;
}
  
/*
** FUNCTION:    PackRCFPState, UnPackRCFPState
** PURPOSE:     when the angular momentum of the subshell is 
**              > 9/2, only 2 electrons are allowed, the final
**              state is packed into an integer using the base 1024.
** INPUT:       
** RETURN:      
** SIDE EFFECT: 
** NOTE:        
*/ 
void UnpackRCFPState(int s, int *j, int *nu, int *J) {
  *J = s%1024;
  *nu = (s = s/1024) % 1024;
  *j = (s/1024) % 1024;
}

int PackRCFPState(int j, int nu, int J) {
  return J + nu*1024 + j*1024*1024;
}

/* 
** FUNCTION:    CoupleOperators
** PURPOSE:     couple two operators to produce a final
**              operator with specified rank.
** INPUT:       {RCFP_OPERATOR *op1, *op2, *op},
**              pointers to the component operators, 
**              and the coupled operator.
**              {int rank},
**              rank of the coupled operator.
** RETURN:      none.
** SIDE EFFECT: 
** NOTE:        
*/
void CoupleOperators(RCFP_OPERATOR *op1,
		     RCFP_OPERATOR *op2,
		     RCFP_OPERATOR *op, int rank) {
  op->rank = rank;
  op->left = op1;
  op->right = op2;
  op->nops = op1->nops + op2->nops;
  op->qm = op1->qm + op2->qm;
  return;
}

/* this routine may be used to calculate any operators including those can
   be calculated by ReducedW, reducedWxA etc. althogh using the latter may
   be slightly faster for these particular cases. */
/* 
** FUNCTION:    ReducedOperator
** PURPOSE:     reduced matrix elements of a general 
**              coupled operator.
** INPUT:       {RCFP_STATE *bra, *ket},
**              pointers to the bra and ket states.
**              {RCFP_OPERATOR *op},
**              pointer to the coupled operator.
** RETURN:      {double},
**              result.
** SIDE EFFECT: 
** NOTE:        this routine may be used to calculate any 
**              operators including those can be calculated by 
**              ReducedW, reducedWxA etc. althogh using the latter 
**              is faster for these particular cases.
*/
double ReducedOperator(RCFP_STATE *bra, RCFP_STATE *ket,
		       RCFP_OPERATOR *op) {
  double coeff;
  RCFP_STATE run;
  int run_nu, Jrun;
  int Jbra, jbra, Jket, jket;
  double w6j1;
  int min_bra, max_bra;
  int min_ket, max_ket;

  coeff = 0.0;
  if (op->nops == 1) return ReducedA(bra, ket, op->qm);
  if (op->nops == 2) 
    return ReducedW(bra, ket, op->rank / 2,
		    op->left->qm, op->right->qm);
  else {
    if (bra->state < 63 && ket->state < 63) {
      Jbra = terms_jj[bra->state].subshellJ;
      Jket = terms_jj[ket->state].subshellJ;
      jbra = terms_jj[bra->state].j;
      jket = terms_jj[ket->state].j;

      if (!Triangle(Jbra, op->rank, Jket)) return coeff;
      if (QSpaceDelta(bra) == 0) return coeff;
      if (QSpaceDelta(ket) == 0) return coeff;
      if (bra->nq - ket->nq - op->qm) return coeff;
      if (IsEven(op->left->nops)) {
	min_bra = rcfp_min_even[bra->state];
	max_bra = rcfp_max_even[bra->state];
      } else {
	min_bra = rcfp_min_odd[bra->state];
	max_bra = rcfp_max_odd[bra->state];
      }
      if (IsEven(op->right->nops)) {
	min_ket = rcfp_min_even[ket->state];
	max_ket = rcfp_max_even[ket->state];
      } else {
	min_ket = rcfp_min_odd[ket->state];
	max_ket = rcfp_max_odd[ket->state];
      }
      if (min_bra != min_ket) return coeff;
      if (max_bra != max_ket) return coeff;
      
      run.nq = ket->nq + op->right->qm;
      run.subshellMQ = run.nq - (jbra + 1)/2;
      run.n = bra->n;
      for (run.state = min_bra; run.state <= max_bra; run.state++) {
	Jrun = terms_jj[run.state].subshellJ;
	if(QSpaceDelta(&run) == 0) continue;
	if((w6j1 = W6j(op->left->rank, op->right->rank, op->rank,
		       Jket, Jbra, Jrun))) {
	  coeff += (ReducedOperator(bra, &run, op->left) *
		    ReducedOperator(&run, ket, op->right)) * w6j1;
	}
      }

    } else if (bra->state > 63 && ket->state > 63) {      
      UnpackRCFPState(bra->state, &jbra, &run_nu, &Jbra);
      UnpackRCFPState(ket->state, &jket, &run_nu, &Jket);
      if (!Triangle(Jbra, op->rank, Jket)) return coeff;
      run.nq = ket->nq + op->right->qm;
      switch (run.nq) {
      case 0:
	min_bra = 0;
	max_bra = 0;
	break;
      case 1:
	min_bra = jbra;
	max_bra = jbra;
	break;
      case 2:
	min_bra = 0;
	max_bra = 2*jbra - 2;
	break;
      default:
	/*
	printf("improper input for ReducedOperator\n");
	*/
	return 0.0;
      }
      
      for (Jrun = min_bra; Jrun <= max_bra; Jrun += 4) {
	run_nu = run.nq;
	if (run.nq == 2 && Jrun == 0) run_nu = 0;
	run.subshellMQ = run.nq - (jbra + 1)/2;
	run.state = PackRCFPState(jbra, run_nu, Jrun);

	if ((w6j1 = W6j(op->left->rank, op->right->rank, op->rank,
		       Jket, Jbra, Jrun))) {
	  coeff += (ReducedOperator(bra, &run, op->left) *
		    ReducedOperator(&run, ket, op->right)) * w6j1;
	}
      }
    }
  
    coeff *= sqrt(op->rank + 1.0);
    if (IsOdd((Jbra + op->rank + Jket)/2)) coeff = -coeff;
    return coeff;
  } 
}

/* 
** FUNCTION:    RCFPTermIndex
** PURPOSE:     determine the index of RCFP_TERM for a state.
** INPUT:       {int j},
**              angular momentum of the shell.
**              {int nu},
**              seneority of the state.
**              {int Nr},
**              possible additional quantum numbers.
**              {int subshellJ},
**              total angular momentum of the state.
** RETURN:      {int},
**              index.
** SIDE EFFECT: 
** NOTE:        j <= 9/2
*/
int RCFPTermIndex(int j, int nu, int Nr, int subshellJ) {
  int no_min[9] = {0, 0, 2, 0, 5, 0, 11, 0, 25};
  int no_max[9] = {1, 0, 4, 0, 10, 0, 24, 0, 62};
  int i;

  if (j <= 9) {
    for (i = no_min[j-1]; i<= no_max[j-1]; i++) {
      if (terms_jj[i].nu == nu &&
	  terms_jj[i].subshellJ == subshellJ &&
	  terms_jj[i].Nr == Nr)
	return i;
    }
  } else {
    i = PackRCFPState(j, nu, subshellJ);
    return i;
  }
  return -1;
}




