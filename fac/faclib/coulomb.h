#ifndef _COULOMB_H_
#define _COULOMB_H_ 1

/*************************************************************
  Header for module "coulomb". 
  This module calculates quatities related to the H-like ions.

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

#include "global.h"
#include "array.h"
#include "angular.h"
#include "grid.h"

#define CBMULTIPOLES   2
#define MAXNCB         (CBMULTIPOLES*(CBMULTIPOLES+3)/2)

void    SetHydrogenicNL(int n, int kl, int nm, int klm);
void    GetHydrogenicNL(int *n, int *kl, int *nm, int *klm);
double  HydrogenicDipole(double z, int n0, int kl0, 
			int n1, int kl1);
double HydrogenicExpectation(double z, int m, int n, int kl);
double HydrogenicSelfEnergy(double z, int n, int k);
double  TRRateHydrogenic(double z, int n0, int kl0,
			int n1, int kl1, int s);
double  CoulombPhaseShift(double z, double e, int kappa);
double  AngularMSub(int lf, int li1, int li2, int q);
double *GetCoulombBethe(int ie2, int ite, int ie1, int m, int k);
double  GetCoulombBetheAsymptotic(double te, double e1);
int     PrepCBIndex(int mode);
int     CoulombBetheTail(int n, double *w, int nkl, double *kl, double *tcb);
int     PrepCoulombBethe(int ne2, int nte, int ne1, double z,
			 double *e2, double *te, double *e1,
			 int nkl, double *kl, 
			 int etype, int ltype, int mode);
int     CoulombBethe(char *s, double z, double te, double e1);
int     InitCoulomb(void);

#endif
