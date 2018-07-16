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

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_ 1

/*************************************************************
  Header for module "config".
  This module generates electron configuations and 
  carries out the angular momentum coupling. 

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

#include "array.h"
#include "global.h"
#include "parser.h"
#include "nucleus.h"
#include "dbase.h"

/*
** STRUCT:      SHELL
** PURPOSE:     a relativistic subshell.
** FIELDS:      {int n},
**              the principle quantum number.
**              {int kappa},
**              the relativistic angular quantum nubmer.
**              {int nq},
**              the occupation number.
** NOTE:        
*/
typedef struct _SHELL_ {
  int n;
  int kappa;
  int nq;
} SHELL;
  
 
/*
** STRUCT:      SHELL_STATE
** PURPOSE:     a shell state after coupling.
** FIELDS:      {int shellJ},
**              the angular momentum of the shell.
**              {int totalJ},
**              the total angular momentum of the shell after coupling.
**              {int nu},
**              the seneority of the state.
**              {int Nr},
**              any additional quantum numbers.
** NOTE:        a SHELL_STATE specify the seniority and the total 
**              angular momentum of a shell with any occupation, along 
**              with the total angular momentum when this shell is 
**              coupled to its next inner shell. if this is the inner 
**              most shell, this angular momentum is the same as the 
**              shell total_j. All angular momenta are represented by 
**              the double of its actual value.
*/
typedef struct _SHELL_STATE_{
  int shellJ;
  int totalJ;
  int nu; 
  int Nr;
} SHELL_STATE;

typedef struct _SHELL_RESTRICTION_ {
  int ns;
  SHELL *shells;
  int op, nq;
} SHELL_RESTRICTION;

/* 
** MACRO:       GetShellJ, GetTotalJ, GetNu, GetNr
** PURPOSE:     convenient macros to access the fields in a SHELL_STATE.
** INPUT:       {SHELL_STATE s},
**              a shell state.
** RETURN:      {int},
**              results.
** SIDE EFFECT: 
** NOTE:        
*/
#define GetShellJ(s) ((s).shellJ)
#define GetTotalJ(s) ((s).totalJ)
#define GetNu(s)     ((s).nu)
#define GetNr(s)     ((s).Nr)

/*
** STRUCT:      CONFIG
** PURPOSE:     a configuration.
** FIELDS:      {int n_electrons},
**              number of electrons of the configuration.
**              {n_shells},
**              number of the subshells.
**              {n_csfs},
**              number of states resulting from coupling.
**              {SHELL *shells},
**              the shell structure.
**              {SHELL_STATE *csfs},
**              a list specifying all states.
** NOTE:        shells and csfs have the shells in reverse order,
**              i.e., the outmost shell is in the beginning of the list.
*/
typedef struct _CONFIG_ {
  int n_electrons;
  int n_shells;
  int n_csfs;
  int nnrs;
  int igroup;
  int icfg;
  double sweight;
  double energy;
  double delta;
  double sth;
  int *nrs;
  int *symstate;
  SHELL *shells;
  SHELL_STATE *csfs; 
} CONFIG;


/*
** STRUCT:      AVERAGE_CONFIG
** PURPOSE:     the mean configuration for the determinaiton 
**              of the central potential.
** FIELDS:      {int n_cfgs},
**              the number of actual configurations which determined 
**              the mean configuration.
**              {int n_shells},
**              the number of subshells in the mean configuration.
**              {int *n, *kappa, *nq},
**              lists specifing all the shells.
** NOTE:        
*/
typedef struct _AVERAGE_CONFIG_ {
  int n_cfgs;
  int n_shells;
  int *n;
  int *kappa;
  double *nq;
  int ng;
  int *kg;
  double *weight;
} AVERAGE_CONFIG;

/*
** STRUCT:      CONFIG_GROUP
** PURPOSE:     a group of configurations.
** FIELDS:      {char name[]},
**              a string identifies the group.
**              {int n_cfgs},
**              number of configurations in the group.
**              {int n_electrons},
**              number of electrons in the configurations.
**              {ARRAY cfg_list},
**              array of all configurations in the group.
** NOTE:        all configurations in a group must have the 
**              same number of electrons.
*/
typedef struct _CONFIG_GROUP_ { 
  int n_cfgs;
  int n_electrons;
  int nmax;
  ARRAY cfg_list;
  double sweight;
  char name[GROUP_NAME_LEN]; 
} CONFIG_GROUP;

/*
** STRUCT:      STATE
** PURPOSE:     a basis state.
** FIELDS:      {int kgroup},
**              which configuration group the state belongs to.
**              which configuration the state belongs to.
**              which state in all the states the configuration 
**              generates.
** NOTE:        if the state is generated by adding one spectator
**              electron to an existing state, then 
**              kgroup = -(k+1), where k is the index of the parent
**              state in the SYMMETRY array. kcfg = ko, where ko is the 
**              index of the spectator orbital in the orbital array.
**              and kstate = tj, where tj is the total angular momentum
**              of the state.
*/
typedef struct _STATE_ {
  int kgroup;
  int kcfg;
  int kstate;
} STATE;

/*
** STRUCT:      SYMMETRY
** PURPOSE:     a j-parity symmetry.
** FIELDS:      {int n_states},
**              number of states in the symmetry.
**              {ARRAY state},
**              an array of states in the symmetry.
** NOTE:        
*/
typedef struct _SYMMETRY_ {
  int n_states;
  ARRAY states;
} SYMMETRY;

int ShellDegeneracy(int g, int nq);
int          ShellsFromString(char *scfg, double *dnq, SHELL **shell);
int          ShellsFromStringNR(char *scfg, double *dnq, SHELL **shell);
int          GetRestriction(char *scfg, SHELL_RESTRICTION **sr, int m);
int          ApplyRestriction(int ncfg, CONFIG *cfg, int nc, SHELL_RESTRICTION *sr);
int          DistributeElectrons(CONFIG **cfg, double *nq, char *scfg);
int          DistributeElectronsNR(CONFIG **cfg, char *scfg);
int          GetConfigOrAverageFromString(CONFIG **cfg, 
					  double **nq, char *scfg);
int          GetConfigFromStringNR(CONFIG **cfg, char *scfg);
int          GetConfigFromString(CONFIG **cfg, char *scfg);
int          GetAverageConfigFromString(int **n, int **kappa, 
					double **nq, char *scfg);
int          Couple(CONFIG *cfg);
int          CoupleOutmost(CONFIG *cfg, CONFIG *outmost, CONFIG *inner);
int          GetSingleShell(CONFIG *cfg);
void         UnpackShell(SHELL *s, int *n, int *kl, int *j, int *nq);
void         PackShell(SHELL *s, int n, int kl, int j, int nq); 
void         UnpackNRShell(int *s, int *n, int *kl, int *nq);
void         PackNRShell(int *s, int n, int kl, int nq); 
int          PackSymState(int s, int k);
void         UnpackSymState(int st, int *s, int *k);
int          GetJ(SHELL *shell);
int          GetL(SHELL *shell);
int          GetNq(SHELL *shell);
void         GetJLFromKappa(int kappa, int *j, int *kl);
int          GetLFromKappa(int kappa); 
int          GetJLFromSymbol(char *s, int *j, int *kl);
int          GetJFromKappa(int kappa);
int          GetKappaFromJL(int j, int kl); 
int          CompareShell(const void *s1, const void *s2);
int          CompareShellInvert(const void *s1, const void *s2);
int          ShellClosed(SHELL *s);
int          ShellToInt(int n, int k);
int          ShellIndex(int n, int kappa, int ns, SHELL *s);
void         IntToShell(int i, int *n, int *k);
void         PackShellState(SHELL_STATE *s, int J, int j, int nu, int Nr);
int ShellNeedNuNr(SHELL *s, SHELL_STATE *st);
int          GetAverageConfig(int ng, int *kg, int ic, double *weight,
			      int n_screen, int *screened_n, 
			      double screened_charge,
			      int screened_kl, AVERAGE_CONFIG *acfg);
int          GroupIndex(char *name);
int          GroupExists(char *name);
int          AddConfigToList(int k, CONFIG *cfg);
int          AddGroup(char *name);
int          RemoveGroup(int k);
CONFIG_GROUP *GetGroup(int k);
CONFIG_GROUP *GetNewGroup(void);
int          GetNumGroups(void);
int          GetNumConfigs(void);
int          ConfigParity(CONFIG *c);
CONFIG       *GetConfig(STATE *s);
CONFIG       *GetConfigFromGroup(int kg, int kc);
int          AddStateToSymmetry(int kg, int kc, int kstate, 
				int parity, int j);
int          AddConfigToSymmetry(int kg, int kc, CONFIG *cfg);
SYMMETRY     *GetSymmetry(int k);
void         DecodePJ(int i, int *p, int *j);
int          SpecSymbol(char *s, int kl);
int          ConstructConfigName(char *s, int n, CONFIG *c);
void         ListConfig(char *fn, int n, int *kg);
int ReadConfig(char *fn, char *c);
int          InGroups(int kg, int ng, int *kgroup);
int          InitConfig(void);
int          ReinitConfig(int m);
int          SetNCG(void);
int ConfigToIList(CONFIG *c, int n, int *s);
CONFIG *ConfigFromIList(int n, int *s);
int ConfigExists(CONFIG *c);
void SetClosedShellNR(int n, int k);
int IsClosedShellNR(int n, int k, int nq);
int IsClosedShellFR(int n, int k, int j, int nq);
int IsClosedComplex(int n, int nq);

#endif
