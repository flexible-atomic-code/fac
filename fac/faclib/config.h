#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_ 1

/*************************************************************
  Header for module "config".
  This module generates electron configuations and 
  carries out the angular momentum coupling. 

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

#include "array.h"
#include "global.h"
#include "nucleus.h"

/* maximum orbital angular momentum for which a symbol is available */
#define MAX_SPEC_SYMBOLS   14 
/* max string length for the level name */
#define LEVEL_NAME_LEN     128
/* max string length for the group name */
#define GROUP_NAME_LEN     32
/* max number of configuration groups */
#define MAX_GROUPS         1024
/* max number of j-parity symmetries */
#define MAX_SYMMETRIES     512
/* number of groups in one array block */
#define CONFIGS_BLOCK      256
/* number os states in one array block */
#define STATES_BLOCK       512

/* a SHELL is a set of qunatum #s, the principal quantum number n, 
   the relativistic angular quantum number kappa. */
typedef struct _SHELL_ {
  int n;
  int kappa;
  int nq;
} SHELL;
  
/* a SHELL_STATE specify the seniority and the total angular momentum of 
   a shell with any occupation, along with the total angular momentum when
   this shell is coupled to its next inner shell. if this is the inner most 
   shell, this angular momentum is the same as the shell total_j. 
   All angular momenta are represented by the double of its actual value.*/  
typedef struct _SHELL_STATE_{
  int shellJ;
  int totalJ;
  int nu; 
  int Nr; /* Nr is the additional quantum # may need to specify the state */
} SHELL_STATE;
#define GetShellJ(s) ((s).shellJ)
#define GetTotalJ(s) ((s).totalJ)
#define GetNu(s)     ((s).nu)
#define GetNr(s)     ((s).Nr)

/* in config, shells are ordered in reverse order, i.e., outer shells come
   first.  */
typedef struct _CONFIG_ {
  int n_electrons;
  int n_shells;
  int n_csfs;
  SHELL *shells;
  SHELL_STATE *csfs; /* all coupled states */
} CONFIG;

/* this is the mean configuration for the 
   determinaiton of central potential */
typedef struct _AVERAGE_CONFIG_ {
  int n_cfgs;
  int n_shells;
  int *n;
  int *kappa;
  double *nq;
} AVERAGE_CONFIG;

typedef struct _CONFIG_GROUP_ {
  char name[GROUP_NAME_LEN];  
  int n_cfgs;
  int n_electrons;
  ARRAY cfg_list;
} CONFIG_GROUP;

typedef struct _STATE_ {
  int kgroup;
  int kcfg;
  int kstate;
} STATE;

typedef struct _SYMMETRY_ {
  int n_states;
  ARRAY states;
} SYMMETRY;

int          Couple(CONFIG *cfg);
int          CoupleOutmost(CONFIG *cfg, CONFIG *outmost, CONFIG *inner);
int          GetSingleShell(CONFIG *cfg);
void         UnpackShell(SHELL *s, int *n, int *kl, int *j, int *nq);
void         PackShell(SHELL *s, int n, int kl, int j, int nq); 
int          GetJ(SHELL *shell);
int          GetL(SHELL *shell);
int          GetNq(SHELL *shell);
void         GetJLFromKappa(int kappa, int *j, int *kl);
int          GetLFromKappa(int kappa); 
int          GetJLFromSymbol(char *s, int *j, int *kl);
int          GetJFromKappa(int kappa);
int          GetKappaFromJL(int j, int kl); 
int          CompareShell(SHELL *s1, SHELL *s2);
int          ShellClosed(SHELL *s);
void         PackShellState(SHELL_STATE *s, int J, int j, int nu, int Nr);
int          GetConfigFromString(CONFIG **cfg, char *scfg);
int          GetAverageConfigFromSTring(int **n, int **kappa, 
					double **nq, char *scfg);
int          GetAverageConfig(int ng, int *kg, double *weight,
			      int n_screen, int *screened_n, 
			      double screened_charge,
			      int screened_kl, AVERAGE_CONFIG *acfg);
int          GroupIndex(char *name);
int          GroupExists(char *name);
int          AddConfigToList(int k, CONFIG *cfg);
int          AddGroup(char *name);
CONFIG_GROUP *GetGroup(int k);
CONFIG_GROUP *GetNewGroup();
int          GetNumGroups();
CONFIG       *GetConfig(STATE *s);
int          AddStateToSymmetry(int kg, int kc, int kstate, 
				int parity, int j);
int          AddConfigToSymmetry(int kg, int kc, CONFIG *cfg);
SYMMETRY     *GetSymmetry(int k);
void         DecodePJ(int i, int *p, int *j);
int          SpecSymbol(char *s, int kl);
int          InGroups(int kg, int ng, int *kgroup);
int          InitConfig();

#endif
