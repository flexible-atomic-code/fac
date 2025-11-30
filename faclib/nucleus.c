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

#include "nucleus.h"
#include "cf77.h"
#include "coulomb.h"
#include "errms.h"
#include "grdcfg.h"
#include "grdcfg1.h"
#include "global.h"

static char *rcsid="$Id: nucleus.c,v 1.14 2005/01/17 05:39:41 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static NUCLEUS atom;
static CXTGT cxtgt;

#define N_CXT 19
#define N_CXS 32
static char _cxtname[N_CXT][N_CXS] = {
				      "H", "He", "H2", "H2O", "CO", "CO2", "O2", "N2", "Ne", "Ar", "Kr", "Xe", "CH4", "C2H6O", "O", "F", "S", "CS2", "SF6"
};

static double _polarizability[N_CXT] = {
					4.5, 1.383, 5.5, 9.82, 12.05, 19.3, 10.4, 11.16, 2.67, 11.08, 16.77, 27.34, 17.24, 36.5, 5.41, 3.76, 19.59, 65.8, 4.49
};

static double _screening[N_CXT] = {
				   1.8, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

static double _ionpot[N_CXT] = {
				0.49973, 0.90357, 0.567, 0.4638, 0.515, 0.5063, 0.44355, 0.5726, 0.7925, 0.5792, 0.5145, 0.4458, 0.46341, 0.38513, 0.50045, 0.64028, 0.49247, 0.370176, 0.563
};

static char _ename[N_ELEMENTS][3] = 
{"H", "He", "Li", "Be", "B", "C", "N", "O", "F",
 "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
 "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", 
 "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
 "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", 
 "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
 "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
 "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
 "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
 "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
 "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
 "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Ux", "Uy"};

static double _emass[N_ELEMENTS] = 
{1.008, 4.003, 6.94, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00,
 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95, 39.10, 
 40.08, 44.96, 47.88, 50.94, 52.00, 54.94, 55.85, 58.93, 58.69, 63.55,
 65.39, 69.72, 72.64, 74.92, 78.96, 79.90, 83.79, 85.47,
 87.62, 88.91, 91.22, 92.91, 95.96, 98.00, 101.1, 102.9,
 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9,
 131.3, 132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0,
 150.4, 152.0, 157.2, 158.9, 162.500001, 164.9, 167.3,
 168.9, 173.0, 175.0, 178.5, 180.9, 183.9, 186.2, 190.2,
 192.2, 195.1, 197.0, 200.5, 204.38, 207.2, 209.0, 
 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0,
 231.0, 238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 
 252.0, 257.0, 258.0, 259.0, 262.0, 267.0, 268.0, 269.0, 270.0,
 277.0, 278., 281.0, 282.0, 285.0, 286.0, 289.0, 289.0, 293.0, 294.0, 294.0, 296.0, 296.0};

static double _arrms[N_ELEMENTS] =
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   3.057, 3.061, 3.122, 3.189, 3.261, 3.365, 3.427, 3.435,
   3.476, 3.544, 3.591, 3.599, 3.642, 3.706, 3.737, 3.788,
   3.775, 3.882, 3.929, 3.997, 4.074, 4.097, 4.140, 4.163,
   4.188, 4.203, 4.220, 4.242, 4.270, 4.324, 4.409, 4.424,
   4.482, 4.494, 4.532, 4.544, 4.614, 4.617, 4.654, 4.680,
   4.743, 4.750, 4.787, 4.804, 4.839, 4.855, 4.877, 4.892,
   4.912, 4.962, 5.084, 5.113, 5.162, 5.060, 5.221, 5.202,
   5.251, 5.226, 5.312, 5.370, 5.342, 5.351, 5.367, 5.339,
   5.413, 5.402, 5.428, 5.436, 5.463, 5.476, 5.501, 5.521,
   5.526, 5.539, 5.655, 5.658, 5.684, 5.670, 5.710, 5.700,
   5.851, 5.744, 5.864, 5.905, 5.815, 5.815, 5.843, 5.850, 5.857};

static double _mserms[N_ELEMENTS] = {
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  3.0055    ,  2.9936    ,  3.057     ,  3.061     ,  3.1224    ,
  3.1889    ,  3.2611    ,  3.365     ,  3.4274    ,  3.4349    ,
  3.4776    ,  3.5459    ,  3.5921    ,  3.6002    ,  3.6452    ,
  3.7057    ,  3.7377    ,  3.7875    ,  3.7757    ,  3.8823    ,
  3.9283    ,  3.9973    ,  4.0742    ,  4.0968    ,  4.14      ,
  4.1629    ,  4.1884    ,  4.2036    ,  4.224     ,  4.243     ,
  4.2694    ,  4.324     ,  4.4091    ,  4.424     ,  4.4809    ,
  4.4945    ,  4.5318    ,  4.5454    ,  4.5944    ,  4.6156    ,
  4.6519    ,  4.6802    ,  4.7423    ,  4.75      ,  4.7859    ,
  4.8041    ,  4.8378    ,  4.855     ,  4.8771    ,  4.8919    ,
  4.9123    ,  4.91548125,  4.85956687,  4.88787188,  4.93128563,
  5.07265   ,  5.0898425 ,  5.17599   ,  5.238471  ,  5.303984  ,
  5.3108    ,  5.302875  ,  5.31358125,  5.31725813,  5.31381881,
  5.3230125 ,  5.331411  ,  5.37975   ,  5.41033488,  5.44389637,
  5.4648    ,  5.48103366,  5.49260437,  5.49895157,  5.51991148,
  5.53461255,  5.51322061,  5.54938223,  5.5893442 ,  5.65759688,
  5.6849025 ,  5.78639063,  5.8118025 ,  5.89141438,  5.88691415,
  5.90087899,  5.91494926,  5.91168879,  5.88264623,  5.86691016,
  5.85471211,  5.84942109,  5.86143996,  5.87414508,  5.89325189,
  5.919     ,  5.92795102,  5.94922305,  5.95825928,  5.97967266,
  5.99007373,  5.97927161,  5.97652411,  5.98987011,  6.02529002,
  6.08752438,  6.13926806,  6.15654237,  6.15927535,  6.16201645,
  6.17454773};
 
static double _errms[N_ELEMENTS][NISO];

SETUPGROUND
SETUPGROUND1

static double _xfermi0 = XFERMI0;
static double _xfermi1 = XFERMI1;
static int _nfermi = NFERMI;
static double *_xfermi = NULL;
static double *_yfermi = NULL;
static double *_rfermi[5] = {NULL,NULL,NULL,NULL,NULL};
static double _afermi = 2.3;

int InitNucleus() {
  int z, k;
  
  for (z = 0; z < N_ELEMENTS; z++) {
    for (k = 0; k < NISO; k++) {
      _errms[z][k] = 0.0;
    }
  }
  SETUPERRMS;
  
  SetupFermi();
  atom.atomic_number = 0.0;
  atom.nepr = 0;
  atom.epr = NULL;
  atom.epv = NULL;
  SetExtraPotential(-1, 0, NULL, NULL);
  return 0;
}

void SetupFermi(void) {
  int i, k;
  
  if (_xfermi != NULL) {
    free(_xfermi);
    free(_yfermi);
    for (k = 0; k <= 4; k++) {
      free(_rfermi[k]);
    }
  }  
  _xfermi = malloc(sizeof(double)*_nfermi);
  _yfermi = malloc(sizeof(double)*_nfermi);
  for (k = 0; k <= 4; k++) {
    _rfermi[k] = malloc(sizeof(double)*_nfermi);
  }
  
  _xfermi[0] = _xfermi0;
  double dx = (_xfermi1-_xfermi0)/(_nfermi-1);
  for (i = 1; i < _nfermi; i++) {
    _xfermi[i] = _xfermi[i-1] + dx;
  }
  for (i = 0; i < _nfermi; i++) {
    _yfermi[i] = dx/(1 + exp(_xfermi[i]));
  }
  for (k = 0; k <= 4; k++) {
    _rfermi[k][0] = 0.0;
    NewtonCotes(_rfermi[k], _yfermi, 0, _nfermi-1, -1, 0);
    if (k < 4) {
      for (i = 0; i < _nfermi; i++) {
	_yfermi[i] *= _xfermi[i];
      }
    }
  }
}

void SetExtraPotential(int m, int n, double *p, char *fn) {
  int i;
  if (m < 0) {
    atom.nep = 0;
    for (i = 0; i < NEP; i++) {
      atom.epm[i] = m;      
    }
    if (atom.nepr > 0) {
      free(atom.epr);
      free(atom.epv);
      atom.epr = NULL;
      atom.epv = NULL;
      atom.nepr = 0;
      atom.cepr = 0;
    }
  } else {
    if (fn != NULL) {
      FILE *f = fopen(fn, "r");
      if (f == NULL) {
	printf("cannot open extra potential file: %s\n", fn);
	return;
      }
      if (atom.nepr > 0) {
	free(atom.epr);
	free(atom.epv);
      }
      char buf[1024];
      atom.nepr = 0;
      while (1) {
	if (NULL == fgets(buf, 1024, f)) break;
	StrTrim(buf, '\0');
	if (buf[0] == '#') {
	  continue;
	}
	atom.nepr++;
      }
      if (atom.nepr == 0) {
	printf("empty extra potential file: %s\n", fn);
	fclose(f);
	return;
      }
      atom.epr = malloc(sizeof(double)*atom.nepr);
      atom.epv = malloc(sizeof(double)*atom.nepr);
      fseek(f, 0, SEEK_SET);
      i = 0;
      while (1) {
	if (NULL == fgets(buf, 1024, f)) break;
	StrTrim(buf, '\0');
	if (buf[0] == '#') {
	  continue;
	}
	sscanf(buf, "%lg %lg", atom.epr+i, atom.epv+i);
	i++;
      }
      fclose(f);
      atom.cepr = m;
      return;
    }
    i = atom.nep;
    if (i == NEP) {
      printf("extra potential terms exceeded max: %d\n", NEP);
      return;
    }
    atom.epm[i] = m;
    if (m >= 0 && n > 0 && p != NULL) {
      memcpy(atom.epp[i], p, sizeof(double)*Min(n, NEPP));
    }
    if (m == 2 || m == 3) {
      double r0 = 3*(atom.atomic_number-atom.epp[i][0]);
      r0 /= (FOUR_PI*atom.epp[i][1]);
      r0 = pow(r0, ONETHIRD)/RBOHR;
      atom.epp[i][1] = r0;
    }
    atom.nep++;
  }
}

char *GetAtomicSymbolTable(void) {
  return (char *) _ename;
}

double *GetAtomicMassTable(void) {
  return _emass;
}

void IntegrateFermi(int nk, double *r, double x) {
  int k, one=1, np = 1, n = _nfermi;
  
  if (x >= _xfermi[0]) {
    for (k = 0; k < nk; k++) {
      UVIP3P(np, n, _xfermi, _rfermi[k], one, &x, &r[k]);
    }
  } else {
    double x1 = x;
    double x0 = _xfermi[0];
    r[0] = x1 - x0;
    if (nk > 1) {
      x1 *= x;
      x0 *= _xfermi[0];
      r[1] = (x1 - x0)/2.0;
      if (nk > 2) {
	x1 *= x;
	x0 *= _xfermi[0];
	r[2] = (x1 - x0)/3.0;
	if (nk > 3) {
	  x1 *= x;
	  x0 *= _xfermi[0];
	  r[3] = (x1 - x0)/4.0;
	  if (nk > 4) {
	    x1 *= x;
	    x0 *= _xfermi[0];
	    r[4] = (x1 - x0)/5.0;
	  }
	}
      }
    }
  }
}

double DiffRRMS(double c, double a, double a2, double a3, double a4, double r2,
		double y0, double y1, double y2, double y3, double y4) {
  double f, c2, c3, c4, r0, r1;
  c2 = c*c;
  c3 = c2*c;
  c4 = c3*c;
  IntegrateFermi(5, atom.rfermi, -c/a);
  y0 -= atom.rfermi[0];
  y1 -= atom.rfermi[1];
  y2 -= atom.rfermi[2];
  y3 -= atom.rfermi[3];
  y4 -= atom.rfermi[4];
  r0 = a2*y2 + 2*a*c*y1 + c2*y0;
  r1 = a4*y4 + 4*a3*c*y3 + 6*a2*c2*y2 + 4*a*c3*y1 + c4*y0;
  
  return r1/r0 - r2;
}

double FermiRMS(double c, double a) {
  double y[5], r[5];
  double a2, a3, a4, c2, c3, c4;
  int i, n;

  if (a < 0) {
    a = _afermi*1e-5/RBOHR/(4*log(3.0));
  }
  
  n = _nfermi-1;
  a2 = a*a;
  a3 = a2*a;
  a4 = a3*a;
  c2 = c*c;
  c3 = c2*c;
  c4 = c3*c;

  IntegrateFermi(5, r, -c/a);
  for (i = 0; i < 5; i++) {
    y[i] = _rfermi[i][n] - r[i];
  }
  double r0 = a2*y[2] + 2*a*c*y[1] + c2*y[0];
  double r1 = a4*y[4] + 4*a3*c*y[3] + 6*a2*c2*y[2] + 4*a*c3*y[1] + c4*y[0];
  return sqrt(r1/r0);
}

double FermiParamC(double rn, double a) {
  int i;

  if (a < 0) {
    a = _afermi*1e-5/RBOHR/(4*log(3.0));
  }
  
  double c0 = 1e-2*rn;
  double c1 = 10*(rn+a);
  double c, r;
  
  for (i = 0; i < 500; i++) {
    if (fabs(c1/c0-1) < 1e-10) break;
    c = 0.5*(c0+c1);
    r = FermiRMS(c, a);
    if (r < rn) {
      c0 = c;
    } else if (r > rn) {
      c1 = c;
    } else {
      break;
    }
  }
  if (i == 500) {
    printf("max iter in finding fermi c param: %d %g %g %g %g %g\n",
	   i, rn, a, c0, c1, c);
  }
  return c;
}

double GraspRRMS(double z, double m) {
  int iz, ia, i, k;
  double r0;
  iz = (int)(z-1);
  ia = (int)(1.5 + m - 2*(iz+1));
  if (ia >= 0 && ia < NISO) {
    r0 = _errms[iz][ia];
    if (r0 > 0) return r0;
  }
  for (i = 1; i < 10; i++) {
    k = ia-i;
    if (k >= 0 && k < NISO) {
      r0 = _errms[iz][k];
      if (r0 > 0) {
	r0 += 0.836*(pow(m,0.333)-pow(m-i,0.333));
	return r0;
      }
    }
    k = ia+i;
    if (k >= 0 && k < NISO) {
      r0 = _errms[iz][k];
      if (r0 > 0) {
	r0 += 0.836*(pow(m,0.333)-pow(m+i,0.333));
	return r0;
      }
    }
  }
  return 0.0;
}

int SetAtom(char *s, double z, double mass, double rn, double a, double rmse) {
  int i;
  char un[3] = "Xx";
  if (s == NULL || strlen(s) == 0) {
    if (z <= 0) {
      printf("atomic symbol and z cannot be both unset\n");
      return -1;
    }
    int iz = (int)z;
    if (iz <= N_ELEMENTS) {
      s = _ename[iz-1];
    } else {
      s = un;
    }
  } else {
    int iz = atoi(s);
    if (iz > 0) {      
      z = (double)iz;
      if (iz <= N_ELEMENTS) {
	s = _ename[iz-1];
      } else {
	s = un;
      }
    }
  }
  strncpy(atom.symbol, s, 2);
  atom.z0 = 0;
  atom.m0 = 0;
  for (i = 0; i < N_ELEMENTS; i++) {
    if (strncasecmp(_ename[i], s, 2) == 0) {
      atom.z0 = i+1;
      if (z <= 0) atom.atomic_number = i+1;
      atom.m0 = _emass[i];
      if (mass <= 0) atom.mass = _emass[i];
      break;
    }
  }
  if (z <= 0 || mass <= 0) {
    if (i == N_ELEMENTS) {
      printf("unknown element must have z and mass set\n");
      return -1;
    }
  }
  if (z > 0) {
    atom.atomic_number = z;
  } 
  if (mass > 0.0) {
    atom.mass = mass;
  }
  atom.rms0 = 1e-5*(0.570 + 0.836*pow(atom.mass,0.3333333333))/RBOHR;
  if (rn < 0.0) {
    if (rn > -1.5) {
      atom.rn = NucleusRRMS(atom.atomic_number);
    } else if (rn > -3.5) {
      if (atom.m0 > 0) {
	atom.rn = GraspRRMS(atom.atomic_number, atom.m0);
      } else {
	atom.rn = GraspRRMS(atom.atomic_number, atom.mass);
      }
      if (atom.rn <= 0) {
	if (rn > -2.5) {
	  atom.rn = NucleusRRMS(atom.atomic_number);
	} else {
	  atom.rn = 0.570 + 0.836*pow(atom.mass,0.3333333333);
	}
      }
    } else {
      int iz = (int)(atom.atomic_number-0.5);
      if (iz >= 0 && iz < N_ELEMENTS) {
	atom.rn = _arrms[iz];
      }
      if (atom.rn <= 0) {
	atom.rn = NucleusRRMS(atom.atomic_number);
      }
    }
    if (atom.rn <= 0) {
      atom.rn = 0.570 + 0.836*pow(atom.mass,0.3333333333);
    }
    if (atom.mass != atom.m0 && atom.m0 > 0) {
      atom.rn += 0.836*(pow(atom.mass,0.3333333)-pow(atom.m0,0.3333333));
    }
  } else {
    atom.rn = rn;
  }

  SetAtomicChargeDist(a, rmse);
  return 0;
}

void SetAtomicChargeDist(double a, double rmse) {
  if (a >= 0) {
    atom.a = a*1e-5/(RBOHR*4*log(3.0));
  } else {
    atom.a = _afermi*1e-5/(RBOHR*4*log(3.0));
  }
  if (atom.a > 0 && atom.rn > 0) {
    if (atom.rn > 100.0) {
      atom.c = (atom.rn-100.0)*1e-5/RBOHR;
      atom.rn = FermiRMS(atom.c, atom.a);
    } else {
      atom.rn *= 1e-5/RBOHR;
      atom.c = FermiParamC(atom.rn, atom.a);
    }
  } else if (atom.rn > 0) {
    atom.rn *= 1e-5/RBOHR;
  }
  atom.rms = atom.rn;

  if (atom.a <= 0 && atom.rn > 0) {
    atom.rn = sqrt(5.0/3.0)*atom.rn;
    atom.z1 = 1.5*atom.atomic_number/atom.rn;
  }
  
  if (atom.rn > 0 && atom.a > 0) {
    int i = _nfermi-1;
    double a2 = atom.a*atom.a;
    double c2 = atom.c*atom.c;
    double ac = atom.a*atom.c;
    IntegrateFermi(3, atom.rfermi, -atom.c/atom.a);    
    atom.b = c2*(_rfermi[0][i]-atom.rfermi[0]);
    atom.b += 2*ac*(_rfermi[1][i]-atom.rfermi[1]);
    atom.b += a2*(_rfermi[2][i]-atom.rfermi[2]);  
    atom.b = atom.atomic_number/atom.b;
    atom.z1 = atom.c*(_rfermi[0][i]-atom.rfermi[0]);
    atom.z1 += atom.a*(_rfermi[1][i] - atom.rfermi[1]);
    atom.z1 *= atom.b;
  }

  if (rmse < 0) {
    int iz = (int)atom.atomic_number;
    if (rmse < 0 && iz <= N_ELEMENTS) {
      rmse = _mserms[iz-1];    
      if (rmse <= 0) {
	atom.rmse = atom.rms*1e5*RBOHR;
      } else {
	atom.rmse = rmse;
      }
    }
  } else {
    atom.rmse = rmse;
  }

  SetExtraPotential(-1, 0, NULL, NULL);
  
  if (atom.atomic_number >= 10 &&
      atom.atomic_number <= 120) {
    INIQED(atom.atomic_number, 9, atom.rn>0, atom.rmse);
  }
}

/* dimensions:
** rn -- L
** z1 -- 1/L
** a  -- L
** b  -- 1/L^2
** c  -- L
*/
void ScaleAtomicChargeDist(double a) {
  atom.rn *= a;
  atom.z1 /= a;
  atom.rms *= a;
  atom.rms0 *= a;
  atom.a *= a;
  atom.b /= a*a;
  atom.c *= a;
}  
  
double NucleusRadius(double z, double m, int md) {
  int iz = (int)(z-0.5);
  double m0 = _emass[iz];
  if (m < 0) m = m0;
  double r;
  double r0 = 0.570 + 0.836*pow(m0,0.3333333);
  if (md >= -1) {
    r = NucleusRRMS(z);
  } else if (md >= -3) {
    r = GraspRRMS(z, m0);
    if (r <= 0) {
      if (md == -2) {
	r = NucleusRRMS(z);
      } else {
	r = r0;
      }
    }
  } else {
    r = _arrms[iz];
    if (r <= 0) {
      r = NucleusRRMS(z);
    }
  }
  if (m > 0 && m != m0) {
    r += 0.836*(pow(m, 0.3333333)-pow(m0, 0.3333333));
  }
  return r;
}
  
NUCLEUS *GetAtomicNucleus() {
  return &atom;
}

void PrintNucleus(int m, char *fn) {
  FILE *f;

  if (fn == NULL || strlen(fn) == 0 || strcmp(fn, "-")==0) {
    f = stdout;
  } else {
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open iso output file: %s\n", fn);
      return;
    }
  }
  if (m == 0) {
    fprintf(f, "atom: %s\n", atom.symbol);
    fprintf(f, "z: %g\n", atom.atomic_number);
    fprintf(f, "mass: %g\n", atom.mass);
    fprintf(f, "rn: %g\n", atom.rn);
    fprintf(f, "rms: %g fm\n", atom.rms*1e5*RBOHR);
    fprintf(f, "rmse: %g fm\n", atom.rmse);
    fprintf(f, "fz1: %g\n", atom.z1);
    fprintf(f, "fa: %g\n", atom.a);
    fprintf(f, "fb: %g\n", atom.b);
    fprintf(f, "fc: %g\n", atom.c);
    int i;
    for (i = 0; i < atom.nep; i++) {
      fprintf(f, "ep: %d %d %g %g\n",
	      i, atom.epm[i], atom.epp[i][0], atom.epp[i][1]);
    }
  } else {
    fprintf(f, "Atomic number:\n");
    fprintf(f, "%15.8E\n", atom.atomic_number);
    fprintf(f, "Mass number (integer) :\n");
    fprintf(f, "%15.8E\n", (double)(int)(0.5+atom.mass));
    fprintf(f, "Fermi distribution parameter a:\n");
    fprintf(f, "%15.8E\n", atom.a*1e5*RBOHR);
    fprintf(f, "Fermi distribution parameter c:\n");
    fprintf(f, "%15.8E\n", atom.c*1e5*RBOHR);
    fprintf(f, "Mass of nucleus (in amu):\n");
    fprintf(f, "%15.8E\n", atom.mass);
    fprintf(f, "Nuclear spin (I) (in units of h / 2 pi):\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "Nuclear dipole moment (in nuclear magnetons):\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "Nuclear quadrupole moment (in barns):\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "RNT:\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "H:\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "HP:\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "NNNP:\n");
    fprintf(f, "%d\n", 0);
  }
  fflush(f);
  if (f != stdout) {
    fclose(f);
  }
}

double GetAtomicMass(void) {
  return atom.mass;
}

double GetAtomicNumber(void) {
  return atom.atomic_number;
}

char *GetAtomicSymbol(void) {
  return atom.symbol;
}

double GetAtomicChargeDist(double r) {
  double x;
  double r3;
  
  if (atom.rn <= 0) return 0.0;
  if (atom.a <= 0) {
    if (r > atom.rn) return 0.0;
    r3 = atom.rn;
    r3 = r3*r3*r3;
    return 3*atom.atomic_number/r3;
  }

  x = (r-atom.c)/atom.a;
  if (x < -10) r3 = atom.b;
  else if (x > 10) r3 = atom.b*exp(-x);
  else r3 = atom.b/(1+exp(x));
  r3 /= atom.a;
  return r3;
}

double GetExtraZ(double r, int i) {
  double z, zi;
  if (i >= atom.nep) {
    if (atom.nepr > 0) {
      if (r < atom.epr[0]) return atom.epv[0];
      if (r > atom.epr[atom.nepr-1]) return atom.epv[atom.nepr-1];
      UVIP3P(3, atom.nepr, atom.epr, atom.epv, 1, &r, &z);
      return z;
    }
    return 0;
  }
  z = 0.0;
  switch (atom.epm[i]) {
  case 0:
  case 100:    
    zi = atom.epp[i][0]*(atom.mass-atom.atomic_number);
    if (atom.epp[i][1] > 0) {
      zi *= exp(-2.68e-4*atom.epp[i][1]*r);
    }
    z += zi;
    break;
  case 1:
  case 101:
    if (r < atom.rms*atom.epp[i][1]) r = atom.rms*atom.epp[i][1];
    zi = 1.48e-8*(0.76+2.79/pow(atom.mass,0.33333))*atom.mass;
    zi *= atom.epp[i][0]*atom.rms*atom.rms/(r*r*r);
    zi *= FINE_STRUCTURE_CONST;
    z += zi;
    break;
  default:
    break;
  }
  return z;
}

double GetAtomicEffectiveZ(double r) {
  double x, y[3], z;
  int np = 3;
  int n = _nfermi;
  int one = 1;

  if (atom.rn <= 0) return (double)(atom.atomic_number);
  if (atom.a <= 0) {    
    if (r > atom.rn) {
      return (double) atom.atomic_number;
    } else {
      x = r/atom.rn;
      z = 3.0 - x*x;
      z = x*z*0.5*(atom.atomic_number);
      return z;
    }
  }
  x = (r - atom.c)/atom.a;
  if (x >= _xfermi[n-1]) {
    return (double)(atom.atomic_number);
  } else {
    IntegrateFermi(3, y, x);
    z = atom.c*atom.c*(y[0] - atom.rfermi[0]);
    z += 2*atom.c*atom.a*(y[1] - atom.rfermi[1]);
    z += atom.a*atom.a*(y[2] - atom.rfermi[2]);
    z += r*(atom.c*(_rfermi[0][n-1] - y[0]) + atom.a*(_rfermi[1][n-1] - y[1]));
    z *= atom.b;
    return z;
  }
}

double GetAtomicR(void) {
  return atom.rn;
}

//H excited state polarizability.
//K. McDowell, J. Chem. Phys. 65, 2518 (1976)
double HPolarizability(double n) {
  double n2 = n*n;
  return 0.5*n2*n2*(2*n2+7.0);
}

int SetCXTarget(char *s0, double a, double b, double e,
		double x, double z, double m) {
  char *p, *p0, s[128];
  char sa[3];
  int i, k, n;
  double z0, m0;
  
  strncpy(cxtgt.symbol, s0, 128);
  StrTrim(cxtgt.symbol, '\0');
  strncpy(s, cxtgt.symbol, 128);
  cxtgt.a = a;
  cxtgt.b = b;
  cxtgt.e = e/HARTREE_EV;
  cxtgt.x = x;
  cxtgt.z = z;
  cxtgt.m = m;  
  if (a > 0 && b > 0 && e > 0 && z > 0 && m > 0) return 0;
  int np = -1;
  if (strlen(cxtgt.symbol) > 0 &&
      (s[0] == 'H' || s[0] == 'D' || s[0] == 'T')) {
    if (s[1] == '\0') {
      np = 1;
    } else if (s[1] == '_') {
      if (strlen(s) > 2) {
	np = atoi(s+2);
      } else {
	np = 2;
      }
    }
    if (np > 0) {
      double n2 = np*np;
      if (cxtgt.e <= 0) {
	cxtgt.e = _ionpot[0]/n2;
      }
      if (cxtgt.m <= 0) {
	if (s[0] == 'D') {
	  cxtgt.m = 2.0;
	} else if (s[0] == 'T') {
	  cxtgt.m = 3.0;
	}
      }
      if (cxtgt.a <= 0) {
	cxtgt.a = HPolarizability((double)np);
      }
      s[0] = 'H';
      s[1] = '\0';
    }
  }
  for (i = 0; i < N_CXT; i++) {
    if (0 == strncmp(s, _cxtname[i], N_CXS)) {
      if (cxtgt.a <= 0) {
	cxtgt.a = _polarizability[i];
      }
      if (cxtgt.e <= 0) {
	cxtgt.e = _ionpot[i];
      }
      if (cxtgt.b <= 0) {
	cxtgt.b = _screening[i];
	if (cxtgt.b <= 0 || np > 0) {
	  cxtgt.b = 0.936 + (cxtgt.e/_ionpot[0])*(_screening[0]-0.936);
	}
      }
    }
  }
  p0 = s;
  p = s;
  k = 0;
  z0 = 0.0;
  m0 = 0.0;
  n = 1;
  while (1) {
    if (*p) {
      if (isupper(*p)) {
	if (k == 0) {
	  sa[k] = *p;
	  p++;
	  k++;
	  continue;
	} else {
	  sa[k] = '\0';
	}
      } else if (islower(*p)) {
	sa[k] = *p;
	p++;
	k++;
	continue;
      } else if (isdigit(*p)) {
	n = atoi(p);
	p++;
	if (n >= 10) p++;
	if (n >= 100) p++;
	continue;
      }
    } else {
      sa[k] = '\0';
    }
    if (sa[k] == '\0') {      
      for (i = 0; i < N_ELEMENTS; i++) {
	if (strncasecmp(_ename[i], sa, 2) == 0) {
	  z0 += (i+1)*n;
	  m0 += _emass[i]*n;
	  break;
	}
      }
      if (i >= N_ELEMENTS) {
	printf("cannot determine atom: %s\n", sa);
	return -1;
      }
      k = 0;
      n = 1;
    }
    if (! *p) {
      break;
    }
  }
  if (cxtgt.z <= 0) {
    cxtgt.z = z0;
  }
  if (cxtgt.m <= 0) {
    cxtgt.m = m0;
  }
  if (cxtgt.z <= 0 || cxtgt.m <= 0) {
    printf("cannot determine cx target: %s\n", s);
    return -1;
  }
  return 0;
}

CXTGT *GetCXTarget() {
  return &cxtgt;
}

void PrintCXTarget(char *fn) {
  FILE *f;

  if (fn == NULL || strlen(fn) == 0 || strcmp(fn, "-")==0) {
    f = stdout;
  } else {
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open cx target output file: %s\n", fn);
      return;
    }
  }

  fprintf(f, "t: %s\n", cxtgt.symbol);
  fprintf(f, "z: %g\n", cxtgt.z);
  fprintf(f, "m: %g\n", cxtgt.m);
  fprintf(f, "a: %g\n", cxtgt.a);
  fprintf(f, "b: %g\n", cxtgt.b);
  fprintf(f, "e: %g\n", cxtgt.e);
  fprintf(f, "x: %g\n", cxtgt.x);
  
  fflush(f);
  if (f != stdout) {
    fclose(f);
  }
}

void SetOptionNucleus(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "nucleus:xfermi")) {
    char buf[1024];
    strncpy(buf, s, 1023);
    int ns = StrSplit(buf, ',');
    int i;
    char *s = buf;
    for (i = 0; i < ns; i++) {
      while(*s == ' ' || *s == '\t') s++;
      if (i == 0) {
	_nfermi = atoi(s);
      } else if (i == 1) {
	_xfermi0 = atof(s);
      } else if (i == 2) {
	_xfermi1 = atof(s);
      }
      while (*s) s++;
      s++;
    }
    SetupFermi();
    return;
  }
  if (0 == strcmp(s, "nucleus:afermi")) {
    _afermi = dp;
    return;
  }
}

int IdxGround(int z, int k, int *c, int md) {
  int i;  
  i = (z*(z-1))/2 + k-1;
  if (md == -10) {
    *c = 0;
    return i;
  }
  if (md == -11) {
    *c = 1;
    return i;
  }
  if (md < 0) md = 2;
  double e = _ipot[i];
  if (e < 1e-6 || ((e-(int)e) < 0.001 && z-k >= md)) {
    *c = 1;
  } else {
    *c = 0;
  }
  return i;
}

int GetGround2J(int z, int k, int md) {
  int i, c;
  i = IdxGround(z, k, &c, md);

  if (c) {
    return _jlev1[i];
  } else {
    return _jlev[i];
  }
}

double GetGroundIP(int z, int k, int md) {
  int i, c;
  i = IdxGround(z, k, &c, md);

  if (c) {
    return _ipot1[i];
  } else {
    return _ipot[i];
  }
}

char *GetGroundLev(int z, int k, int md) {
  int i, c;
  i = IdxGround(z, k, &c, md);

  if (c) {
    return _glev1[i];
  } else {
    return _glev[i];
  }
}

int GetGroundParity(int z, int k, int md) {
  int i, c;
  i = IdxGround(z, k, &c, md);

  if (c) {
    return _plev1[i];
  } else {
    return _plev[i];
  }
}

char *GetGroundCfg(int z, int k, int md) {
  int i, c;
  i = IdxGround(z, k, &c, md);

  if (c) {
    return _gcfg1[i];
  } else {
    return _gcfg[i];
  }
}
