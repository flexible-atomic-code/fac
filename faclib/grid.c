#include "grid.h"

static char *rcsid="$Id: grid.c,v 1.3 2001/10/12 18:49:19 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

int AddPW(int *nkl0, double *kl, double *logkl, 
	  int maxkl, int n, int step) {
  int i;
  for (i = *nkl0; i < n+(*nkl0); i++) {
    if (i >= MAXNKL) {
      printf("Maximum partial wave grid points reached: "); 
      printf("%d > %d in constructing grid\n",  i, MAXNKL);
      abort();
    }
    kl[i] = kl[i-1] + step;
    logkl[i] = log(kl[i]);
    if ((int)(kl[i]) > maxkl) break;
  }
  (*nkl0) = i;
  return 0;
}

int SetPWGrid(int *nkl0, double *kl, double *logkl, 
	      int maxkl, int *ns, int *n, int *step) {
  int i, m, k, j;

  if ((*ns) > 0) {
    for (i = 0; i < (*ns); i++) {
      AddPW(nkl0, kl, logkl, maxkl, n[i], step[i]);
    }
    k = step[(*ns)-1]*2;
    j = 2;
  } else {
    (*ns) = -(*ns);
    if ((*ns) == 0) (*ns) = 8;
    AddPW(nkl0, kl, logkl, maxkl, (*ns), 1);
    k = 2;
    j = 2;
  }   

  m = kl[(*nkl0)-1];
  while (m+k <= maxkl) {
    AddPW(nkl0, kl, logkl, maxkl, j, k);
    m = kl[(*nkl0)-1];
    if (k < 50) k *= 2;
    else k += 50;
  }
  kl[(*nkl0)] = maxkl+1;
  return (*nkl0);
}

int SetTEGridDetail(double *te, double *logte, int n, double *x) {
  int i;
  
  if (n > MAXNTE) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  for (i = 0; i < n; i++) {
    te[i] = x[i];
    logte[i] = log(te[i]);
  }
  return n;
}

int SetTEGrid(double *te, double *logte, int n, double emin, double emax) {
  int i;
  double del;

  if (n < 1) {
    te[0] = -1.0;
    return 0;
  }

  if (emin < 0.0) {
    te[0] = emin;
    return n;
  }

  if (n > MAXNTE) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  if (n == 1) {
    te[0] = emin;
    if (logte) logte[0] = log(emin);
    return n;
  }

  if (n == 2) {
    te[0] = emin;
    te[1] = emax;
    if (logte) {
      logte[0] = log(emin);
      logte[1] = log(emax);
    }
    return n;
  }

  if (emax < emin) {
    printf("emin must > 0 and emax < emin in SetTEGrid\n");
    return -1;
  }
  
  del = emax - emin;
  del /= n-1.0;
  te[0] = emin;
  if (logte) logte[0] = log(emin);
  for (i = 1; i < n; i++) {
    te[i] = te[i-1] + del;
    if (logte) logte[i] = log(te[i]);
  }
  
  return n;
}
  
int SetEGridDetail(double *e, double *log_e, int n, double *xg) {
  int i;
  
  for (i = 0; i < n; i++) {
    e[i] = xg[i];
    log_e[i] = log(e[i]);
  }

  return n;
}

int SetEGrid(double *e, double *log_e, 
	     int n, double emin, double emax, double eth) {
  double del, et;
  int i;

  if (n < 1) {
    e[0] = -1.0;
    return 0;
  }
  if (emin < 0.0) {
    e[0] = emin;
    return 0;
  }

  if (emax < emin) {
    printf("emin must > 0 and emax < emin in SetEGrid\n");
    return -1;
  }

  et = fabs(eth);
  if (et > 1E-30) {
    emin += et;
    emax += et;
  }
  
  e[0] = emin;
  log_e[0] = log(emin);
  e[n-1] = emax;
  log_e[n-1] = log(emax);
  del = (log_e[n-1] - log_e[0])/(n-1.0);
  del = exp(del);
  for (i = 1; i < n-1; i++) {
    e[i] = e[i-1]*del;
    log_e[i] = log(e[i]);
  }

  if (eth > 1E-30) {
    for (i = 0; i < n; i++) {
      e[i] -= eth;
    }
  }
  return n;
}
