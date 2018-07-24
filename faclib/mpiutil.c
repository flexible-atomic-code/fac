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

#include "mpiutil.h"
#include "parser.h"
#include "stdarg.h"
#include "init.h"
#include "radial.h"
#include "cf77.h"

static int _initialized = 0;
static LOCK *_plock = NULL;
static LOCK *_mpilock = NULL;
static volatile long long _cwid = -1;
static MPID mpi = {0, 1, 0};
static double _tlock = 0, _tskip = 0;
static long long _nlock = 0;
#pragma omp threadprivate(mpi,_tlock,_tskip, _nlock)

int SkipWMPI(int w) {
  int r = 0;
#ifdef USE_MPI
  if (mpi.nproc > 1) {
    r = w%mpi.nproc != mpi.myrank;
  }
#endif
  return r;
}
  
int SkipMPI() {
  int r = 0;
#if USE_MPI == 1
  if (mpi.nproc > 1) {
    if (mpi.wid%mpi.nproc != mpi.myrank) {
      r = 1;
    } 
    mpi.wid++;
  }
  return r;
#elif USE_MPI == 2
  if (mpi.nproc > 1) {
#ifdef OMP_STAT
    double t0 = WallTime();
#endif    
    SetLockNT(_mpilock);
    mpi.wid++;      
    if (mpi.wid <= _cwid) r = 1;
    else _cwid = mpi.wid;
    ReleaseLock(_mpilock);
#ifdef OMP_STAT
    double t1 = WallTime();
    _tskip += t1-t0;
#endif
  }
  return r;
#else
  return 0;
#endif
}
  
void MPISeqBeg() {
#if USE_MPI == 1
  if (MPIReady()) {
    int myrank;
    int k;
    MPI_Status s;
    
    myrank = MPIRank(NULL);
    if (myrank > 0) {
      k = -1;
      MPI_Recv(&k, 1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD, &s);
      if (k != myrank-1) {
	printf("Error in MPISeqBeg %d %d\n", myrank, k);
      }
    }
  }
#endif
}

void MPISeqEnd() {
#if USE_MPI == 1
  if (MPIReady()) {
    int myrank;
    int nproc;
    
    myrank = MPIRank(&nproc);
    if (myrank < nproc-1) {
      MPI_Send(&myrank, 1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif
}

void MPrintf(int ir, char *format, ...) {
  va_list args;

  va_start(args, format);
#ifdef USE_MPI
  if (MPIReady()) {
    int myrank;
    int nproc;
    myrank = MPIRank(&nproc);
    if (ir < 0) {
      if (_plock) SetLockNT(_plock);
      printf("Rank=%d, ", myrank);
      vprintf(format, args);
      if (_plock) ReleaseLock(_plock);
    } else {
      if (myrank == ir%nproc) {	
	if (ir >= nproc) {
	  if (_plock) SetLockNT(_plock);
	  printf("Rank=%d, ", myrank);
	  vprintf(format, args);
	  if (_plock) SetLockNT(_plock);
	} else {
	  vprintf(format, args);
	}
      }
    }
  } else {
    vprintf(format, args);
  }
#else
  vprintf(format, args);
#endif
  va_end(args);
  fflush(stdout);
}

int MPIRank(int *np) {
  if (np) *np = mpi.nproc;
  return mpi.myrank;
}

int MyRankMPI() {
  return mpi.myrank;
}

int NProcMPI() {
  return mpi.nproc;
}

long long WidMPI() {
  return mpi.wid;
}

long long CWidMPI() {
  return _cwid;
}

void ResetWidMPI(void) {
#if USE_MPI == 2
  if (!MPIReady()) {
    printf("openmp not initialized\n");
    Abort(1);
  }
#endif
  _cwid = -1;  
#pragma omp parallel
  {
  mpi.wid = 0;
  }
}

void SetWidMPI(long long w) {
  mpi.wid = w;
}

MPID *DataMPI() {
  return &mpi;
}

int MPIReady() {
  return _initialized;
}

void SetLockMPI() {
  if (_mpilock) {
    SetLock(_mpilock);
  }
}

void ReleaseLockMPI() {
  if (_mpilock) {
    ReleaseLock(_mpilock);
  }
}

void InitializeMPI(int n, int m) {
#ifdef USE_MPI
  if (_initialized) {
    printf("MPI system already initialized\n");
    return;
  }
#if USE_MPI == 1
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.nproc);
  MPI_Initialized(&_initialized);
  if (!_initialized) {
    printf("cannot initialize MPI\n");
    exit(1);
  }
#elif USE_MPI == 2
  D1MACH(1);
  if (n > 0) {
    int nm = omp_get_thread_limit();
    if (n > nm) {
      printf("OMP thread number exceeds limit: %d > %d\n", n, nm);
      n = nm;
    }
    omp_set_num_threads(n);
    if (n > 1) {
      _mpilock = (LOCK *) malloc(sizeof(LOCK));
      if (0 != InitLock(_mpilock)) {
	printf("cannot initialize mpilock in InitializeMPI\n");
	free(_mpilock);
	_mpilock = NULL;
	Abort(1);
      }
    } else {
      _mpilock = NULL;
    }
  } else if (n == 0) {
    omp_set_num_threads(1);
  }  
#pragma omp parallel
  {
    mpi.wid = 0;
    mpi.myrank = omp_get_thread_num();
    mpi.nproc = omp_get_num_threads();
    _tlock = 0;
    _tskip = 0;
  }
  _initialized = 1;
  if (mpi.nproc == 1) {
    RemoveMultiLocks();
    if (m == 0) RemoveOrbitalLock();
  }
#endif
  _plock = (LOCK *) malloc(sizeof(LOCK));
  if (0 != InitLock(_plock)) {
    printf("cannot initialize lock in InitializeMPI\n");
    free(_plock);
    _plock = NULL;
  }
#if USE_MPI == 2
  if (m == 0) CopyPotentialOMP(1);
#endif
#endif
}

void SetLockWT(LOCK *x) {
  double t0 = WallTime();
  SetLockNT(x);
  double t1 = WallTime();
  _tlock += t1-t0;
  _nlock++;
}

double TimeSkip() {
  return _tskip;
}

double TimeLock() {
  return _tlock;
}

long long NumLock() {
  return _nlock;
}

void FinalizeMPI() {
#if USE_MPI == 1
  MPI_Finalize();
#endif
  if (_plock) {
    DestroyLock(_plock);
    free(_plock);
    _plock = NULL;
  }
  if (_mpilock) {
    DestroyLock(_mpilock);
    free(_mpilock);
    _mpilock = NULL;
  }
}

void Abort(int r) {
#if USE_MPI == 1
  MPI_Abort(MPI_COMM_WORLD, r);
#else  
  exit(r);
#endif
}

double WallTime() {
#if USE_MPI == 2
  return omp_get_wtime();
#else
  return ((double)clock())/CLOCKS_PER_SEC;
#endif
}

BFILE *BFileOpen(char *fn, char *md, int nb) {  
  BFILE *bf;  
  bf = malloc(sizeof(BFILE));  
  bf->p = 0;
  bf->w = &bf->p;
  bf->n = 0;
  bf->eof = 0;
  bf->nbuf = nb>=0?nb:RBUFL;
  if (bf->nbuf == 0) {
    bf->nr = 1;
    bf->mr = 0;
    bf->buf = NULL;
    bf->f = fopen(fn, md);
    if (bf->f == NULL) {
      free(bf);
      return NULL;
    }
    bf->fn = NULL;
    return bf;
  }
#if USE_MPI == 1
  bf->mr = MPIRank(&bf->nr);
  if (bf->mr == 0) {
    bf->f = fopen(fn, md);
    if (bf->f == NULL) bf->p = -1;
  } else {
    bf->f = NULL;
  }
  if (bf->nr > 1) {
    MPI_Bcast(&bf->p, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (bf->p < 0) {
    free(bf);
    return NULL;
  }
  if (bf->nbuf > 0) {
    bf->buf = malloc(bf->nbuf);
  } else {
    bf->buf = NULL;
  }
#elif USE_MPI == 2
  bf->mr = MPIRank(&bf->nr);
  bf->f = fopen(fn, md);
  if (bf->f == NULL) {
    free(bf);
    //printf("cannot open file: %s\n", fn);
    return NULL;
  }
  bf->p = 0;
  bf->n = 0;
  if (bf->nbuf > 0) {
    bf->buf = malloc(bf->nbuf*bf->nr);
    bf->w = malloc(bf->nr*sizeof(int));
    int i;
    for (i = 0; i < bf->nr; i++) bf->w[i] = 0;
    InitLock(&bf->lock);
  } else {
    bf->buf = NULL;
  }
#else
  bf->nr = 1;
  bf->mr = 0;
  if (bf->nbuf > 0) {
    bf->buf = malloc(bf->nbuf);
  } else {
    bf->buf = NULL;
  }
  bf->f = fopen(fn, md);
  if (bf->f == NULL) {
    free(bf);
    return NULL;
  }
#endif

  bf->fn = malloc(strlen(fn)+1);
  strcpy(bf->fn, fn);

  return bf;
}

int BFileClose(BFILE *bf) {
  int r = 0;
#if USE_MPI == 1
  if (bf == NULL) return 0;
  if (bf->nr <= 1) {
    r = fclose(bf->f);
  } else {
    if (bf->mr == 0) {
      r = fclose(bf->f);
    }  
    free(bf->buf);
  }
#elif USE_MPI == 2
  if (bf == NULL) return 0;
  BFileFlush(bf);
  if (bf->nr <= 1) {
    r = fclose(bf->f);
  } else {
    r = fclose(bf->f);
    if (bf->nbuf > 0) {
      free(bf->buf);
      bf->buf = NULL;
      free(bf->w);
      bf->w = NULL;
      DestroyLock(&bf->lock);
    }
  }
#else
  if (bf == NULL) return 0;
  BFileFlush(bf);
  if (bf->nbuf > 0) {
    free(bf->buf);
    bf->buf = NULL;
  }
  r = fclose(bf->f);
#endif
  
  free(bf->fn);
  free(bf);

  return r;
}

size_t BFileRead(void *ptr, size_t size, size_t nmemb, BFILE *bf) {
#if USE_MPI == 1
  if (bf->nr <= 1) {
    return fread(ptr, size, nmemb, bf->f);
  }
  if (size > bf->nbuf) {
    if (bf->mr == 0) {
      printf("buffer size %d smaller than data size %d\n", (int)bf->nbuf, (int)size);
    }
    Abort(1);
  }
  int nb = bf->n - bf->p;
  int n = nb/size;
  int nr=0, nm=0, nn=0, nread;  
  nread = 0;
  while (nmemb) {
    if (n >= nmemb) {
      nr = size*nmemb;
      memcpy(ptr, bf->buf+bf->p, nr);
      bf->p += nr;
      nread += nmemb;
      return nread;
    } else if (n > 0) {
      nm = size*n;
      memcpy(ptr, bf->buf+bf->p, nm);
      ptr += nm;
      bf->p += nm;
      nb -= nm;
      nread += n;
      nmemb -= n;
    }
    if (bf->eof) break;    
    if (nb > 0) {
      memmove(bf->buf, bf->buf+bf->p, nb);    
    }
    bf->p = 0;
    bf->n = nb;
    if (bf->mr == 0) {
      nn = bf->nbuf - bf->n;
      nr = fread(bf->buf+bf->n, 1, nn, bf->f);
      if (nr < nn) {
	bf->eof = 1;
      }
      bf->n += nr;
    }
    MPI_Bcast(&bf->n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (bf->n > nb) {
      MPI_Bcast(bf->buf+nb, bf->n-nb, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&bf->eof, 1, MPI_INT, 0, MPI_COMM_WORLD);
    nb = bf->n - bf->p;
    n = nb/size;
  }
  return nread;
#else
  return fread(ptr, size, nmemb, bf->f);
#endif
}

char *BFileGetLine(char *s, int size1, BFILE *bf) {
#if USE_MPI == 1
  if (bf->nr <= 1) {
    return fgets(s, size1, bf->f);
  }
  int n;
  int size = size1-1;

  n = BFileRead(s, size, 1, bf);
  if (n == 1) {
    bf->p -= size;
  } else {
    size = bf->n-bf->p;
    if (size == 0) return NULL;
    memcpy(s, bf->buf+bf->p, size);
  }
  int i;
  for (i = 0; i < size; i++, bf->p++) {
    if (s[i] == '\n') {
      bf->p++;
      i++;
      break;
    }
  }
  if (bf->p == bf->nbuf) {
    bf->p = 0;
    bf->n = 0;
  }
  s[i] = '\0';
  return s;
#else
  return fgets(s, size1, bf->f);
#endif
}

void BFileRewind(BFILE *bf) {
#if USE_MPI == 1
  if (bf->nr <= 1) {
    rewind(bf->f);
    return;
  }
  bf->p = 0;
  bf->n = 0;
  bf->eof = 0;
  if (bf->mr == 0) rewind(bf->f);
#else
  rewind(bf->f);
#endif
}

size_t BFileWrite(void *ptr, size_t size, size_t nmemb, BFILE *bf) {
  int n, m, k, mr;
  char *buf;
  
  if (bf->buf == NULL) {
    n = fwrite(ptr, size, nmemb, bf->f);
    return n;
  }
  
#if USE_MPI == 2
  mr = MPIRank(NULL);
#else
  mr = 0;
#endif

  buf = bf->buf + bf->nbuf*mr;
  m = size*nmemb;
  k = bf->nbuf - bf->w[mr];
  n = 0;
  
  if (m >= k) {
#if USE_MPI == 2
    if (bf->nr > 1) SetLock(&bf->lock);
#endif
    if (bf->w[mr] > 0) {
      n = fwrite(buf, 1, bf->w[mr], bf->f);
      bf->w[mr] = 0;
    }
    n = fwrite(ptr, size, nmemb, bf->f);
#if USE_MPI == 2
    if (bf->nr > 1) ReleaseLock(&bf->lock);
#endif
  } else {
    memcpy(buf+bf->w[mr], ptr, m);
    bf->w[mr] += m;
    n = nmemb;
  }
  return n;
}

int BFileSeek(BFILE *bf, long offset, int w) {
  BFileFlush(bf);
  return fseek(bf->f, offset, w);
}

long BFileTell(BFILE *bf) {
  BFileFlush(bf);
  return ftell(bf->f);
}

int BFileFlush(BFILE *bf) {
  int i;
  if (bf->buf != NULL) {
    for (i = 0; i < bf->nr; i++) {
      if (bf->w[i] > 0) {
	fwrite(bf->buf+i*bf->nbuf, 1, bf->w[i], bf->f);
	bf->w[i] = 0;
      }
    }
  }
  return fflush(bf->f);
}

int CompareRandIdx(const void *p1, const void *p2) {
  RANDIDX *r1, *r2;
  r1 = (RANDIDX *) p1;
  r2 = (RANDIDX *) p2;
  if (r1->r < r2->r) return -1;
  else if (r1->r > r2->r) return 1;
  return 0;
}

int ComparePtrIdx(const void *p1, const void *p2) {
  PTRIDX *r1, *r2;
  r1 = (PTRIDX *) p1;
  r2 = (PTRIDX *) p2;
  if (r1->r < r2->r) return -1;
  else if (r1->r > r2->r) return 1;
  return 0;
}

RANDIDX *RandList(int n) {
  int i;
  RANDIDX *w;
  const long seed = 0xFFFF;
  w = (RANDIDX *) malloc(sizeof(RANDIDX)*n);
  srand48(seed);
  for (i = 0; i < n; i++) {
    w[i].r = drand48();
    w[i].i = i;
  }
  qsort(w, n, sizeof(RANDIDX), CompareRandIdx);
  return w;
}

void RandIntList(int n, int *k) {
  RANDIDX *w = RandList(n);
  int i;
  for (i = 0; i < n; i++) {
    w[i].i = k[w[i].i];
  }
  for (i = 0; i < n; i++) {
    k[i] = w[i].i;    
  }
  free(w);
}

void ArgSort(int n, double *r, int *k) {
  RANDIDX *rid;
  int i;
  rid = malloc(sizeof(RANDIDX)*n);
  for (i = 0; i < n; i++) {
    rid[i].i = i;
    rid[i].r = r[i];
  }
  qsort(rid, n, sizeof(RANDIDX), CompareRandIdx);
  for (i = 0; i < n; i++) {
    k[i] = rid[i].i;
  }
  free(rid);
}

void PtrSort(int n, char **r, int *k) {
  PTRIDX *rid;
  int i;
  rid = malloc(sizeof(PTRIDX)*n);
  for (i = 0; i < n; i++) {
    rid[i].i = i;
    rid[i].r = r[i];
  }
  qsort(rid, n, sizeof(PTRIDX), ComparePtrIdx);
  for (i = 0; i < n; i++) {
    k[i] = rid[i].i;
  }
  free(rid);
}
