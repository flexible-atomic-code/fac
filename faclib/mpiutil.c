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

static int _initialized = 0;
static LOCK *_plock = NULL;
static volatile long long _cwid = -1;
static MPID mpi = {0, 1, 0};
#pragma omp threadprivate(mpi)

int SkipMPI() {
  int r = 0;
#if USE_MPI == 1
  if (mpi.nproc > 1) {
    if (mpi.wid%mpi.myrank != 0) {
      r = 1;
    } 
    mpi.wid++;
  }
  return r;
#elif USE_MPI == 2
#pragma omp critical  
  {
    if (mpi.nproc > 1) {
      mpi.wid++;      
      if (mpi.wid <= _cwid) r = 1;
      else _cwid = mpi.wid;
    }
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
      if (_plock) SetLock(_plock);
      printf("Rank=%d, ", myrank);
      vprintf(format, args);
      if (_plock) ReleaseLock(_plock);
    } else {
      if (myrank == ir%nproc) {	
	if (ir >= nproc) {
	  if (_plock) SetLock(_plock);
	  printf("Rank=%d, ", myrank);
	  vprintf(format, args);
	  if (_plock) SetLock(_plock);
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

MPID *DataMPI() {
  return &mpi;
}

int MPIReady() {
  return _initialized;
}

void InitializeMPI(int n) {
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
  if (n > 0) {
    int nm = omp_get_thread_limit();
    if (n > nm) {
      printf("OMP thread number exceeds limit: %d > %d\n", n, nm);
      n = nm;
    }
    omp_set_num_threads(n);
  } else if (n == 0) {
    omp_set_num_threads(1);
  }  
#pragma omp parallel
  {
    mpi.wid = 0;
    mpi.myrank = omp_get_thread_num();
    mpi.nproc = omp_get_num_threads();
  }
  _initialized = 1;
#endif
  _plock = (LOCK *) malloc(sizeof(LOCK));
  if (0 != InitLock(_plock)) {
    printf("cannot initialize lock in InitializeMPI\n");
    free(_plock);
    _plock = NULL;
  }
#if USE_MPI == 2
  CopyPotentialOMP(1);
#endif
#endif
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
  if (bf->nr > 1) {
    bf->buf = malloc(bf->nbuf);
  } else {
    bf->buf = NULL;
  }
#else
  bf->nr = 1;
  bf->mr = 0;
  bf->buf = NULL;
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
#else
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
