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

/* Minimal main program -- everything is loaded from the library */

#include <Python.h>
#include "sysdef.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

extern DL_EXPORT(int) Py_Main();

int main(int argc, char *argv[]) {
#if USE_MPI == 1
  int rc;
  
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    fprintf(stderr, "MPI Initialization failed: error code %d\n",
	    rc);
    abort();
  }
#endif /* USE_MPI */

  Py_Main(argc, argv);
  
#if USE_MPI == 1
  MPI_Finalize();
#endif
  
  return 0;
}
