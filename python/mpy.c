/* Minimal main program -- everything is loaded from the library */

#include <Python.h>
#include "sysdef.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

extern DL_EXPORT(int) Py_Main();

int main(int argc, char *argv[]) {
#ifdef USE_MPI
  int rc;
  
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    fprintf(stderr, "MPI Initialization failed: error code %d\n",
	    rc);
    abort();
  }
#endif /* USE_MPI */

  Py_Main(argc, argv);
  
#ifdef USE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
