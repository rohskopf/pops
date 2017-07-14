#include <stdlib.h>
#include <iostream>
#include "pops.h"
#include "mpi.h"

using namespace POPS_NS;

int main(int argc, char **argv)
{

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  /* Begin a POPS instance */
  POPS *pops = new POPS(argc, argv);

  /* Delete the memory */
  delete pops;

  /* Close MPI */
  int MPI_Comm_free(MPI_Comm *comm);
  MPI_Finalize();

  return EXIT_SUCCESS;
}

