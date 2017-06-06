

#include <iostream>
#include <iomanip>
#include <vector>
#include "popsmemory.h"
#include "popsinput.h"
#include "objfun.h"
#include "ga.h"
#include "popstimer.h"
#include "mpi.h"

using namespace std;

using namespace POPS_NS;

POPS::POPS(int narg, char **arg)
{

  /************************** Set up MPI settings **************************/

  int color,key,global,local;
  MPI_Comm comm;

  // Split the communicators so that multiple instances of LAMMPS can be run
  MPI_Comm_rank(MPI_COMM_WORLD, &global);
  color = global / 1; // Change "1" to 2 in order to use 2 procs per instance, etc..
  key = global; 
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
  MPI_Comm_rank(comm,&local);

  //  Get the number of processes.
  procs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group
  //  Get the individual process ID.
  rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  nepp = atoi(arg[2]); // Num elements per proc

  /************************** Initial Screen Output **************************/
  if (rank == 0)
  {
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << " +                             POPS 0.0                            +" << std::endl;
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << " Running on " << procs << " procs" << std::endl;
  }
  
  popstimer = new PopsTimer(this);
  if (rank == 0) std::cout << " Job started at " << popstimer->DateAndTime() << std::endl;

  /************************** Create LAMMPS Pointer **************************/
  char *args1[] = {
    (char *) "lmp",
    (char *) "-screen",
    (char*) "none",
  0};

  lmp = new LAMMPS_NS::LAMMPS(3,args1,comm);

  /************************** Extract Optional Random Number Seed **************************/
  bool time_seed;
  /* Extract optional random number seed */
  seed;
  if (narg == 3) {
    char *seed_char = arg[1]; // obtain a value from the input to seed the random number generator
    seed = atoi(seed_char); // convert the seed value to an int to be used for srand()
    time_seed = false;
  }
  else {
    seed = 0;
    time_seed = true;
  }

  if (time_seed)
  {
    // Calculate current number of seconds since the new year
    time_t now;
    struct tm newyear;
    double seconds;
    time(&now);  // get current time; same as: now = time(NULL)
    newyear = *localtime(&now);
    newyear.tm_hour = 0; newyear.tm_min = 0; newyear.tm_sec = 0;
    newyear.tm_mon = 0;  newyear.tm_mday = 1;
    seconds = difftime(now,mktime(&newyear));
    // Add this number of seconds to the seed
    seed = seed + seconds;
  }

  /************************** Proceed with POPS **************************/

  // Dynamically allocated pointers
  create();

  // Grab user input to initialize settings
  initialize();

  // Run the GA
  int gen_num = popsinput->gen_num_in;
  g = 1;
  //if (rank == 0) 
  //{
    ga->initialize();;
  //}
  double *xpop_r = NULL;
  popsmemory->allocate(sub_xpop_r, nepp);
  for (int i=0; i<nepp; i++)
  {
    sub_xpop_r[i] = 0.0;
  }
  while (g <= gen_num)
  {

  //double *xpop_r = NULL;
  //double *sub_xpop_r = (double *)malloc(sizeof(double) * nepp);
  //if (rank == 0) 
  //{
    //double *xpop_r = NULL;
    xpop_r = ga->phenotype();
    //std::cout << xpop_r[0] << std::endl;
    //MPI_Scatter(xpop_r, nepp, MPI::DOUBLE, sub_xpop_r, nepp, MPI::DOUBLE, 0, MPI_COMM_WORLD);
  //}

  //double *sub_xpop_r = (double *)malloc(sizeof(double) * nepp);
    
  // Scatter the random numbers from the root process to all processes in the MPI world
  //cout << "Scatter parameters" << endl; 
  MPI_Scatter(xpop_r, nepp, MPI::DOUBLE, sub_xpop_r, nepp, MPI::DOUBLE, 0, MPI_COMM_WORLD);

  // Get the sub error quantities per proc
  vector<double> sub_Z_vec = objfun->calcfitness(sub_xpop_r, nepp);
  vector<double> sub_Zf2_vec = objfun->sub_Zf2_vec;

  // Convert vectors to arrays for MPI_Gather
  int sub_Z_vec_size = sub_Z_vec.size();
  double sub_Z_vec_arr[sub_Z_vec_size];
  copy(sub_Z_vec.begin(), sub_Z_vec.end(), sub_Z_vec_arr);
  int sub_Zf2_vec_size = sub_Zf2_vec.size();
  double sub_Zf2_vec_arr[sub_Zf2_vec_size];
  copy(sub_Zf2_vec.begin(), sub_Zf2_vec.end(), sub_Zf2_vec_arr);

  // Allocate memory for Z_vec_arr
  popsmemory->allocate(Z_vec_arr, popsinput->pop_size_in); 
  popsmemory->allocate(Zf2_vec_arr, popsinput->pop_size_in);

  // Gather the sub_Z_vec vectors into Z_vec
  //cout << "Gather parameters" << endl;
  MPI_Gather(sub_Z_vec_arr,popsinput->pop_size_in/procs, MPI::DOUBLE, Z_vec_arr, popsinput->pop_size_in/procs, MPI::DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(sub_Zf2_vec_arr,popsinput->pop_size_in/procs, MPI::DOUBLE, Zf2_vec_arr, popsinput->pop_size_in/procs, MPI::DOUBLE, 0, MPI_COMM_WORLD);
  //cout << "Parameters gathered" << endl;
  // Resize Z_vec to be same size as Z_vec_arr
  Z_vec.resize(popsinput->pop_size_in);
  Zf2_vec.resize(popsinput->pop_size_in);
  // Convert Z_vec_arr to Z_vec
  if (rank == 0)
  {
    for (int i=0; i<popsinput->pop_size_in; i++)
    {
      //cout << Z_vec_arr[i] << endl;
      Z_vec[i] = Z_vec_arr[i];
      Zf2_vec[i] = Zf2_vec_arr[i];
    }

    // Run the rest of the GA
    ga->selection();
    ga->crossover();
    ga->mutation();
  }

  // Deallocate memory for Z_vec_arr
  //cout << "Deallocate Z_vec_arr" << endl;
  popsmemory->deallocate(Z_vec_arr);
  popsmemory->deallocate(Zf2_vec_arr);

  g++;
  }
  popsmemory->deallocate(sub_xpop_r);

  // Delete dynamically allocated pointers
  finalize();

  if (rank == 0) std::cout << std::endl << " Job finished at " 
        << popstimer->DateAndTime() << std::endl;
}

void POPS::create()
{
    popsmemory = new PopsMemory(this);
    popsinput = new PopsInput(this);
    objfun = new Objfun(this);
    ga = new Ga(this);
}

void POPS::initialize()
{
  popsinput->readinput();
  popsinput->readconfigs();
}

POPS::~POPS()
{
  delete lmp;
  delete popstimer;

}

void POPS::finalize()
{

    delete popsmemory;
    delete popsinput;
    delete objfun;
    delete ga;
}

