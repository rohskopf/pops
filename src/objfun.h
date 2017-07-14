

#pragma once

#include <vector>
#include <string>
#include "mpi.h"

#include <iostream>
#include <new>
#include <cstdlib>
#include "popspointers.h"

using namespace std;

namespace POPS_NS
{
  class Objfun: protected PopsPointers
  {
  public:
    Objfun(class POPS *);
    ~Objfun();

    vector<double> calcfitness(double sub_xpop_r[], int nepp);
    double *params;

    double calcz(double params[], int ind);
    double *x;
    double *f;
    double *e; // energy pointer for compute extraction
    double *s; // stress vector for stress tensor extraction
    double *press; // pressure vector for pressure extraction
    double *fi; // force vector for storage from *f
    double *ei; // energies that will store the "e" scalars
    double *ei_s; // energies that will store the "e" scalars
    double *pi; // energies that will store the "e" scalars
    double *pi_s; // energies that will store the "e" scalars

    vector<double> sub_Zf2_vec;
    vector<double> sub_Ze2_vec;
    vector<double> sub_Ze2s_vec;
    vector<double> sub_Zs2_vec;
    vector<double> sub_Z2fmax_vec;
    
    

  };
}

