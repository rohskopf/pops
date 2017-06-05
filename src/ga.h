

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
  class Ga: protected PopsPointers
  {
  public:
    Ga(class POPS *);
    ~Ga();
    
    void initialize();
    char **pop;
    double seed; 
    int unknowns_tot;
    vector<double> Z_vec;
    vector<int> ssl;
    vector<double> d_vec; 
    vector<double> intervals;
    int length; 

    double* phenotype();
    double **xpop; 
    double *xpop_r; 
    double *subxpop;

    void selection();
    double **sorted_xpop;
    char **sorted_pop;
    char **parents;
    char **rest;
    char **elite;

    void crossover();
    int *parent1;
    int *parent2;
    int *child1;
    int *child2;
    char **children;

    void mutation();

    FILE * history = fopen("HISTORY" ,"w");
    FILE * individuals = fopen("INDIVIDUALS", "w");
    

  };
}

