

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
  class PopsInput: protected PopsPointers
  {
  public:
    PopsInput(class POPS *);
    ~PopsInput();
    void readinput();
    void readconfigs();

    /* readinput() variables */
    int gen_num_in;
    int extrainds_in;
    int pop_size_in;
    double mut_rate_in;
    double eliteperc_in; 
    int unknowns_in;
    int unknowns_extra_in;
    int M_in;
    int N_in;
    double w_f_in;
    double w_e_in;
    double w_e_s_in;
    double w_s_in;
    double w_p_in;
    double w_p_s_in;
    // Declare vector inputs
    vector<double> intervals_in;
    vector<int> ssl_in;
    vector<double> d_vec_in;
    vector<double> sym_vec_in;
    // Declare string inputs
    string sym_coeffs_in;
    string commands; // lammps commands
    string nontab_commands; // non-tabulated parameter commands

    /* readconfigs() group variables */
    int N_tot;
    int *N_vec;
    double *eio;
    double *eio_s;
    double *sio;
    double *pio;
    double *pio_s;
    double *x_lengths;
    double *y_lengths;
    double *z_lengths;
    double *xy_vec;
    double *xz_vec;
    double *yz_vec;

    /* readconfigs() per-atom variables */
    double *x_all;
    char *atomtypes;
    double *fio;
    
    

  };
}

