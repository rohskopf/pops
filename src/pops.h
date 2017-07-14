
/* Declaration of pointers used in the whole program. */

#pragma once
#include <string>
#include <vector>
#include "mpi.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"


namespace POPS_NS
{
    class POPS
    {
    public:


        LAMMPS_NS::LAMMPS *lmp;
        class PopsMemory *popsmemory;
        class PopsInput *popsinput;
        class PopsTimer *popstimer;
        class Objfun *objfun;
        class Ga *ga;
        POPS(int, char **);
        ~POPS();
        void create();
        void initialize();
        void finalize();
        int procs;
        int rank;
        int nepp; // num elements per proc

        int seed;

        double *Z_vec_arr;
        std::vector<double> Z_vec;
        double *Zf2_vec_arr;
        std::vector<double> Zf2_vec;
        double *Ze2_vec_arr;
        std::vector<double> Ze2_vec;
        double *Ze2s_vec_arr;
        std::vector<double> Ze2s_vec;
        double *Zs2_vec_arr;
        std::vector<double> Zs2_vec;
        double *Z2fmax_vec_arr;
        std::vector<double> Z2fmax_vec;

        int g; //generation number
        int gen_num;
        double *sub_xpop_r;

        std::string mode;
    };
}

