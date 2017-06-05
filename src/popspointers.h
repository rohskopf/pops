

#pragma once

#include "mpi.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"

#include "pops.h"

namespace POPS_NS
{
    class PopsPointers
    {
    public:
        PopsPointers(POPS *ptr) :
            pops(ptr),
            popsmemory(ptr->popsmemory),
            popsinput(ptr->popsinput),
            objfun(ptr->objfun),
            popstimer(ptr->popstimer),
            lmp(ptr->lmp)
            {}

        virtual ~PopsPointers() {}

    protected:
        POPS *pops;
        PopsMemory *&popsmemory;
        PopsInput *&popsinput;
        PopsTimer *&popstimer;
        Objfun *&objfun;
        LAMMPS_NS::LAMMPS *&lmp;
    };
}

