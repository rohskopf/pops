
#pragma once

#include <string>
#include "popspointers.h"

#if defined(WIN32) || defined(_WIN32)
#include <Windows.h>
#else
#include <time.h>
#include <sys/time.h>
#endif

namespace POPS_NS
{
    class PopsTimer : protected PopsPointers
    {
    public:
        PopsTimer(class POPS *);
        ~PopsTimer();

        void reset();
        double elapsed();
        void print_elapsed();
        std::string DateAndTime();

    private:
#if defined(WIN32) || defined(_WIN32)
        LARGE_INTEGER time_ref;
        LARGE_INTEGER frequency;
#else
        timeval time_ref;
#endif
    };
}

