#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// driver.cpp
//
//===========================================================================
#include <iostream>
#include <fstream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"

#include "AmrHydro.H"

#include "HydroIBC.H"
#include "memusage.H" 

#ifdef CH_USE_PETSC
#include "petsc.h"
#endif

//===========================================================================
// 2D Amr driver example
//
//===========================================================================
int
main(int argc, char* argv[])
{
    int ierr = 0;

#ifdef CH_USE_PETSC
    ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    CHKERRQ(ierr);
#else
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif
#endif // end petsc conditional

    { // Begin nested scope

#ifdef CH_MPI
        MPI_Barrier(Chombo_MPI::comm);
#endif
        int rank, number_procs;
#ifdef CH_MPI
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
        MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
        rank = 0;
        number_procs = 1;
#endif

        if (argc < 2)
        {
            std::cerr << " usage: " << argv[0] << " <input_file>\n";
            exit(0);
        }
        char* in_file = argv[1];
        ParmParse pp(argc - 2, argv + 2, NULL, in_file);
        ParmParse ppMain("main");

        std::string poutBaseName = "pout";
        ppMain.query("poutBaseName", poutBaseName);
        setPoutBaseName(poutBaseName);

        RealVect domainSize;
        Vector<Real> domSize(SpaceDim);
        ppMain.getarr("domain_size", domSize, 0, SpaceDim);
        domainSize = RealVect(D_DECL(domSize[0], domSize[1], domSize[2]));

        AmrHydro amrObject;

        // ---------------------------------------------
        // set IBC -- this includes initial conditon
        // and boundary conditions
        // ---------------------------------------------
        HydroIBC* thisIBC = new HydroIBC();

        amrObject.setIBC(thisIBC);

        amrObject.setDomainSize(domainSize);

        // set up initial grids, initialize data, etc.
        amrObject.initialize();

        int maxStep;
        Real maxTime;
        // Real startTime;
        ppMain.get("maxTime", maxTime);
        ppMain.get("maxStep", maxStep);

        amrObject.run(maxTime, maxStep);

        if (thisIBC != NULL)
        {
            delete thisIBC;
            thisIBC = NULL;
        }

    } // end nested scope

    CH_TIMER_REPORT();
    dumpmemoryatexit();

#ifdef CH_USE_PETSC
    ierr = PetscFinalize();
    CHKERRQ(ierr);
#else
#ifdef CH_MPI
    MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional

    return ierr;
}
