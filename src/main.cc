#include "Filter.h"
#include "LinearElasticity.h"
#include "MMA.h"
#include "MPIIO.h"
#include "TopOpt.h"
#include "mpi.h"
#include <petsc.h>

#include "Hyperoptimization.h"

/*
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

Updated: June 2019, Niels Aage
Copyright (C) 2013-2019,

Disclaimer:
The authors reserves all rights but does not guaranty that the code is
free from errors. Furthermore, we shall not be liable in any event
caused by the use of the program.
 */

static char help[] = "3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n";

int main(int argc, char* argv[]) {

    // Error code for debugging
    PetscErrorCode ierr;

    // Initialize PETSc / MPI and pass input arguments to PETSc
    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
    TopOpt* opt = new TopOpt();

    // STEP 2: THE PHYSICS
    LinearElasticity* physics = new LinearElasticity(opt->da_nodes);

    // STEP 3: THE FILTERING
    Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);

    // STEP 4: VISUALIZATION USING VTK
    MPIIO* output = new MPIIO(opt->da_nodes, 3, "ux, uy, uz", 3, "x, xTilde, xPhys");
    // STEP 5: THE OPTIMIZER MMA
    // MMA*     mma;
    PetscInt itr = 0;
    // opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !
    // mma->SetAsymptotes(0.2, 0.65, 1.05);

    // STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
    ierr = filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);
    CHKERRQ(ierr);

    // STEP 7: OPTIMIZATION LOOP
    PetscPrintf(PETSC_COMM_WORLD, "\n\n######################## Hyperoptimization ########################\n");
    LagrangianMultiplier lagmult(filter, opt);

    /** @todo Figure out how to pass stuff in from the wrapper! */
    PetscScalar temperature = 10;
    PetscScalar NHChainOrder = 10;
    PetscFloat dt = 0.001;

    Hyperoptimization solver;
    solver.init(physics,
                opt,
                filter,
                // data,
                lagmult,
                temperature,
                opt->xPhys, /** @todo initialize the positions properly */
                NHChainOrder,
                opt->maxItr,
                dt);

    PetscPrintf(PETSC_COMM_WORLD, "Initialized, starting loop\n");

    solver.runDesignLoop();

    PetscPrintf(PETSC_COMM_WORLD, "Done loop!\n");

    PetscPrintf(PETSC_COMM_WORLD, "\n\n###################################################################\n");


    // // Write restart WriteRestartFiles
    // opt->WriteRestartFiles(&itr, mma);
    // physics->WriteRestartFiles();

    // Dump final design
    // output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhys, itr + 1);

    // STEP 7: CLEAN UP AFTER YOURSELF
    // delete mma;
    delete output;
    delete filter;
    delete opt;
    delete physics;

    // Finalize PETSc / MPI
    PetscFinalize();
    return 0;
}
