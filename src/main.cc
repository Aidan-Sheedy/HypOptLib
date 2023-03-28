#include "Filter.h"
#include "LinearElasticity.h"
#include "MMA.h"
#include "MPIIO.h"
#include "TopOpt.h"
#include "mpi.h"
#include <petsc.h>

#include <random>

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
    // PetscInt itr = 0;
    // opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !
    // mma->SetAsymptotes(0.2, 0.65, 1.05);

    // STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
    // ierr = filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);
    // CHKERRQ(ierr);

    // ierr = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), &(opt->gx[0]), opt->dfdx, opt->dgdx[0],
    //                                                         opt->xPhys, opt->Emin, opt->Emax, opt->penal,
    //                                                         opt->volfrac);
    // CHKERRQ(ierr);

    // STEP 7: OPTIMIZATION LOOP
    PetscPrintf(PETSC_COMM_WORLD, "\n\n######################## Hyperoptimization ########################\n");

    /* Setup random starting positions */

    // Vec initialPositions;

    // PetscCall(VecCreateSeq(PETSC_COMM_WORLD, opt->n, &initialPositions));
    // // PetscCall(VecSetSizes(this->evenNoseHooverMass, PETSC_DECIDE, NHChainOrder/2));
    // PetscCall(VecSetFromOptions(this->evenNoseHooverMass));
    // PetscCall(VecDuplicate(this->evenNoseHooverMass, &(this->oddNoseHooverMass)));
// #if 0
    PetscScalar *initialValues;
    PetscCall(VecGetArray(opt->x, &initialValues));

    PetscInt localSize;
    PetscCall(VecGetLocalSize(opt->x, &localSize));

    std::default_random_engine generator;
    std::uniform_real_distribution<PetscScalar> distribution;

    if (0.5 < opt->volfrac)
    {
        distribution = std::uniform_real_distribution<PetscScalar>(2 * opt->volfrac - 1, 1);
    }
    else
    {
        distribution = std::uniform_real_distribution<PetscScalar>(0, 2 * opt->volfrac);
    }

    for (PetscInt i = 0; i < localSize; i++)
    {
        initialValues[i] = distribution(generator); //opt->volfrac;
    }

    PetscCall(VecRestoreArray(opt->x, &initialValues));
// #endif
    LagrangianMultiplier lagmult(filter, opt);

    /** @todo Figure out how to pass stuff in from the wrapper! */
    PetscScalar temperature = 0.1;//0;//5;
    PetscScalar NHChainOrder = 10;
    PetscScalar dt = 0.002;//0.01;//0.001;//0.0002;

    Hyperoptimization solver;
    PetscCall(solver.init(physics,
                opt,
                filter,
                // data,
                lagmult,
                temperature,
                opt->x, /** @todo initialize the positions properly */
                NHChainOrder,
                opt->maxItr,
                dt,
                2000,
                false));

    PetscPrintf(PETSC_COMM_WORLD, "Initialized, starting design loop\n");

    double t1 = MPI_Wtime();
    solver.runDesignLoop();
    double t2 = MPI_Wtime();

    PetscPrintf(PETSC_COMM_WORLD, "Done loop!\nTotal Runtime: %f\nAverage Iteration Time: %f\n***Note these timings include file saving***", t2 - t1, (t2-t1)/opt->maxItr);
    // PetscPrintf(PETSC_COMM_WORLD, "Done loop!\nT2 - T1: %f\nTime 1: %f\nTime 2: %f\n***Note these timings include file saving***", t2-t1, t1, t2);

    PetscPrintf(PETSC_COMM_WORLD, "\n\n###################################################################\n");


    // // Write restart WriteRestartFiles
    // opt->WriteRestartFiles(&itr, mma);  
    // physics->WriteRestartFiles();

    // Dump final design
    // output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhys, itr + 1);

    // STEP 7: CLEAN UP AFTER YOURSELF
    // delete mma;
    // delete output;
    delete filter;
    delete opt;
    delete physics;

    // Finalize PETSc / MPI
    PetscFinalize();
    return 0;
}
