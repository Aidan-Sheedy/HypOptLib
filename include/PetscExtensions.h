
#pragma once

#include <petsc.h>
#include <vector>

class PetscExtensions
{
    public:
        /* Delete all constructors, all access should be static */
        PetscExtensions() = delete;
        PetscExtensions(const PetscExtensions&) = delete;

        static PetscErrorCode VecParallelFromSequential(Vec sequential, Vec *parallel);

        static PetscErrorCode VecSequentialFromParallel(Vec parallel, Vec *sequential);

        static PetscErrorCode VecGetOffProcessIndex(Vec parallelVector, PetscInt index, PetscScalar *result);

        static PetscErrorCode VecParallelFromStdVector(Vec parallel, std::vector<PetscScalar> stdVector);

        static PetscErrorCode VecLeftShift(Vec input, Vec *output);

        static PetscErrorCode PrintVecMinMaxMean(Vec vector, const char * name);

        static PetscErrorCode PrintVec(Vec vector, const char * name);

};
