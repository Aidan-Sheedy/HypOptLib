/***************************************************************************//**
 * @file PetscExtensions.cc
 *
 * Implementation of useful functions not included in the standard Petsc
 * library.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
 ******************************************************************************/

#include "PetscExtensions.h"

PetscErrorCode PetscExtensions::VecParallelFromSequential(Vec sequential, Vec *parallel)
{
    PetscErrorCode errorStatus = 0;

    PetscInt sequentialSize;
    PetscCall(VecGetSize(sequential, &sequentialSize));

    VecScatter  scatter;
    IS          from;
    IS          to;
    PetscInt    idx_from[sequentialSize];
    PetscInt    idx_to[sequentialSize];

    for (PetscInt i = 0; i < sequentialSize; i++)
    {
        idx_from[i] = i;
        idx_to[i]   = i;
    }

    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, sequentialSize, idx_from, PETSC_COPY_VALUES, &from));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, sequentialSize, idx_to,   PETSC_COPY_VALUES, &to));

    PetscCall(VecScatterCreate(sequential, from,   *parallel, to, &scatter));
    PetscCall(VecScatterBegin(scatter, sequential, *parallel, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter,   sequential, *parallel, INSERT_VALUES, SCATTER_FORWARD));

    /* Cleanup */
    PetscCall(ISDestroy(&from));
    PetscCall(ISDestroy(&to));
    PetscCall(VecScatterDestroy(&scatter));

    return errorStatus;
}

PetscErrorCode PetscExtensions::VecSequentialFromParallel(Vec parallel, Vec *sequential)
{
    PetscErrorCode errorStatus = 0;

    PetscInt parallelSize;
    PetscCall(VecGetSize(parallel, &parallelSize));

    Vec localSequential;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, parallelSize, &localSequential));

    VecScatter  scatter;      /* scatter context */
    IS          from;
    IS          to;     /* index sets that define the scatter */
    PetscInt    idx_from[parallelSize];
    PetscInt    idx_to[parallelSize];

    for (PetscInt i = 0; i < parallelSize; i++)
    {
        idx_from[i] = i;
        idx_to[i]   = i;
    }

    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, parallelSize, idx_from, PETSC_COPY_VALUES, &from));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, parallelSize, idx_to,   PETSC_COPY_VALUES, &to));

    PetscCall(VecScatterCreate(parallel, from, localSequential, to, &scatter));
    PetscCall(VecScatterBegin(scatter, parallel, localSequential, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter,   parallel, localSequential, INSERT_VALUES, SCATTER_FORWARD));

    /* Cleanup */
    PetscCall(ISDestroy(&from));
    PetscCall(ISDestroy(&to));
    PetscCall(VecScatterDestroy(&scatter));

    PetscCall(VecDuplicate(localSequential, sequential));
    PetscCall(VecCopy(localSequential, *sequential));

    PetscCall(VecDestroy(&localSequential));

    return errorStatus;
}

PetscErrorCode PetscExtensions::VecGetOffProcessIndex(Vec parallelVector, PetscInt index, PetscScalar *result)
{
    PetscErrorCode errorStatus = 0;

    Vec         sequentialDestination; /* destination vector */
    VecScatter  scatter;      /* scatter context */
    IS          from, to;     /* index sets that define the scatter */
    PetscScalar *values;
    PetscInt    idx_to = 0;

    /* Create sequential vector */
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, 1, &sequentialDestination));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, 1, &index, PETSC_COPY_VALUES, &from));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, 1, &idx_to, PETSC_COPY_VALUES, &to));

    /* Perform scatter */
    PetscCall(VecScatterCreate(parallelVector, from, sequentialDestination, to, &scatter));
    PetscCall(VecScatterBegin(scatter, parallelVector, sequentialDestination, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, parallelVector, sequentialDestination, INSERT_VALUES, SCATTER_FORWARD));

    /* Cleanup */
    PetscCall(ISDestroy(&from));
    PetscCall(ISDestroy(&to));
    PetscCall(VecScatterDestroy(&scatter));

    /* Save result */
    PetscCall(VecGetArray(sequentialDestination, &values));
    *result = values[0];
    PetscCall(VecRestoreArray(sequentialDestination, &values));

    PetscCall(VecDestroy(&sequentialDestination));

    return errorStatus;
}

PetscErrorCode PetscExtensions::VecParallelFromStdVector(std::vector<PetscScalar> stdVector, Vec parallel)
{
    PetscErrorCode errorStatus = 0;

    PetscInt parallelGlobalSize = 0;
    PetscCall(VecGetSize(parallel, &parallelGlobalSize));

    if (parallelGlobalSize != stdVector.size())
    {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscScalar *parallelArray;
    PetscInt localParallelLowerIndex;
    PetscInt localParallelUpperIndex;
    PetscCall(VecGetArray(parallel, &parallelArray));
    PetscCall(VecGetOwnershipRange(parallel, &localParallelLowerIndex, &localParallelUpperIndex));

    PetscInt localParallelSize = localParallelUpperIndex - localParallelLowerIndex;

    if (stdVector.size() < localParallelLowerIndex + localParallelSize )
    {
        return PETSC_ERR_INT_OVERFLOW;
    }

    for (PetscInt i = 0; i < localParallelSize; i++)
    {
        parallelArray[i] = stdVector[i + localParallelLowerIndex];
    }

    PetscCall(VecRestoreArray(parallel, &parallelArray));

    return errorStatus;
}

PetscErrorCode PetscExtensions::VecLeftShift(Vec input, Vec *output)
{
    PetscErrorCode errorStatus = 0;

    Vec sequentialInput;
    Vec sequentialoutput;

    // PetscCall(VecDuplicate(input, &output));

    PetscCall(VecSequentialFromParallel(input, &sequentialInput));
    PetscCall(VecSequentialFromParallel(*output, &sequentialoutput));

    PetscInt inputSize;
    PetscScalar *outputArray;
    const PetscScalar *inputArray;

    PetscCall(VecGetSize(sequentialInput, &inputSize));
    PetscCall(VecGetArray(sequentialoutput, &outputArray));
    PetscCall(VecGetArrayRead(sequentialInput, &inputArray));

    for (PetscInt i = 0; i < inputSize - 1; i++)
    {
        outputArray[i] = inputArray[i+1];
    }
    outputArray[inputSize-1] = 0;

    PetscCall(VecRestoreArray(sequentialoutput, &outputArray));
    PetscCall(VecRestoreArrayRead(sequentialInput, &inputArray));

    PetscCall(VecParallelFromSequential(sequentialInput, &input));
    PetscCall(VecParallelFromSequential(sequentialoutput, output));

    PetscCall(VecDestroy(&sequentialInput));
    PetscCall(VecDestroy(&sequentialoutput));

    return errorStatus;
}

PetscErrorCode PetscExtensions::PrintVecMinMaxMean(Vec vector, const char * name)
{
    PetscReal minPos;
    PetscReal maxPos;
    PetscScalar meanPos;

    PetscCall(VecMax(vector, NULL, &maxPos));
    PetscCall(VecMin(vector, NULL, &minPos));
    PetscCall(VecMean(vector, &meanPos));

    PetscPrintf(PETSC_COMM_WORLD, "\n-------------------------------------------\n%s", name);
    PetscPrintf(PETSC_COMM_WORLD, "\nmax: %.20e, min: %.20e, mean: %.20e\n", maxPos, minPos, meanPos);

    return 0;
}

PetscErrorCode PetscExtensions::PrintVec(Vec vector, const char * name)
{
    PetscErrorCode errorStatus = 0;
    Vec temp;
    PetscCall(PetscExtensions::VecSequentialFromParallel(vector, &temp));

    const PetscScalar *tempArray;
    PetscInt tempSize;
    PetscCall(VecGetSize(temp, &tempSize));
    PetscCall(VecGetArrayRead(temp, &tempArray));

    PetscPrintf(PETSC_COMM_WORLD, "\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
    PetscPrintf(PETSC_COMM_WORLD, "Printing %s:\n", name);
    for (PetscInt i = 0; i < tempSize; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%d: %f\n", i,tempArray[i]);
    }
    PetscPrintf(PETSC_COMM_WORLD, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");

    PetscCall(VecRestoreArrayRead(temp, &tempArray));
    PetscCall(VecDestroy(&temp));

    return errorStatus;
}