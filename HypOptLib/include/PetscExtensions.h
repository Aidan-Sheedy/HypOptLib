/***************************************************************************//**
 * @file PetscExtensions.h
 *
 * Contains a class with useful functions not included in the standard Petsc
 * library.
 *
 * @author Aidan Sheedy
 *
 * Copyright (C) 2024 Aidan Sheedy
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 ******************************************************************************/

#pragma once

#include <petsc.h>
#include <vector>

class PetscExtensions
{
    public:
        /** Delete all constructors, all access should be static */
        PetscExtensions() = delete;

        /** Delete all constructors, all access should be static */
        PetscExtensions(const PetscExtensions&) = delete;

        /**
         * Creates a parallel vector from a given sequential vector.
         *
         * The parallel vector must already be created using either VecCreate followed by VecSetFromOptions,
         * or from a VecDuplicate call.
         *
         * @param sequential the sequential vector to prallelize.
         * @param parallel [out] the parallel vector to populate.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode VecParallelFromSequential(Vec sequential, Vec *parallel);

        /**
         * Creates a sequential vector from a given parallel vector.
         *
         * The sequential vector will be created in the function, so a fresh Vec object can be passed.
         *
         * @param parallel the parallel vector to sequentialize.
         * @param sequential [out] the sequential vector to populate.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode VecSequentialFromParallel(Vec parallel, Vec *sequential);

        /**
         * Fetches a value from a parallel vector that is not stored in the local vector.
         *
         * This may not be the most efficient way of doing so, but it is fast enough not to bottleneck
         * the code with large finite element analysis meshes to solve. It works as follows:
         *
         * 1. Get a local sequential copy of the parallel vector of length 1
         * 2. Perform scatter operation for only the given index to retrieve
         * 3. Clean up and return the result
         *
         * @param parallelVector the parallel vector to access.
         * @param index the index to retrieve from
         * @param result [out] resulting value from the desired index.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode VecGetOffProcessIndex(Vec parallelVector, PetscInt index, PetscScalar *result);

        /**
         * Creates a parallel vector from a given C++ std vector.
         *
         * The parallel vector must already be created using either VecCreate followed by VecSetFromOptions,
         * or from a VecDuplicate call. This function simply populates the parallel vector with the contents of
         * the std vector. This function is parallelized.
         * 
         * @note This version of the overloaded function uses a templated vector type which is cast to PetscScalar.
         * It should only be used if the vector type can safely be cast to PetscScalar.
         *
         * @param stdVector the std vector to prallelize.
         * @param parallel [out] the parallel vector to populate.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        template<typename T> 
            static PetscErrorCode VecParallelFromStdVector(std::vector<T> stdVector, Vec parallel)
            {
                std::vector<PetscScalar> copy;
                for (auto value : stdVector)
                {
                    copy.push_back((PetscScalar)value);
                }
                return VecParallelFromStdVector(copy, parallel);
            }

        /**
         * Creates a parallel vector from a given C++ std vector.
         *
         * The parallel vector must already be created using either VecCreate followed by VecSetFromOptions,
         * or from a VecDuplicate call. This function simply populates the parallel vector with the contents of
         * the std vector. This function is parallelized.
         *
         * @param stdVector the std vector to prallelize.
         * @param parallel [out] the parallel vector to populate.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode VecParallelFromStdVector(std::vector<PetscScalar> stdVector, Vec parallel);


        /**
         * Shifts all elements in the input vector down one index.
         *
         * Note that this is an explicit shift, not a rotate. The last element in the vector is set to 0.
         * This is a parallelized implementation, and is not optimized for sequential vectors.
         *
         * @param input the vector to shift.
         * @param output [out] the shifted vector.
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode VecLeftShift(Vec input, Vec *output);

        /**
         * Convenience function to print the minimum, maximum, and mean values of a vector.
         *
         * @param vector the vector to print values of.
         * @param name name of the vector to print to print (for debug purposes).
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode PrintVecMinMaxMean(Vec vector, const char * name);

        /**
         * Convenience function to print all elements in a vector.
         *
         * @param vector the vector to print values of.
         * @param name name of the vector to print to print (for debug purposes).
         *
         * @returns 0 on success, PetscError otherwise.
         */
        static PetscErrorCode PrintVec(Vec vector, const char * name);
};
