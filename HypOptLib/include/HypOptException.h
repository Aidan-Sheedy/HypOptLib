/***************************************************************************//**
 * @file HypOptException.h
 *
 * This file implements the custom Hyperopimization exception class.
 *
 * @author Aidan Sheedy
 *
 * @todo THIS FILE NEEDS LICENSE INFORMATION
 *
******************************************************************************/

#pragma once

#include <exception>

/**
 * Exception class intended for hyperoptimization related errors.
 *
 * Petsc exception handling is used for the majority of the design loop, but
 * C++ style exceptions are useful especially when propagating errors to the Python
 * wrapper.
 */
class HypOptException: public std::exception
{
    public:
        /** Basic exception constructor. */
        HypOptException(const char* msg) : msg(msg) {}

        /** Exception message function. */
        virtual const char* what() const throw()
        {
            return msg;
        }

        private:
            const char* msg;
};
