/***************************************************************************//**
 * @file HypOptException.h
 *
 * This file implements the custom Hyperopimization exception class.
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
