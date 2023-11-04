
#pragma once

#include <exception>

class HypOptException: public std::exception
{
    public:
        HypOptException(const char* msg) : msg(msg) {}

        virtual const char* what() const throw()
        {
            return msg;
        }

        private:
            const char* msg;
};
