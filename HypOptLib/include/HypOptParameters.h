
#pragma once

#include <petsc.h>

typedef struct
{
    Vec position;
    Vec velocity;
    Vec evenNoseHooverPosition;
    Vec evenNoseHooverVelocity;
    Vec oddNoseHooverPosition;
    Vec oddNoseHooverVelocity;
} HypOptParameters;
