/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2014- Stanford University and the Authors                   */
/*                                                                       */
/* Authors: Robert McGibbon                                              */
/* Contributors:                                                         */
/*                                                                       */
/* MDTraj is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as        */
/* published by the Free Software Foundation, either version 2.1         */
/* of the License, or (at your option) any later version.                */
/*                                                                       */
/* This library is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU Lesser General Public License for more details.                   */
/*                                                                       */
/* You should have received a copy of the GNU Lesser General Public      */
/* License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.*/
/*=======================================================================*/

/*
Implementatin for a C struct (C-style object oriented programming) to compute
the mean, second, third and fourth moments of a bunch of numbers in a single
pass.

Compute the mean, variance, 3rd and 4th central moments with a
single pass through the data and O(1) storage

The statandard way to compute the 2nd moment
    \frac{1}{N} \sum_{i=1}^N (X_i - \bar{X})^2

would use a single pass through the data to calculate the mean, and then
a second pass to accumulate the deviations (differences from the mean).

But with some clever work on, you can do it all with one pass

http://www.johndcook.com/standard_deviation.html
https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
http://people.xiph.org/~tterribe/notes/homs.html
*/

#include "moments.h"


/* Clear the object, and forget all data previously seen */
void moments_clear(moments_t *self) {
    self->_n = 0;
    self->_u = 0.0;
    self->_M2 = 0.0;
    self->_M3 = 0.0;
/*
    self->_M4 = 0.0;
*/
}

/* Push a number or a set of numbers onto the RunningMoments */
void moments_push(moments_t *self, double x) {
    int n1;
    double delta, term1, delta_n;  //, delta_n2;

    n1 = self->_n;
    self->_n += 1;
    delta = x - self->_u;
    delta_n = delta / self->_n;
    /* delta_n2 = delta_n*delta_n; */
    term1 = delta * delta_n * n1;
    self->_u += delta_n;
    /*
    self->_M4 += term1 * delta_n2 * (self->_n*self->_n - 3*self->_n + 3) + \
        6 * delta_n2 * self->_M2 - 4 * delta_n * self->_M3;
    */
    self->_M3 += term1 * delta_n * (self->_n - 2) - 3 * delta_n * self->_M2;
    self->_M2 += term1;
}

/* Mean of the numbers which have been pushed onto the stack */
double moments_mean(moments_t *self) {
    return self->_u;
}

/* Second centralmoment of the numbers which have been pushed onto the stack */
double moments_second(moments_t *self) {
    return self->_M2 / self->_n;
}

/* third central moment of the numbers which have been pushed onto the stack */
double moments_third(moments_t *self) {
    return self->_M3 / self->_n;
}

/*
// Fourth central moment of the numbers which have been pushed onto the stack
double moments_fourth(moments_t *cls) {
    return cls->_M4 / cls->_n;
}
*/
