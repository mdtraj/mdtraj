/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2012-2013 Stanford University and the Authors               */
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
Header file for a C struct (C-style object oriented programming)
that computes the mean, second, third and fourth moments of a bunch of numbers
in a single pass.

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

#ifndef _MDTRAJ_MOMENTS_H_
#define _MDTRAJ_MOMENTS_H_

typedef struct {
    int _n;
    double _u;
    double _M2;
    double _M3;
    /* double _M4; */
} moments_t;

/* Push a number onto the object */
void moments_push(moments_t *cls, double x);

/* Clear the object, and forget all data previously seen */
void moments_clear(moments_t *cls);

/* Mean of the numbers which have been pushed onto the stack */
double moments_mean(moments_t *cls);

/* Second centralmoment of the numbers which have been pushed onto the stack */
double moments_second(moments_t *cls);

/* Third central moment of the numbers which have been pushed onto the stack */
double moments_third(moments_t *cls);

/*
// Fourth central moment of the numbers which have been pushed onto the stack
double moments_fourth(moments_t *cls);
*/

#endif
