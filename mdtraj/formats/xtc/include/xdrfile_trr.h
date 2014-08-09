/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: xdrfile_trr.h,v 1.2 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */

#ifndef _xdrfile_trr_h
#define _xdrfile_trr_h

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "xdrfile.h"
  
  /* All functions return exdrOK if succesfull. 
   * (error codes defined in xdrfile.h).
   */  
   
  /* This function returns the number of atoms in the xtc file in *natoms */
  extern int read_trr_natoms(char *fn,int *natoms);
  extern int read_trr_nframes(char* fn, unsigned long *nframes);
  
  /* Read one frame of an open xtc file. If either of x,v,f,box are
     NULL the arrays will be read from the file but not used.  */
  extern int read_trr(XDRFILE *xd,int natoms,int *step,float *t,float *lambda,
		      matrix box,rvec *x,rvec *v,rvec *f);

  /* Write a frame to xtc file */
  extern int write_trr(XDRFILE *xd,int natoms,int step,float t,float lambda,
		       matrix box,rvec *x,rvec *v,rvec *f);

  
#ifdef CPLUSPLUS
}
#endif

#endif
