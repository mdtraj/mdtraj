/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: trr2xtc.c,v 1.3 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"


/* This program tests reading and writing to XDR files */

static void _die(char *msg, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die(msg) _die(msg,__LINE__,__FILE__)

static void _die_r(char *msg, int result, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr,"result = %d\n", result);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die_r(msg,res) _die_r(msg,res,__LINE__,__FILE__)


void ReadWrite(char *rfile, char *wfile, int in_xtcBool, int out_xtcBool, int in_trrBool, int out_trrBool)
{
  XDRFILE *xd_read, *xd_write;
  int result_xtc, result_trr;
  int natoms_xtc, natoms_trr;
  int step_xtc, step_trr;
  float time_xtc, time_trr;
  matrix box_xtc, box_trr;
  rvec *x_xtc, *x_trr, *v_trr, *f_trr;
  float prec_xtc = 1000.0;
  float lambda_trr = 0.0;
  
      
  xd_read = xdrfile_open(rfile, "r");
  if (NULL == xd_read)
    die("Opening xdrfile for reading");
  
  /* Test whether output file exists */
  if ((xd_write = xdrfile_open(wfile,"r")) != NULL) {
    xdrfile_close(xd_write);
    die("Output file exists.");
  }
  
  /* Output file does not exist. Now we can open it for writing */
  xd_write = xdrfile_open(wfile, "w");
  if (NULL == xd_write)
    die("Opening xdrfile for writing");
  
  
  /* .xtc -> .xtc */
  if(in_xtcBool && out_xtcBool)
    {
      result_xtc = read_xtc_natoms(rfile, &natoms_xtc);
      if (exdrOK != result_xtc)
	die_r("read_xtc_natoms",result_xtc);
      
      x_xtc = calloc(natoms_xtc, sizeof(x_xtc[0]));

      while(1)
	{
	  result_xtc = read_xtc(xd_read, natoms_xtc, &step_xtc, &time_xtc,
				box_xtc, x_xtc, &prec_xtc);
	    
	  if (result_xtc == 0) // if not reach the end of file, write it to the output.xtc file
	    {
	      if (exdrOK != result_xtc)
		die_r("Reading xtc file", result_xtc);
	      
	      result_xtc = write_xtc(xd_write, natoms_xtc, step_xtc, time_xtc,
				     box_xtc, x_xtc, prec_xtc);
	      
	      if (result_xtc != 0)
		die_r("Writing xtc file", result_xtc);
	    }
	  else
	    break;
	}
    }
  

  /* .xtc -> .trr */
  if(in_xtcBool && out_trrBool)
    {
      result_xtc = read_xtc_natoms(rfile, &natoms_xtc);
      if (exdrOK != result_xtc)
	die_r("read_xtc_natoms",result_xtc);
      
      x_xtc = calloc(natoms_xtc, sizeof(x_xtc[0]));
      
      while(1)
	{
	  result_xtc = read_xtc(xd_read, natoms_xtc, &step_xtc, &time_xtc,
				box_xtc, x_xtc, &prec_xtc);
	  
	  if (result_xtc == 0) // if not reach the end of file, write it to the output.trr file
	    {
	      if (exdrOK != result_xtc)
		die_r("Reading xtc file", result_xtc);
	      
	      result_trr = write_trr(xd_write, natoms_xtc, step_xtc, time_xtc, lambda_trr,
				     box_xtc, x_xtc, NULL, NULL);
	      
	      if (0 != result_trr)
		die_r("Writing trr file",result_trr);
	      
	    }
	  else
	    break;
	}
    }
  
  
  /* .trr -> .trr */
  if(in_trrBool && out_trrBool)
    {
      result_trr = read_trr_natoms(rfile, &natoms_trr);
      
      if (exdrOK != result_trr)
	die_r("read_trr_natoms",result_trr);
      
      x_trr = calloc(natoms_trr, sizeof(x_trr[0]));
      v_trr = calloc(natoms_trr, sizeof(v_trr[0]));
      f_trr = calloc(natoms_trr, sizeof(f_trr[0]));
      
      while (1)
	{
	  result_trr = read_trr(xd_read, natoms_trr, &step_trr, &time_trr, &lambda_trr,
				box_trr, x_trr, v_trr, f_trr);
	      
	  int ii_trr, jj_trr, x_ck=0, v_ck=0, f_ck=0;
	  int x_ck_bool=0, v_ck_bool=0, f_ck_bool=0;
	  
	  for (ii_trr = 0; ii_trr < natoms_trr; ii_trr++)
	    {
	      for(jj_trr = 0; jj_trr < DIM; jj_trr++)
		{
		  if (x_trr[ii_trr][jj_trr] == 0)
		    x_ck++;
		  if (v_trr[ii_trr][jj_trr] == 0)
		    v_ck++;
		  if (f_trr[ii_trr][jj_trr] == 0)
		    f_ck++;
		}
	    }
	  
	  if (x_ck == natoms_trr*DIM)
	    x_ck_bool = 1;
	  if (v_ck == natoms_trr*DIM)
	    v_ck_bool = 1;
	  if (f_ck == natoms_trr*DIM)
	    f_ck_bool = 1;
	      
	  if (result_trr == 0) // if not reach the end of file, write it to the output.trr file
	    {
	      if (exdrOK != result_trr)
		die_r("Reading trr file",result_trr);
	      
	      if(v_ck_bool && f_ck_bool)
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, NULL, NULL);
	      else if(v_ck_bool)
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, NULL, f_trr);
	      else if(f_ck_bool)
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, v_trr, NULL);
	      else
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, v_trr, f_trr);
	      
	      if (0 != result_trr)
		die_r("Writing trr file",result_trr);
	      
	    }
	  else
	    break;
	}
    }
  
  
  /* .trr -> .xtc */
  if(in_trrBool && out_xtcBool)
    {
      result_trr = read_trr_natoms(rfile, &natoms_trr);
      
      if (exdrOK != result_trr)
	die_r("read_trr_natoms",result_trr);
      
      x_trr = calloc(natoms_trr, sizeof(x_trr[0]));
      v_trr = calloc(natoms_trr, sizeof(v_trr[0]));
      f_trr = calloc(natoms_trr, sizeof(f_trr[0]));
      
      while(1)
	{
	  result_trr = read_trr(xd_read, natoms_trr, &step_trr, &time_trr, &lambda_trr,
				box_trr, x_trr, v_trr, f_trr);
	  
	  if (result_trr == 0) // if not reach the end of file, write it to the output.trr file
	    {
	      if (exdrOK != result_trr)
		die_r("Reading trr file", result_trr);
	      
	      result_xtc = write_xtc(xd_write, natoms_trr, step_trr, time_trr,
				     box_trr, x_trr, prec_xtc);
	      
	      if (result_xtc != 0)
		die_r("Writing xtc file", result_xtc);
		}
	  else
	    break;
	}
    }
  
  xdrfile_close(xd_read);
  xdrfile_close(xd_write);
  
}



int main(int argc, char *argv[])
{
  int inFileBool = 0;
  int outFileBool = 0;
  int in_xtcBool = 0;
  int out_xtcBool = 0;
  int in_trrBool = 0;
  int out_trrBool = 0;
  char *rfile=NULL, *wfile=NULL;
  

  if(argc != 5)
    {
      fprintf(stderr,"Usage: %s -i inFile -o outFile\n",argv[0]);
      exit(1);
    }
  else
    {
      int ii = 1;

      while(ii < argc)
	{
	  if(strcmp(argv[ii], "-i") == 0)         // if (argv[ii] == "-i")
	    {
	      ii++;
	      inFileBool = 1;
	 
	      if(strstr(argv[ii], ".xtc") != NULL)
		{
		  in_xtcBool = 1;
		  in_trrBool = 0;
		}
	      if(strstr(argv[ii], ".trr") != NULL)
		{
		  in_trrBool = 1;
		  in_xtcBool = 0;
		}

	      rfile = argv[ii];
	    }

	  if(strcmp(argv[ii], "-o") == 0)
	    {
	      ii++;
	      outFileBool = 1;
	 
	      if(strstr(argv[ii], ".xtc") != NULL)
		{		
		  out_xtcBool = 1;
		  out_trrBool = 0;
		}
	      if(strstr(argv[ii], ".trr") != NULL)
		{
		  out_trrBool = 1;
		  out_xtcBool = 0;
		}

	      wfile = argv[ii];
	    }

	  ii++;
	}
    }

  if(!inFileBool || !outFileBool)
    {
      perror("Usage : ./ReadWrite -i inFile -o outFile");
      exit(1);
    }

  ReadWrite(rfile, wfile, in_xtcBool, out_xtcBool, in_trrBool, out_trrBool);
  

  return 0;
}
