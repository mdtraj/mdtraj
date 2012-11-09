/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: xdrfile_c_test.c,v 1.7 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */

/* Get HAVE_RPC_XDR_H, F77_FUNC from config.h if available */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifdef HAVE_UNISTD
#include "unistd.h"
#endif

/* get fixed-width types if we are using ANSI C99 */
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#elif (defined HAVE_INTTYPES_H)
#  include <inttypes.h>
#endif

#ifdef HAVE_RPC_XDR_H
#  include <rpc/rpc.h>
#  include <rpc/xdr.h>
#endif

#include <time.h>
#include <float.h>
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

static void _die_r(char *msg, int result,int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr,"result = %d\n",result);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die_r(msg,res) _die_r(msg,res,__LINE__,__FILE__)

static void test_xtc()
{
	char *testfn = "test.xtc";
	XDRFILE *xd;
	int result,i,j,k,nframes=13;
	int natoms2,natoms1=173;
	int step2,step1=1993;
	float time2,time1=1097.23;
	matrix box2,box1;
	rvec *x2,*x1;
	float prec2,prec1=1000;
	float toler=1e-3;
	
	printf("Testing xtc functionality:");
	for(i=0; (i<DIM); i++)
		for(j=0; (j<DIM); j++)
			box1[i][j] = (i+1)*3.7 + (j+1);
	x1 = calloc(natoms1,sizeof(*x1));
	if (NULL == x1)
		die("Allocating memory for x1 in test_xtc");
	
	for(i=0; (i<natoms1); i++)
		for(j=0; (j<DIM); j++)
			x1[i][j] = (i+1)*3.7 + (j+1);
	xd = xdrfile_open(testfn,"w");
	if (NULL == xd)
		die("Opening xdrfile for writing");
	for(k=0; (k<nframes); k++)
		{
			result = write_xtc(xd,natoms1,step1+k,time1+k,box1,x1,prec1);
			if (0 != result)
				die_r("Writing xtc file",result);
		}
	xdrfile_close(xd);
	
	result = read_xtc_natoms(testfn,&natoms2);
	if (exdrOK != result)
		die_r("read_xtc_natoms",result);
	if (natoms2 != natoms1)
		die("Number of atoms incorrect when reading trr");
	x2 = calloc(natoms2,sizeof(x2[0]));
	if (NULL == x2)
		die("Allocating memory for x2");
		
		
	xd = xdrfile_open(testfn,"r");
	if (NULL == xd)
		die("Opening xdrfile for reading");
	
	k = 0;
	do
		{
			result = read_xtc(xd,natoms2,&step2,&time2,box2,x2,&prec2);
			if (exdrENDOFFILE != result)
				{
					if (exdrOK != result)
						die_r("read_xtc",result);
					if (natoms2 != natoms1)
						die("natoms2 != natoms1");
					if (step2-step1 != k)
						die("incorrect step");
					if (fabs(time2-time1-k) > toler)
						die("incorrect time");
					if (fabs(prec2-prec1) > toler)
						die("incorrect precision");
					for(i=0; (i<DIM); i++)
						for(j=0; (j<DIM); j++)
							if (fabs(box2[i][j] - box1[i][j]) > toler)
								die("box incorrect");
					for(i=0; (i<natoms1); i++)
						for(j=0; (j<DIM); j++)
							if (fabs(x2[i][j] - x1[i][j]) > toler)
								die("x incorrect");
				}
			k++;
		} while (result == exdrOK);
		
	xdrfile_close(xd);
#ifdef HAVE_UNISTD
	unlink(testfn);
#endif
	printf(" PASSED\n");
}

static void test_trr()
{
	char *testfn = "test.trr";
	XDRFILE *xd;
	int result,i,j,k,nframes=13;
	int natoms2,natoms1=173;
	int step2,step1=1993;
	float time2,time1=1097.23;
	matrix box2,box1;
	rvec *x2,*x1;
	float lambda2,lambda1=0.4;
	float toler=1e-3;
	
	printf("Testing trr functionality:");
	for(i=0; (i<DIM); i++)
		for(j=0; (j<DIM); j++)
			box1[i][j] = (i+1)*3.7 + (j+1);
	x1 = calloc(natoms1,sizeof(*x1));
	if (NULL == x1)
		die("Allocating memory for x1 in test_xtc");
	
	for(i=0; (i<natoms1); i++)
		for(j=0; (j<DIM); j++)
			x1[i][j] = (i+1)*3.7 + (j+1);
	xd = xdrfile_open(testfn,"w");
	if (NULL == xd)
		die("Opening trr file for writing");
	for(k=0; (k<nframes); k++)
		{
			result = write_trr(xd,natoms1,step1+k,time1+k,lambda1,box1,x1,NULL,NULL);
			if (0 != result)
				die_r("Writing trr file",result);
		}
	xdrfile_close(xd);
	
	result = read_trr_natoms(testfn,&natoms2);
	if (exdrOK != result)
		die_r("read_trr_natoms",result);
	if (natoms2 != natoms1)
		die("Number of atoms incorrect when reading trr");
	x2 = calloc(natoms2,sizeof(x2[0]));
	if (NULL == x2)
		die("Allocating memory for x2");
		
		
	xd = xdrfile_open(testfn,"r");
	if (NULL == xd)
		die("Opening trr file for reading");
	
	for(k=0; (k<nframes); k++) 
		{
			result = read_trr(xd,natoms2,&step2,&time2,&lambda2,box2,x2,NULL,NULL);
			if (exdrOK != result)
				die_r("read_xtc",result);
			if (natoms2 != natoms1)
				die("natoms2 != natoms1");
			if (step2-step1 != k)
				die("incorrect step");
			if (fabs(time2-time1-k) > toler)
				die("incorrect time");
			if (fabs(lambda2-lambda2) > toler)
				die("incorrect lambda");
			for(i=0; (i<DIM); i++)
				for(j=0; (j<DIM); j++)
					if (fabs(box2[i][j] - box1[i][j]) > toler)
						die("box incorrect");
			for(i=0; (i<natoms1); i++)
				for(j=0; (j<DIM); j++)
					if (fabs(x2[i][j] - x1[i][j]) > toler)
						die("x incorrect");
		} 
	
	xdrfile_close(xd);
#ifdef HAVE_UNISTD
	unlink(testfn);
#endif
	printf(" PASSED\n");
}

static void test_basic()
{
	float test_ii;      // 7 significant digits
	double test_jj;     // 16 significant digits
#define EPSILON_1 1e-7
#define EPSILON_2 1e-4
	
	printf("Testing basic xdrfile library:");
	for (test_ii = 1.0e1; test_ii < 1.0e2; (test_ii = test_ii + pow(M_PI, 0.00011)))
		{
		
#define BUFLEN 37
	XDRFILE *xfp;
	int i, j, k, len, ncoord = BUFLEN/3;
	char ptr[BUFLEN], *buf = "abcdefghijklmnopqrstuvwxyz";
	char *testfn = "test.xdr";
	unsigned char uptr[BUFLEN];
	short sptr[BUFLEN], sptr2[BUFLEN];
	unsigned short usptr[BUFLEN], usptr2[BUFLEN];
	int iptr[BUFLEN], iptr2[BUFLEN];
	unsigned int uiptr[BUFLEN], uiptr2[BUFLEN];
	float fptr[BUFLEN], fptr2[BUFLEN];
	double dptr[BUFLEN], dptr2[BUFLEN];
	char optr[BUFLEN], optr2[BUFLEN];
#define NPREC 1

	float  fprec[] = {234.45};
	double dprec[] = {234.45};

	/* Can not write a string that's on the stack since all data is
	   treated as variables.
	 */
	len = strlen(buf) + 1;
	if (len >= BUFLEN)
		die("Increase BUFLEN");
	strcpy(ptr, buf);
	strcpy((char *) uptr, buf);
	/* Initiate float arrays */
	for (i = 0; (i < BUFLEN); i++) 
	{
		fptr[i] = cos(i * 13.0 / M_PI);
		dptr[i] = sin(i * 13.0 / M_PI);
	}
	/* Initiate opaque array */
	memcpy(optr, dptr, BUFLEN);

	/*************************************/
	/*           WRITING BIT             */
	/*************************************/

	if ((xfp = xdrfile_open("test.xdr", "w")) == NULL)
		die("Can not open file for writing");

	if (xdrfile_write_char(ptr, len, xfp) != len)
		die("Writing char string");
	if (xdrfile_write_uchar(uptr, len, xfp) != len)
		die("Writing uchar string");
	if (xdrfile_write_short(sptr, BUFLEN,xfp) != BUFLEN)
		die("Writing short array");
	if (xdrfile_write_ushort(usptr, BUFLEN,xfp) != BUFLEN)
		die("Writing ushort array");
	if (xdrfile_write_int(iptr, BUFLEN,xfp) != BUFLEN)
		die("Writing int array");
	if (xdrfile_write_uint(uiptr, BUFLEN,xfp) != BUFLEN)
		die("Writing uint array");
	if (xdrfile_write_float(fptr, BUFLEN,xfp) != BUFLEN)
		die("Writing float array");
	if (xdrfile_write_double(dptr, BUFLEN,xfp) != BUFLEN)
		die("Writing double array");
	if (xdrfile_write_string(buf, xfp) != len)
		die("Writing string");
	if (xdrfile_write_opaque(optr, BUFLEN,xfp) != BUFLEN)
		die("Writing opaque");
	for (k = 0; (k < NPREC); k++) {
		if (xdrfile_compress_coord_float(fptr, ncoord, fprec[k], xfp)
				!= ncoord)
			die("Writing compress_coord_float");
		if (xdrfile_compress_coord_double(dptr, ncoord, dprec[k], xfp)
				!= ncoord)
			die("Writing compress_coord_double");
	}
	if (xdrfile_close(xfp) != 0)
		die("Can not close xdr file");

	/*************************************/
	/*          READING BIT              */
	/*************************************/
	if ((xfp = xdrfile_open(testfn, "r")) == NULL)
		die("Can not open file for reading");

	if ((xdrfile_read_char(ptr, len, xfp)) != len)
		die("Not the right number of chars read from string");
	if (strcmp(ptr, buf) != 0)
		printf("did not read the expected chars");
	if (xdrfile_read_uchar(uptr, len, xfp) != len)
		die("Not the right number of uchars read from string");
	if (strcmp((char *) uptr, buf) != 0)
		printf("did not read the expected uchars");
	if (xdrfile_read_short(sptr2, BUFLEN,xfp) != BUFLEN)
		die("Reading short array");
	for (i = 0; (i < BUFLEN); i++)
		if (sptr2[i] != sptr[i]) {
			fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, sptr[i],
					sptr2[i]);
			die("Comparing short array");
		}
	if (xdrfile_read_ushort(usptr2, BUFLEN,xfp) != BUFLEN)
		die("Reading ushort array");
	for (i = 0; (i < BUFLEN); i++)
		if (usptr2[i] != usptr[i]) {
			fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, usptr[i],
					usptr2[i]);
			die("Comparing ushort array");
		}
	if (xdrfile_read_int(iptr2, BUFLEN,xfp) != BUFLEN)
		die("Reading int array");
	for (i = 0; (i < BUFLEN); i++)
		if (iptr2[i] != iptr[i]) {
			fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, iptr[i],
					iptr2[i]);
			die("Comparing int array");
		}
	if (xdrfile_read_uint(uiptr2, BUFLEN,xfp) != BUFLEN)
		die("Reading uint array");
	for (i = 0; (i < BUFLEN); i++)
		if (uiptr2[i] != uiptr[i]) {
			fprintf(stderr, "i: %5d, wrote: %10d, read: %10d\n", i, uiptr[i],
					uiptr2[i]);
			die("Comparing uint array");
		}
	if (xdrfile_read_float(fptr2, BUFLEN,xfp) != BUFLEN)
		die("Reading float array");
	for (i = 0; (i < BUFLEN); i++)
		if (fptr2[i] != fptr[i]) {
			fprintf(stderr, "i: %5d, wrote: %12g, read: %12g\n", i, fptr[i],
					fptr2[i]);
			die("Comparing float array");
		}
	if (xdrfile_read_double(dptr2, BUFLEN,xfp) != BUFLEN)
		die("Reading double array");
	for (i = 0; (i < BUFLEN); i++)
		if (dptr2[i] != dptr[i]) {
			fprintf(stderr, "i: %5d, wrote: %12g, read: %12g\n", i, dptr[i],
					dptr2[i]);
			die("Comparing double array");
		}
	if (xdrfile_read_string(ptr, BUFLEN,xfp) != len)
		die("Reading string");
	if (strcmp(ptr, buf) != 0)
		die("Comparing strings");
	if (xdrfile_read_opaque(optr2, BUFLEN,xfp) != BUFLEN)
		die("Reading opaque array");
	for (i = 0; (i < BUFLEN); i++)
		if (optr2[i] != optr[i]) {
			fprintf(stderr, "i: %5d, wrote: %2d, read: %2d\n", i, optr[i],
					optr2[i]);
			die("Comparing opaque array");
		}
	for (k = 0; (k < NPREC); k++) {
		float ff, fx;
		double dd, dx;
		int nc = ncoord;
		if (xdrfile_decompress_coord_float(fptr2, &nc, &ff, xfp) != ncoord)
			die("Reading compress_coord_float");
		if (fabs(ff - fprec[k]) > EPSILON_1) 
			{
				printf("Found precision %f, expected %f\n", ff, fprec[k]);
				die("Float precision");
			}
		if (ff <= 0)
			ff = 1000;


		for (i = 0; (i < ncoord); i++)
			for (j = 0; (j < 3); j++) {
				fx = rint(fptr[3 * i + j] * ff) / ff;
				if (fabs(fx - fptr2[3 * i + j]) > EPSILON_1) {
					printf(	"prec: %10g, i: %3d, j: %d, fx: %10g, fptr2: %12g, fptr: %12g\n",
							ff, i, j, fx, fptr2[3 * i + j], fptr[3 * i + j]);
					die("Reading decompressed float coordinates");
				}
			}
		if (xdrfile_decompress_coord_double(dptr2, &nc, &dd, xfp) != ncoord)
			die("Reading compress_coord_double");

		if (fabs(dd - dprec[k]) > EPSILON_2)
			die("Double precision");

		for (i = 0; (i < ncoord); i++)
			for (j = 0; (j < 3); j++) {
				dx = rint(dptr[3 * i + j] * dd) / dd;
				if (fabs(dx - dptr2[3 * i + j]) > EPSILON_2) {
					printf(	"prec: %10g, i: %3d, j: %d, dx: %10g, dptr2: %12g, dptr: %12g\n",
							dd, i, j, dx, dptr2[3 * i + j], dptr[3 * i + j]);
					die("Reading decompressed double coordinates");
				}
			}
	}

	if (xdrfile_close(xfp) != 0)
		die("Can not close xdr file");
		} 
#ifdef HAVE_UNISTD
	unlink(testfn);
#endif
	printf(" PASSED\n");
}

int main(int argc, char *argv[]) 
{
	/* Test basic stuff */
	test_basic();
	/* Now test writing a complete xtc file */
	test_xtc();

	test_trr();
		
	return 0;
}
