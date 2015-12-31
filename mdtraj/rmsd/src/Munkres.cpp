// This code was downloaded from https://github.com/jfrelinger/cython-munkres-wrapper
// accessed March 13, 2014.
//
// Copyright (c) 2012, Jacob Frelinger
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//////////////////////////////////////////////////////////////////

/*
 * Munkres.cpp
 *
 *  Created on: Sep 29, 2010
 *      Author: jolly
 */

#include "Munkres.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
#include <algorithm>

//using namespace std;
//using std::vector;

path_item::path_item(int i, int j, path_type p_or_s) {
	row = i;
	col = j;
	type = p_or_s;
}

path_item::~path_item() {

}

Munkres::Munkres() {
	// TODO Auto-generated constructor stub

}

Munkres::~Munkres() {
	// TODO Auto-generated destructor stub
}

void Munkres::solve(double* icost, int* answer, int m, int n) {
	rows = m;
	cols = n;
	cost = new double*[rows];
	starred = new bool*[rows];

	primed = new bool*[rows];
	covered_rows = new bool[rows];
	covered_cols = new bool[cols];

	for (int i = 0; i < rows; i++) {
		covered_rows[i] = false;
	}
	for (int i = 0; i < cols; i++) {
		covered_cols[i] = false;
	}

	for (int i = 0; i < rows; i++) {
		cost[i] = new double[cols];
		starred[i] = new bool[cols];
		primed[i] = new bool[cols];

		for (int j = 0; j < cols; j++) {
			cost[i][j] = icost[(i * cols) + j];
			starred[i][j] = 0;
			primed[i][j] = 0;
		}
	}

	smallest = std::min(rows, cols);
	largest = std::max(rows, cols);

	if (rows > cols) {
		step0();
	} else {
		k = min_uncovered();

		step1();
	}

	int index = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			answer[index] = starred[i][j];
			index++;
		}
	}

	for (int i = 0; i < rows; i++) {
		delete cost[i];
		delete primed[i];
		delete starred[i];
	}
	delete cost;
	delete primed;
	delete starred;
	delete covered_rows;
	delete covered_cols;

}

void Munkres::step0() {
	int minimum;
	for (int j = 0; j < cols; j++) {
		minimum = cost[0][j];
		for (int i = 0; i < rows; i++) {
			if (minimum > cost[i][j]) {
				minimum = cost[i][j];
			}
		}
		for (int i = 0; i < rows; i++) {
			cost[i][j] = cost[i][j] - minimum;
		}
	}
	k = min_uncovered();
	step2();
}

void Munkres::step1() {
	/*
	 * subtract the smallest element of each row from that row.
	 * goto step 2
	 */
	for (int i = 0; i < rows; i++) {
		double * a = cost[i];
		double m = std::numeric_limits<double>::infinity();
		for (int cii = 0; cii < cols; cii++) {
			if (m > a[cii]) {
				m = a[cii];
			}
		}

		for (int cii = 0; cii < cols; cii++) {
			a[cii] = a[cii] - m;
		}

	}
	step2();
}
void Munkres::step2() {
	/*
	 * find a zero, if now starred zeros in row or column star z
	 * repeat for each element.
	 *
	 * goto step 3
	 */
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			if (cost[i][j] == 0) {
				if (!is_starred_in_row_col(i, j)) {
					starred[i][j] = 1;
				}
			}
		}
	}
	step3();
}
void Munkres::step3() {
	/* cover each coulumn containing a starred zero
	 * if size covered columns we're done.
	 * else goto step 4
	 */
	int cov_count = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (starred[i][j] == 1) {
				cover_col(j);
				cov_count += 1;
			}
		}
	}
	if (cov_count != smallest) {
		step4();
	}

}
void Munkres::step4() {
	/* find an uncovered zero and prime it
	 * if now starred zeros in row goto step 5
	 * else cover row and uncover colum of starred zero
	 *
	 * once no uncovered exist goto step 6.
	 */
	bool done = false;
	int i = 0;
	int j = 0;
	int a;
	while (!done) {
		if (find_zero(&i, &j)) {
			if (!is_covered(i, j)) {
				prime(i, j);
				a = starred_in_row(i);
				if (a == -1) // if no starred zeros
				{
					done = true;
					step5(i, j);
				} else {
					uncover_col(a);
					cover_row(i);
				}

			}
		} else {
			done = true;
			step6( min_uncovered());
		}
	}
}
void Munkres::step5(int i, int j) {
	/* take a primed zero, and construct a list of...
	 * 1. a starred zero in it's column (if it exists)
	 * 2. if there's a starred zero there will be a primed zero in its row.
	 *
	 * then
	 * unstar the starred,
	 * star all the primes
	 * erase all primes
	 * uncover everything
	 * return to step 3.
	 */
	std::vector<path_item> path;
	path.push_back(path_item(i, j, PRIMED));
	bool done = false;
	int row = 0;
	int col = j;
	while (!done) {
		row = find_starred_zero_in_col(col);
		if (row != -1) {
			path.push_back(path_item(row, col, STARRED));
			col = find_primed_zero_in_row(row);
			path.push_back(path_item(row, col, PRIMED));
		} else {
			done = true;
		}
	}

	for (unsigned int i = 0; i < path.size(); i++) {
		path_item item = path[i];
		if (item.type == PRIMED) // primed so we star
		{
			starred[item.row][item.col] = 1;
		} else { // we're starred so we unstar
			starred[item.row][item.col] = 0;
		}
	}
	// remove all primes
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			primed[i][j] = 0;
		}
	}
	for (int i = 0; i < rows; i++) {
		covered_rows[i] = 0;
	}
	// uncover all covered lines
	for (int i = 0; i < rows; i++) {
		covered_rows[i] = 0;
	}
	for (int i = 0; i < cols; i++) {
		covered_cols[i] = 0;
	}
	step3();
}
void Munkres::step6(double val) {
	/* take a value and add it to ever covered row
	 * then subtract it from every uncovered column.
	 * return to step 4
	 */

	for (int i = 0; i < rows; i++) {
		if (is_covered_row(i)) {
			for (int j = 0; j < cols; j++) {
				cost[i][j] += val;
			}
		}
	}
	for (int i = 0; i < cols; i++) {
		if (!is_covered_col(i)) { // uncovered column
			for (int j = 0; j < rows; j++) {
				cost[j][i] -= val;
			}
		}
	}
	step4();
}

bool Munkres::is_starred_in_row_col(int row, int col) {
	bool * a = starred[row];

	for (int r = 0; r < cols; r++) {
		if (a[r] != 0) {
			return true;
		}
	}

	for (int r = 0; r < row; r++) {
		if (starred[r][col] != 0) {
			return true;
		}
	}

	// fall through to false

	return false;
}

int Munkres::starred_in_row(int row) {
	// find a starred value in a row
	bool * a = starred[row];

	for (int i = 0; i < cols; i++) {
		if (a[i] == 1) {
			return i;
		}
	}
	return -1;
}

void Munkres::cover_col(int col) {
	//cover a column
	covered_cols[col] = 1;
}

void Munkres::uncover_col(int col) {
	//uncover a column
	covered_cols[col] = 0;
}

void Munkres::cover_row(int row) {
	// cover a row
	covered_rows[row] = 1;
}

void Munkres::uncover_row(int row) {
	// uncover a row
	covered_rows[row] = 0;
}

bool Munkres::is_covered(int row, int col) {
	// check if a position in a covered row or column
	if ((covered_rows[row] == 1) || (covered_cols[col] == 1)) {
		return true;
	} else {
		return false;
	}
}

bool Munkres::is_covered_col(int col) {
	// check if a column is covered
	if (covered_cols[col] == 1) {
		return true;
	} else {
		return false;
	}
}

bool Munkres::is_covered_row(int row) {
	// check if a row is covered
	if (covered_rows[row] == 1) {
		return true;
	} else {
		return false;
	}
}

void Munkres::prime(int row, int col) {
	// prime a postion
	primed[row][col] = 1;
}

bool Munkres::find_zero(int* row, int* col) {
	// find a zero thats uncovered
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (cost[i][j] == 0) {
				if (!is_covered(i, j)) {
					*row = i;
					*col = j;
					return true;
				}
			}
		}
	}
	return false;
}

double Munkres::min_uncovered() {
	// find the minumum uncovered value in the cost matrix.

	double min = std::numeric_limits<double>::infinity();

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (!is_covered(i, j)) {
				if (cost[i][j] < min) {
					min = cost[i][j];
				}
			}
		}
	}
	return min;
}

int Munkres::find_starred_zero_in_col(int col) {
	// given a column find a starred zero in it, otherwise return -1
	for (int i = 0; i < rows; i++) {
		if (starred[i][col] == true)
			return i;
	}
	return -1;
}

int Munkres::find_primed_zero_in_row(int row) {
	// given a row, find if it has a primed zero in it, otherwise return -1
	for (int i = 0; i < cols; i++) {
		if (primed[row][i] == true)
			return i;
	}
	return -1;
}

