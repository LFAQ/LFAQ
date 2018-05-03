/*
#	This part is done based on the source code in http://www.rob-mcculloch.org/
#
#  This program is a free software; you can redistribute it and / or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/Licenses/
*/
#ifndef GUARD_info_h
#define GUARD_info_h


//data
class dinfo {
public:
	dinfo() { p = 0; n = 0; x = 0; y = 0; }
	size_t p;  //number of vars
	size_t n;  //number of observations
	double *x; // jth var of ith obs is *(x + p*i+j)
	double *y; // ith y is *(y+i) or y[i]
};

//prior and mcmc
class pinfo
{
public:
	//pinfo() { pbd = 1.0; pb = .5; alpha = .95; beta = .5; tau = 1.0; sigma = 1.0; }
	pinfo() {  pb = .5; alpha = .95; beta = .5; tau = 1.0; sigma = 1.0; }

	//mcmc info
	double pb;  //prob of birth
	double alpha;
	double beta;
	double tau;
	double sigma;
};

//sufficient statistics for 1 node
class sinfo
{
public:
	sinfo() { n = 0; sy = 0.0; sy2 = 0.0; }
void clear() { n = 0; sy = 0.0; sy2 = 0.0; }
	size_t n;
	double sy;
	double sy2;
};

#endif
