/*
#  Copyright(C) 2015-2018  all rights reserved

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

#include<vector>
#include <string>
using namespace std;

class BootLineaRegress
{
public:
	void mf_linrearRegress(vector<double> x, vector<double> y, double& a, double &b, double &DevSQ, double &sd, double &RSS, double &Maxd, double &Mind, double& Meand);
	vector<double> mf_loadVectorData(string path);
	void mf_Normalization(vector<double>v, vector<double>&Normv, double &sd);
	double mf_mean(vector<double> v);
	void mf_BootstrapLinearRegression(vector<double>x, vector<double> y, int nboot, double& a, double &b, double& SdOfa, double &SdOfb, double &DevSQ, double &sd, double &RSS, double &Maxd, double &Mind, double& Meand);
	void mf_Predict(vector<double> intensity, double a, double b, vector<double>& mols);
	double mf_Predict(double intensity, double a, double b);

};
