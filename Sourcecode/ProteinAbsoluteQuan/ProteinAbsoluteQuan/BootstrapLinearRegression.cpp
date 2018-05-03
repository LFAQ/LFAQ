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

#include "stdafx.h"
#include"BootstrapLinearRegression.h"

vector<double> BootLineaRegress::mf_loadVectorData(string path)
{

	ifstream infile(path);
	vector<double> vecTemp;
	if (!infile)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " + path + "\n");
		exit(1);
	}

	string lineBuffer;
	double dTemp;
	while (getline(infile, lineBuffer))
	{
		dTemp = atof(lineBuffer.c_str());
		vecTemp.push_back(dTemp);
	}
	return vecTemp;
}
void BootLineaRegress::mf_Normalization(vector<double>v,vector<double>&Normv, double &sd)
{
	Normv.clear();
	double mean = 0.0;
	sd = 0.0;
	int length = v.size();
	for (int i = 0; i < length; i++)
	{
		mean += v.at(i);

	}
	mean = mean / length;

	for (int i = 0; i < length; i++)
	{
		sd += (v.at(i) - mean)*(v.at(i) - mean);
	}
	sd = sqrt(sd / (length-1));

	for (int i = 0; i < length; i++)
	{
		Normv.push_back((v.at(i) - mean) / sd);
	}

}
void BootLineaRegress::mf_linrearRegress(vector<double> x, vector<double> y, double& a, double &b, double &DevSQ, double &sd, double &RSS, double &Maxd, double &Mind, double& Meand)
{
	double meanx = 0.0, meany = 0.0;
	if (x.size() != y.size())
	{
		cout << "Error:\tThe dimension of x and y must be same when doing linear regression.";
		flog.mf_Input("Error:\tThe dimension of x and y must be same when doing linear regression.\n");
		flog.mf_Destroy();
		exit(1);
	}
	int length = x.size();
	for (int i = 0; i < length; i++)
	{
		meanx += x.at(i);
		meany += y.at(i);
	}
	meanx = meanx / length;
	meany = meany / length;

	double tempSum1 = 0.0, tempSum2 = 0.0;
	for (int i = 0; i < length; i++)
	{
		tempSum1 += (x.at(i) - meanx)*(y.at(i) - meany);
		tempSum2 += (x.at(i) - meanx)*(x.at(i) - meanx);
	}
	if (tempSum2 == 0.0)
	{
		a = meany / meanx;
		b = 0;
	}
	else
	{
		a = tempSum1 / tempSum2;
		b = meany - a*meanx;
	}

	double PredictedxMean=0.0;
	for (int i = 0; i < length; i++)
	{
		PredictedxMean += a*x.at(i) + b;
	}
	PredictedxMean = PredictedxMean / length;

	DevSQ = 0.0;    //DevSQ definition: https://support.office.com/en-us/article/DEVSQ-function-8b739616-8376-4df5-8bd0-cfe0a6caf444
	sd = 0.0;
	RSS = 0.0;
	Maxd = 0.0;
	Mind = 0.0;
	Meand = 0.0;
	double tempd;
	for (int i = 0; i < length; i++)
	{
		DevSQ += (PredictedxMean - a*x.at(i) - b)*(PredictedxMean - a*x.at(i) - b);
		RSS += ((a*x.at(i) + b) - meany)*((a*x.at(i) + b) - meany);
		tempd = abs(x.at(i)*a + b - y.at(i));
		if (tempd> Maxd)
		{
			Maxd = tempd;
		}
		if (tempd < Mind)
		{
			Mind = tempd;
		}
		Meand += tempd;

	}
	sd = sqrt(DevSQ / length);
	Meand = Meand / length;
}
double BootLineaRegress::mf_mean(vector<double> v)
{
	double Tempmean = 0.0;
	for (size_t i = 0; i < v.size(); i++)
	{
		Tempmean += v.at(i);
	}
	Tempmean = Tempmean / v.size();
	return Tempmean;
}
void BootLineaRegress::mf_BootstrapLinearRegression(vector<double>x, vector<double> y, int nboot, double& a, double &b, double& SdOfa, double &SdOfb, double &DevSQ, double &sd, double &RSS, double &Maxd, double &Mind, double& Meand)
{
	int length = x.size();
	if (x.size() != y.size())
	{
		cout << "Error:\tThe size of x does not equal to the size of y in bootstrap linear regression.\n";
		flog.mf_Input("Error:\tThe size of x does not equal to the size of y in bootstrap linear regression.\n");
		flog.mf_Destroy();
		exit(1);
	}
	int irand;
	vector<double> trainx;
	vector<double> trainy;
	vector<double> va;
	vector<double> vb;
	vector<double> vDevSQ;
	vector<double> vsd;
	vector<double> vRSS;
	vector<double> vMaxd;
	vector<double> vMind;
	vector<double> vMeand;
	srand((unsigned)time(NULL));
	for (int i = 0; i < nboot; i++)
	{
		trainx.clear();
		trainy.clear();
		for (int j = 0; j < length; j++)
		{
			irand = rand() % length;
			trainx.push_back(x.at(irand));
			trainy.push_back(y.at(irand));
		}
		mf_linrearRegress(trainx, trainy, a, b, DevSQ, sd, RSS, Maxd, Mind, Meand);
		va.push_back(a);
		vb.push_back(b);
		vDevSQ.push_back(DevSQ);
		vsd.push_back(sd);
		vRSS.push_back(RSS);
		vMaxd.push_back(Maxd);
		vMind.push_back(Mind);
		vMeand.push_back(Meand);
	}

	a = mf_mean(va);
	b = mf_mean(vb);
	DevSQ = mf_mean(vDevSQ);
	sd = mf_mean(vsd);
	RSS = mf_mean(vRSS);
	Maxd = mf_mean(vMaxd);
	Mind = mf_mean(vMind);
	Meand = mf_mean(vMeand);

	double dMeanOfa=0.0;
	double dMeanOfb=0.0;
	for (size_t i = 0; i < va.size(); i++)
	{
		dMeanOfa += va.at(i);
		dMeanOfb += vb.at(i);	
	}
	dMeanOfa = dMeanOfa / va.size();
	dMeanOfb = dMeanOfb / vb.size();

	SdOfa = 0.0;
	SdOfb = 0.0;
	for (size_t i = 0; i < va.size(); i++)
	{
		SdOfa += (va.at(i) - dMeanOfa)*(va.at(i) - dMeanOfa);
		SdOfb += (vb.at(i) - dMeanOfb)*(vb.at(i) - dMeanOfb);
	}
	SdOfa = sqrt(SdOfa / va.size());
	SdOfb = sqrt(SdOfb / vb.size());
}

void BootLineaRegress::mf_Predict(vector<double> intensity, double a, double b, vector<double>& mols)
{
	mols.clear();
	double tempd;
	for (size_t i = 0; i < intensity.size(); i++)
	{
		tempd = intensity.at(i)*a + b;
		mols.push_back(tempd);
	}

}
double BootLineaRegress::mf_Predict(double intensity, double a, double b)
{
	return intensity*a + b;
}
