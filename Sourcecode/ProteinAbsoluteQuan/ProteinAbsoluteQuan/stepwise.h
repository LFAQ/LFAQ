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

#include"BasicClass.h"
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
using namespace std;

extern int UniquePeptidesTrainThreshold;

class Cstepwise
{
private:

	double  **m_px, m_dthresh1, m_dthresh2, m_dEps;
	double  *m_pxx;// The average value of indenpendent variables and dependent variables
	double **m_pTrainXForCV;  //Storing X in train dataset for Cross Validation£»
	double * m_pTrainPredictY;
	int *m_psortIndex;
	int *m_pLeaveXIndex;
	double *m_pCoefficient; //regression coefficients
	double *m_pv;//sum of squares of partial regression and residual sum of squares
	double *m_ps;  // the standard deviation of regression coefficients and the estimated SD
	double m_dC;  //multiple correlation coefficients
	double m_dF;  //the value of F-test
	double *m_pye;// m_NumberOfLimitedSamples estimated values of the conditional expectation of independent variables
	double *m_pyr;//m_NumberOfLimitedSamples residuals of observed variables
	double **m_pr;  //normalized coefficient correlation matrix
	int m_iTestSampleNumber;//the number of samples in test dataset for cross validation
	double m_dCorrelation;
	double m_dBestCorrelation;

	double m_dRMSE;
	double m_dMeanY;
	double m_dSdeY;
	int m_iNumberOfFeature;
	double **m_pTestXY;      //Storing X in test dataset for Cross Validation and prediction.
	double *m_pTestPredictY;
	double *m_pNativeTestPredictY; //to save the predicted y before discretization.

	map<string, string> m_mapProteinIDAndSequence;

public:
	int m_iNumberOfTrainProteins;
	int m_iNumberOfTestProteins;
	int m_NumberOfTrainSamples;
	int m_NumberOfLimitedSamples;  //the number of samples limited by unique peptides threshold.
	Cstepwise()
	{
		m_NumberOfTrainSamples = 0;
		m_NumberOfLimitedSamples = 0;
	}

	void mf_StepwiseInit(vector<CProtein> proteins);

	bool mf_Normalize(double **p, int m, int n);// Normalization For X and Y in int **x; 
	bool mf_Normalize(int n);

	void mf_StepwiseRegression( double, double, double); 
	void mf_showXY( );  
	void mf_saveStepwiseXY(string path);

	bool mf_ConstructXY( vector<CProtein> proteins);
	bool mf_ConstructTestX();
	void mf_saveRegressionResult( string path);  

	void mf_Predict(vector<CProtein> &proteins, string regressionResultPath);
	bool mf_loadModel(string m_strRegressionResult);

	bool mf_ConstructTestX( vector<CProtein> proteins);
	void mf_saveStepwiseX( vector<CProtein> proteins, string path);
	void mf_LoadTestXLFAQpep( string LFAQpepPathByOtherMethod);

	void mf_SaveFeatures(string strFeaturePath);
	void mf_Prediction(vector<CProtein> &proteins);  
	void mf_saveTestXY( string path);
	~Cstepwise();
};
