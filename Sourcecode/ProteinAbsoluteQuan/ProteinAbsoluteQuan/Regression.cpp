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
#include "Regression.h"


void CRegression::mf_StepwiseRegression(vector<CProtein> &proteins, CQuantificationParam trainparam, string strExperimentName)
{
	Cstepwise stepwise;
	stepwise.mf_StepwiseInit(proteins);
	stepwise.mf_ConstructXY(proteins);

#ifdef Debug_Mode
	// Save TrainXY for testing
	string CheckXYpath = "CheckXY.txt";
	stepwise.mf_saveStepwiseXY(CheckXYpath);
#endif

	double f1 = 0.0, f2 = 0.0;
	boost::math::fisher_f_distribution<double> F1_distribution(1, stepwise.m_NumberOfLimitedSamples - 1);
	boost::math::fisher_f_distribution<double> F2_distribution(1, stepwise.m_NumberOfLimitedSamples);
	f1 = boost::math::quantile(F1_distribution, trainparam.m_dAlpha1);
	f2 = boost::math::quantile(F2_distribution, trainparam.m_dAlpha2);
	cout << "\tf1= " << f1 << " \tf2= " << f2 << endl;
	flog.mf_Input("\tf1= "+fDouble2String( f1) + " \tf2= " +fDouble2String(f2)+"\n");

	string strRegressionResultPath = trainparam.m_strRegressionResult + strExperimentName + ".txt";
	stepwise.mf_StepwiseRegression( f1, f2, 1.0e-32); 
	stepwise.mf_saveRegressionResult( strRegressionResultPath);
	stepwise.mf_Predict( proteins, strRegressionResultPath);
	string strSelectedFeaturePath = trainparam.m_strQuantiResultPath + "\\SelectedFeatures" + strExperimentName + ".txt";
	stepwise.mf_SaveFeatures(strSelectedFeaturePath);
}

void CRegression::mf_BARTRegression(const vector<CProtein> &vTrainProteins, vector<CProtein> & vTestProteins, CQuantificationParam trainparam, string strExperimentName)
{

	CBART bartRegression;
	bartRegression.m_strQuantiResultPath = trainparam.m_strQuantiResultPath;
	bartRegression.m_strSelectedFeaturePath = bartRegression.m_strQuantiResultPath + "\\SelectedFeatures" + strExperimentName + ".txt";


	if (trainparam.m_bIfOptimizeParameters)
	{
		// CrossValidation for parameters
		//bartRegression.mf_ChooseBestParameter(vTestProteins, trainparam);
		//bartRegression.mf_ChooseBestParameterByCVCorrelation(proteins, trainparam);
	}
	else
	{
		bartRegression.m_dBestAlpha = trainparam.m_dAlpha;
		bartRegression.m_dBestBeta = trainparam.m_dBeta;
		bartRegression.m_iBestNumberOfTrees = trainparam.m_iNumberOfTrees;
	}

	bartRegression.mf_ConstructXY(vTrainProteins, strExperimentName);

	bartRegression.mf_ConstructTestX(vTestProteins, strExperimentName);
	bartRegression.mf_BartRegressionRun(vTestProteins, 100, 100, bartRegression.m_iBestNumberOfTrees, 1.0, 3.0, 2.0, bartRegression.m_dBestAlpha, bartRegression.m_dBestBeta);
	string strRegressionResultPath = trainparam.m_strRegressionResult + strExperimentName + ".txt";
	bartRegression.mf_saveRegressionResult(strRegressionResultPath);

}

void CRegression::mf_RegressionRun(const vector<CProtein> &vTrainProteins, vector<CProtein> & vTestProteins, CQuantificationParam trainparam, string strExperimentName)
{
	if (trainparam.m_strRegressionMethod == "stepwise")
	{
		mf_StepwiseRegression(vTestProteins, trainparam, strExperimentName);
	}
	else if (trainparam.m_strRegressionMethod == "BART")
	{
		mf_BARTRegression(vTrainProteins, vTestProteins, trainparam, strExperimentName);
	}
}

