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
#include"QuantificationParam.h"
#include "funs.h"


extern int UniquePeptidesTrainThreshold;

class CBART
{
public:
	CBART();
	bool mf_BARTIteration( vector<CProtein> &proteins,  CQuantificationParam trainparam, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 3.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);
	bool mf_BARTAutoRegressionRun( vector<CProtein> &proteins, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 3.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);
	bool mf_CalculatePeptidesCVofProteins(vector<CProtein> proteins, double &FirstQ, double &median, double &ThirdQ);
	void mf_xitofile(std::ofstream& fnm, xinfo xi);
	bool mf_ConstructXY(vector<CProtein> proteins, string strExperimentName);
	bool mf_ConstructTestX(vector<CProtein> proteins, string strExperimentName);
	bool mf_ConstructTestXByDigestion(vector<CProtein> proteins);
	bool mf_BartRegressionRun( vector<CProtein>&proteins, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 1.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);
	void mf_CorrectedPeptidesIntensity(vector<CProtein>& proteins);
	void mf_saveRegressionResult( string path);  
	int mf_GetNumberOfTrainProteins(vector<CProtein>& proteins);


	void mf_SelectFeature( vector<tree> trees, size_t burn,map<size_t, size_t> &FeatureCount);
	void mf_saveTestXY( string path);
	void mf_saveTrainXY( string path);

	// for ten-fold cross-validation
	void mf_ChooseBestParameter( vector<CProtein> &proteins, CQuantificationParam trainparam);

	void mf_MakeIndexPartition(size_t SampleNumber,vector<size_t> &IndexPartation); //make index partation for dividing dataset
	void mf_ConstructCVData(  vector<CProtein> proteins, vector<size_t>IndexPartation, size_t fold, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double>&TestyForCV);
	void mf_TrainAndPredictForCV(vector<CProtein> &proteins, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double> &PredictedY, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 3.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);

	double m_dBestBurn;
	double m_dBestAlpha;
	double m_dBestBeta;
	size_t m_iBestNumberOfTrees;
	string m_strQuantiResultPath;
	string m_strSelectedFeaturePath;

private:
	std::vector<double> m_vTrainx;
	std::vector<double> m_vTrainy;
	size_t m_iFeatureNumber;
	size_t m_iTrainSampleNumber;
	double m_dAutoCorrelation;

	double m_dMiny;
	double m_dMaxy;
	sinfo m_cAllys;

	double m_dMeany;
	double	m_dStdy;

	double * m_pTestx;
	std::vector<double> m_vTesty;
	std::vector<double> m_vTestPredicty;
	vector<double> m_vTrainPredicty;
	dinfo dip; //data information for prediction
	size_t m_iTestSampleNumber;
	size_t m_iNumberOfTrainProteins;
	size_t m_iNumberOfTestProteins;

	xinfo xi;
	size_t nc = 100; //100 equally spaced cutpoints from min to max.

	//trees
	std::vector<tree> m_vTrees;

	//prior and mcmc
	pinfo pi;
	dinfo di;

	//for sigma draw
	double m_dRss;  //residual sum of squares
	double m_dRestemp; //a residual

	double m_dPearsonCorrelationUPS2PredictedWithMols;

};