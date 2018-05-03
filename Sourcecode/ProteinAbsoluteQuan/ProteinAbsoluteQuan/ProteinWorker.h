/*
#  Copyright(C) 2015-2018 all rights reserved

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

#include<windows.h>
#include"ProteinDigestion.h"
#include"Regression.h"
#include"BootstrapLinearRegression.h"
#include<boost/random.hpp>
#include <chrono>
#include"ExtremeDistribution.h"

// this class is responsible for implementing protein quantification 
class CProteinWorker
{
public:
	CProteinWorker()
	{
		m_bIfiBAQIntensityExist = false; 
		m_iPeptidesNumberOfProteins = 0;		
	}
	
	bool mf_LoadProteinFasta( string strProteinFasterFilePath, FastaType fastatype);
	bool mf_LoadProteins(CQuantificationParam Params, int iExperimentIndex);
	bool mf_LoadTraindProteins(CQuantificationParam Params);
	bool mf_LoadProteinLFAQpep( string path);
	bool mf_LoadUPS2Mols( string path, map<string, double>&mapUPS2Id2Mols);


	bool mf_AnalysisPeptidesCVForSameUPS2Mols( string path, vector<CProtein>& vecProteins, map<string, double>&mapUPS2Id2Mols);

	int mf_GetExperimentNames(CQuantificationParam params);
	void mf_CalculateExperimentalQfactor();
	void mf_CalcExperiQfactorWithMaxIdentiNum();
	void mf_CalcExperiQfactorWithMaxTheoryNum();
	void mf_CalcExperiQfactorWithMaxSequenceLenth();
	void mf_CalcExperiQfactorWithlogMaxIdentiNum();
	void mf_CalcExperiQfactorWithlogMaxTheoryNum();
	void mf_CalcExperiQfactorWithlogMaxSequenceLenth();
	void mf_DetermineTrainSet();

	void mf_CalculateProteinIntensity( int iExperimentIndex);
	void mf_CalculateProteiniBAQ();
	void mf_CalculateProteinLFAQpep( int iExperimentIndex);
	void mf_CorrectProteinPeptides(vector<CProtein>& proteins);
	void mf_PredictMolsByStandAndLFAQ(vector<CProtein>& proteins,double& REsq);
	void mf_CalculateTop3( int iExperimentIndex);

	void mf_SaveProteinsNativePeptides( string path);
	void mf_saveProteinLFAQpep(string ProteinLFAQpepResultPath);
	void mf_SaveMergedProtein(const string& MergedProteinResultPath, const CMergedProteins& mergedProteins);
    void mf_SavePeptidesQfactors(string path);
	void mf_SaveQfactorsOfUPSPeptides(string path);
	void mf_MergeProteinLFAQpep(string ProteinQsoreresultPath);
	void mf_MergeRegressionResult( string RegressionResultPath);
	void mf_ShowUPS2AnalysisResult(int iExperimentIndex);

	
	void mf_PredictMolsByStandProteinsAndLFAQ(int iExperimentIndex);
	void mf_PredictMolsByStandProteinsAndiBAQ(int iExperimentIndex);
	void mf_PredictMolsByStandProteinsAndTopN(int iExperimentIndex);


	void mf_Clear();

	void CalExperQfactorWithMontoKarlo();
	void CalExperQfactorWithEExponentDist();
	void CalExperQfactorWithExtremDist();

	vector<string> m_vecStrExperiments;  // experiment names in order;
	CQuantificationParam m_Params;
	vector<CProtein> m_vecProteins;
	vector<CProtein> m_vTrainProteins;
	map<string, int> m_mapLowSpikedInUPS2;


private:
	int m_iPeptidesNumberOfProteins;  //The number of peptides of proteins loaded
	map<string, string> m_mapProteinsIdAndSequence;

	bool mf_GetProteinInfoAttributeColumns(const map<string, int> &mapAttrtibuteAndcolumns, int &iProteinIDColumnNum,\
		int &iProteinFullNameColumnNum,int &iPeptideSequenceColumnNum,vector<int> &vecIntensitiesColumn,\
		int &iiBAQExistColumn, vector<int> &veciBAQColumns);
	bool mf_GetTrainProteinInfoAttributeColumns(const map<string, int> &mapAttrtibuteAndcolumns, int &iProteinIDColumnNum, \
		int &iProteinFullNameColumnNum, int &iPeptideSequenceColumnNum, int &iCombinedIntensityColumn);
	int m_iNumberOfExperiments;
	bool m_bIfiBAQIntensityExist;

};
