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

// ProteinAbsoluteQuan.cpp : the main file of the program
#include "stdafx.h"
#include"ProteinWorker.h"


int _tmain(int argc, _TCHAR* argv[]) 
{
	clock_t begin, end;
	begin = clock();
	if (argc <= 1)
	{
		cout << "ERROR: Usage: ProteinAbsoluteQuan [parameters file]" << endl;
		cout << "Input parameters file\'s Format:" << endl;
		cout << "xxx.params" << endl;
		cout << "Version -0.1" << endl;
		system("pause");
		return 1;
	}
	string ParaFilePath = Unicode2Multibyte(argv[1]);
	if (!CheckFilePath(ParaFilePath))
	{
		return 2;
	}
	CProteinWorker ProteinWorker;
	string NativePeptidesPath;
	string strProteinLFAQpepPath;
	string strPredictedMolsPath;
	string strResultPath = GetResultPath(ParaFilePath);
	
	if (!flog.mf_Init(strResultPath,ios::app))
	{
		return -1;
	}

	string strPrologue = AddDecorateStar("The beginning of ProteinAbsoluteQuan module");
	cout << strPrologue;
	flog.mf_Input(strPrologue);

	//set parameters
	ProteinWorker.m_Params.mf_setparameters( ParaFilePath);

	//load data
	ProteinWorker.mf_LoadProteinFasta(ProteinWorker.m_Params.m_strProteinFastaPath, ProteinWorker.m_Params.m_eFastaType);
	int iNumberOfExperiments = ProteinWorker.mf_GetExperimentNames(ProteinWorker.m_Params);
	
	ProteinWorker.mf_LoadTraindProteins(ProteinWorker.m_Params);

	for (int i = 0; i < iNumberOfExperiments; i++)
	{

		flog.mf_Input("In experiment " + ProteinWorker.m_vecStrExperiments.at(i) + "\n");
		ProteinWorker.mf_LoadProteins(ProteinWorker.m_Params, i);
		//quantify the proteins  
		ProteinWorker.mf_CalculateProteinIntensity(i);
		// save results
		NativePeptidesPath = ProteinWorker.m_Params.m_strQuantiResultPath \
            + "\\PeptideFeatures" + ProteinWorker.m_vecStrExperiments.at(i) \
            + ".txt";
		ProteinWorker.mf_SaveProteinsNativePeptides( NativePeptidesPath);
		strProteinLFAQpepPath = ProteinWorker.m_Params.m_strProteinLFAQpepPath \
            + ProteinWorker.m_vecStrExperiments.at(i) + ".txt";

		if (ProteinWorker.m_Params.m_bIfCotainStandardProtein)
		{			
			ProteinWorker.mf_PredictMolsByStandProteinsAndLFAQ(i);
			ProteinWorker.mf_PredictMolsByStandProteinsAndiBAQ(i);
			ProteinWorker.mf_PredictMolsByStandProteinsAndTopN(i);

			ProteinWorker.mf_ShowUPS2AnalysisResult(i);
			string strQfactorsOfUPSPeptidesPath = ProteinWorker.m_Params.m_strQuantiResultPath + "\\QfactorsOfUPSPeptides_"+\
				ProteinWorker.m_vecStrExperiments[i]+".txt";
			ProteinWorker.mf_SaveQfactorsOfUPSPeptides(strQfactorsOfUPSPeptidesPath);
		}

		ProteinWorker.mf_saveProteinLFAQpep( strProteinLFAQpepPath);
	}

	string AllProteinLFAQpepPath = ProteinWorker.m_Params.m_strQuantiResultPath \
            +"\\ProteinMergedResults.txt";
	ProteinWorker.mf_MergeProteinLFAQpep( AllProteinLFAQpepPath);
	string strAllRegressionResultPath = ProteinWorker.m_Params.m_strQuantiResultPath + "\\MergedRegressionResults.txt";
	ProteinWorker.mf_MergeRegressionResult( strAllRegressionResultPath);

	end = clock();
		
	int iRunTime = (end - begin) / CLOCKS_PER_SEC;
	string strRunTime = fInt2String(iRunTime);
	cout << "The ProteninAbsoluteQuan module took " << strRunTime<< " seconds" << endl;
	flog.mf_Input("The ProteninAbsoluteQuan module took " + strRunTime + " seconds\n");
	string strEnd = AddDecorateStar("The end of ProteninAbsoluteQuan module");
	cout << strEnd;
	flog.mf_Input(strEnd);
	flog.mf_Destroy();
	return 0;
}

