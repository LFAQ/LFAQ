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


#include"stdafx.h"
#include"ProteinWorker.h"


// load information of proteins from ProteinsInfo.txt
bool CProteinWorker::mf_LoadProteins(CQuantificationParam trainParam, int iExperimentIndex)
{
	cout << "Loading proteins from " << trainParam.m_strProteinsPath << "\n";
	flog.mf_Input("Loading proteins from " + trainParam.m_strProteinsPath + "\n");
	//open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, trainParam.m_strProteinsPath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << trainParam.m_strProteinsPath << endl;
		flog.mf_Input("Error:\tCannot open " + trainParam.m_strProteinsPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}

	//temporary variable
	char Buffer[BUFFERLENGTH];
	string strLine;
	char *pstr;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	int iProteinIDColumnNum = 0, iProteinFullNameColumnNum=0, iPeptideSequenceColumnNum = 0;
	int iiBAQExistColumn = 0;
	vector<int> vecIntensitiesColumn;
	vector<int> veciBAQColumns;

	// get the columes of peptides attributes according to the first row
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	strLine = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns, "\t");
	if (!mf_GetProteinInfoAttributeColumns(mapAttrtibuteAndcolumns, iProteinIDColumnNum,\
		iProteinFullNameColumnNum, iPeptideSequenceColumnNum,vecIntensitiesColumn,\
		iiBAQExistColumn, veciBAQColumns))
	{
		flog.mf_Destroy();
		exit(1);
	}

	//string strProteinIDTemp;
	//string strProteinFullNameTemp;
	string strPeptideSequencesTemp;
	vector<string> vecPeptideSequencesTemp;
	vector<string>::iterator vecPeptideSequenceTempIter;
	string strPepIntensitiesTemp;
	vector<double> vecPepIntensitiesTemp;
	vector<double>::iterator vecPepIntensitiesTempIter;
	string striBAQExistTemp;

	vector<string> vecStrTemps;
	CProtein ProteinTemp;
	map<string, string>::iterator mapProteinIDAndSequenceIter;
	map<string, double> mapPeptideAndIntensitys;
	map<string, double>::iterator mapPeptideAndIntensityIter;
	string strWithFlankingRegionTemp;
	map<string, string> mapSequenceAndFlankingRegions;
	map<string, string>::iterator mapSequenceAndFlankingRegionsIter;

	CPeptideAttribute peptideFeaturesTemp;
	PepAttributeWorker PeptideAttributeWorkerTemp;
	string strTemp;
	int begin, end;

	m_vecProteins.clear();
	m_iPeptidesNumberOfProteins = 0;

	fgets(Buffer, BUFFERLENGTH, pFile);
	strLine = Buffer;
	while (!feof(pFile))
	{
		if (strLine == "\n")
		{
			fgets(Buffer, BUFFERLENGTH, pFile);
			strLine = Buffer;
			continue;
		}
		vecStrTemps.clear();
		vecStrTemps = split(strLine, "\t");
		if (vecStrTemps.size() != mapAttrtibuteAndcolumns.size())
		{
			flog.mf_Input("Warning:\tThe columns number of line " + strLine + " do not equal to the columns number of the header line\n");
		}
		ProteinTemp.Clear();

		ProteinTemp.m_strProteinID=vecStrTemps.at(iProteinIDColumnNum);
		

		mapProteinIDAndSequenceIter = m_mapProteinsIdAndSequence.find(ProteinTemp.m_strProteinID);
		if (mapProteinIDAndSequenceIter != m_mapProteinsIdAndSequence.end())
		{
			ProteinTemp.m_strProteinSequence = mapProteinIDAndSequenceIter->second;
		}
		else
		{
			cout << "\tWarning:\tCannot find protein " << ProteinTemp.m_strProteinID << " in fasta file." << endl;
			flog.mf_Input("\tWarning:\tCannot find protein " + ProteinTemp.m_strProteinID + " in fasta file.\n");
			fgets(Buffer, BUFFERLENGTH, pFile);
			strLine = Buffer;
			continue;
		}

		ProteinTemp.m_strProteinFullName=vecStrTemps.at(iProteinFullNameColumnNum);
		strPeptideSequencesTemp = vecStrTemps.at(iPeptideSequenceColumnNum);
		vecPeptideSequencesTemp.clear();
		
		begin = 0;
		end = strPeptideSequencesTemp.find(";", begin);
		while (end != strPeptideSequencesTemp.npos)
		{   //many peptides
			strTemp = strPeptideSequencesTemp.substr(begin, end - begin);
			vecPeptideSequencesTemp.push_back(strTemp);
			begin = end + 1;
			end = strPeptideSequencesTemp.find(";", begin);
		}

		strPepIntensitiesTemp = vecStrTemps.at(vecIntensitiesColumn[iExperimentIndex]);	
		vecPepIntensitiesTemp.clear();
		begin = 0;
		end = strPepIntensitiesTemp.find(";", begin);
		while (end != strPepIntensitiesTemp.npos)
		{   //many peptides
			strTemp = strPepIntensitiesTemp.substr(begin, end - begin);
			vecPepIntensitiesTemp.push_back(atof(strTemp.c_str()));
			begin = end + 1;
			end = strPepIntensitiesTemp.find(";", begin);
		}
		if (vecPeptideSequencesTemp.size() != vecPepIntensitiesTemp.size())
		{
			cout << "The column \"vPeptidesSequence\" does not correspond with \"Intensity\" in line " << strLine << endl;
			flog.mf_Input("The column \"vPeptidesSequence\" does not correspond with \"Intensity\" in line " + strLine + "\n");
			flog.mf_Destroy();
			exit(-1);
		}
		mapPeptideAndIntensitys.clear();
		mapSequenceAndFlankingRegions.clear();
		for (int i = 0; i < vecPeptideSequencesTemp.size(); i++)
		{
			if (vecPepIntensitiesTemp[i] == 0.0)
			{
				continue;
			}
			strWithFlankingRegionTemp = ProteinTemp.mf_GetPeptidesAdjacentSequence(vecPeptideSequencesTemp[i]);
			if (strWithFlankingRegionTemp == "NULL")
			{
				cout << "\tWarning:\tCannot find peptide " << vecPeptideSequencesTemp[i] << "'s adjacent Sequence!\n";
				flog.mf_Input("\tWarning:\tCannot find peptide " + vecPeptideSequencesTemp[i] + "'s adjacent Sequence!\n");
				continue;
			}
			mapSequenceAndFlankingRegions.insert(pair<string, string>(vecPeptideSequencesTemp[i], strWithFlankingRegionTemp));
			mapPeptideAndIntensitys.insert(pair<string, double>(vecPeptideSequencesTemp[i],vecPepIntensitiesTemp[i]));
		}
		if (mapPeptideAndIntensitys.size() == 0){
			fgets(Buffer, BUFFERLENGTH, pFile);
			strLine = Buffer;
			continue;
		}

		m_iPeptidesNumberOfProteins += mapPeptideAndIntensitys.size();;
		mapPeptideAndIntensityIter = mapPeptideAndIntensitys.begin();
		for (; mapPeptideAndIntensityIter != mapPeptideAndIntensitys.end(); mapPeptideAndIntensityIter++)
		{
			ProteinTemp.m_vPeptidesSequences.push_back(mapPeptideAndIntensityIter->first);
			mapSequenceAndFlankingRegionsIter = mapSequenceAndFlankingRegions.find(mapPeptideAndIntensityIter->first);
			assert(mapSequenceAndFlankingRegionsIter != mapSequenceAndFlankingRegions.end());
			ProteinTemp.m_vPeptidesSequencesWithFlankingRegion.push_back(mapSequenceAndFlankingRegionsIter->second);

			peptideFeaturesTemp.Clear();
			peptideFeaturesTemp.m_vecFeatures = PeptideAttributeWorkerTemp.mf_GetAttributeFromSequence(mapPeptideAndIntensityIter->first,\
				mapSequenceAndFlankingRegionsIter->second);
			ProteinTemp.m_vPeptidesFeatures.push_back(peptideFeaturesTemp);
			ProteinTemp.m_vdPeptidesMW.push_back(peptideFeaturesTemp.m_vecFeatures.at(2));  //The third attribute is MW
			ProteinTemp.m_vPeptidesNativeIntensity.push_back(mapPeptideAndIntensityIter->second);
			ProteinTemp.m_vPeptidesCorrectedIntensity.push_back(mapPeptideAndIntensityIter->second); //assign the same value with m_vPeptidesNativeIntensity for iterration
			ProteinTemp.m_vPeptideExperimentalQFactors.push_back(0.0);
		}

		ProteinTemp.m_iPeptidesNumber = ProteinTemp.m_vPeptidesSequences.size();

		striBAQExistTemp=vecStrTemps[iiBAQExistColumn];
		m_bIfiBAQIntensityExist = atoi(striBAQExistTemp.c_str());
		if (striBAQExistTemp == "1")
		{
			ProteinTemp.m_dMaxQuantiBAQ = atof(vecStrTemps[veciBAQColumns[iExperimentIndex]].c_str());
		}
		else
		{
			ProteinTemp.m_dMaxQuantiBAQ = 0.0;
		}

		m_vecProteins.push_back(ProteinTemp);

		fgets(Buffer, BUFFERLENGTH, pFile);
		strLine = Buffer;
	}

	if (m_vecProteins.size() == 0)
	{
		cout << "Error:\tPlease check the fasta file or the input file!" << endl;
		flog.mf_Input("Error:\tPlease check the fasta file or the input file!\n");
		flog.mf_Destroy();
		exit(1);
	}
	cout << "\t" << m_vecProteins.size() << " proteins with " << m_iPeptidesNumberOfProteins << " peptides in experiment "<<m_vecStrExperiments.at(iExperimentIndex) <<" loaded"<< endl;
	flog.mf_Input("\t" +fSize_t2String(m_vecProteins.size())+ " proteins with " +fInt2String(m_iPeptidesNumberOfProteins)+ " peptides in experiment " + m_vecStrExperiments.at(iExperimentIndex)+" loaded\n");
	flog.mf_Input("Done\n");
	return 1;
}

bool CProteinWorker::mf_LoadTraindProteins(CQuantificationParam Params)
{
	cout << "Loading proteins from " << Params.m_strProteinsPath << "\n";
	flog.mf_Input("Loading proteins from " + Params.m_strProteinsPath + "\n");
	//open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, Params.m_strProteinsPath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << Params.m_strProteinsPath << endl;
		flog.mf_Input("Error:\tCannot open " + Params.m_strProteinsPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}

	//temporary variable
	char Buffer[BUFFERLENGTH];
	string strLine;
	char *pstr;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	int iProteinIDColumnNum = 0, iProteinFullNameColumnNum = 0, iPeptideSequenceColumnNum = 0;
	int iCombinedIntensitiesColumn;

	// get the columes of peptides attributes according to the first row
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	strLine = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns, "\t");
	if (!mf_GetTrainProteinInfoAttributeColumns(mapAttrtibuteAndcolumns, iProteinIDColumnNum, \
		iProteinFullNameColumnNum, iPeptideSequenceColumnNum, iCombinedIntensitiesColumn))
	{
		flog.mf_Destroy();
		exit(1);
	}
	string strPeptideSequencesTemp;
	vector<string> vecPeptideSequencesTemp;
	vector<string>::iterator vecPeptideSequenceTempIter;
	string strPepIntensitiesTemp;
	vector<double> vecPepIntensitiesTemp;
	vector<double>::iterator vecPepIntensitiesTempIter;

	vector<string> vecStrTemps;
	CProtein ProteinTemp;
	map<string, string>::iterator mapProteinIDAndSequenceIter;
	map<string, double> mapPeptideAndIntensitys;
	map<string, double>::iterator mapPeptideAndIntensityIter;
	string strWithFlankingRegionTemp;
	map<string, string> mapSequenceAndFlankingRegions;
	map<string, string>::iterator mapSequenceAndFlankingRegionsIter;
	CPeptideAttribute peptideFeaturesTemp;
	PepAttributeWorker PeptideAttributeWorkerTemp;
	string strTemp;
	int begin, end;

	m_vTrainProteins.clear();
	m_iPeptidesNumberOfProteins = 0;
	fgets(Buffer, BUFFERLENGTH, pFile);
	strLine = Buffer;
	while (!feof(pFile))
	{
		if (strLine == "\n")
		{
			fgets(Buffer, BUFFERLENGTH, pFile);
			strLine = Buffer;
			continue;
		}
		vecStrTemps.clear();
		vecStrTemps = split(strLine, "\t");
		if (vecStrTemps.size() != mapAttrtibuteAndcolumns.size())
		{
			flog.mf_Input("Warning:\tThe columns number of line " + strLine + " do not equal to the columns number of the header line\n");
		}
		ProteinTemp.Clear();

		ProteinTemp.m_strProteinID = vecStrTemps.at(iProteinIDColumnNum);

		mapProteinIDAndSequenceIter = m_mapProteinsIdAndSequence.find(ProteinTemp.m_strProteinID);
		if (mapProteinIDAndSequenceIter != m_mapProteinsIdAndSequence.end())
		{
			ProteinTemp.m_strProteinSequence = mapProteinIDAndSequenceIter->second;
		}
		else
		{
			cout << "\tWarning:\tCannot find protein " << ProteinTemp.m_strProteinID << " in fasta file." << endl;
			flog.mf_Input("\tWarning:\tCannot find protein " + ProteinTemp.m_strProteinID + " in fasta file.\n");
			fgets(Buffer, BUFFERLENGTH, pFile);
			strLine = Buffer;
			continue;
		}

		ProteinTemp.m_strProteinFullName = vecStrTemps.at(iProteinFullNameColumnNum);
		strPeptideSequencesTemp = vecStrTemps.at(iPeptideSequenceColumnNum);
		vecPeptideSequencesTemp.clear();
		begin = 0;
		end = strPeptideSequencesTemp.find(";", begin);
		while (end != strPeptideSequencesTemp.npos)
		{   //many peptides
			strTemp = strPeptideSequencesTemp.substr(begin, end - begin);
			vecPeptideSequencesTemp.push_back(strTemp);
			begin = end + 1;
			end = strPeptideSequencesTemp.find(";", begin);
		}
		strPepIntensitiesTemp = vecStrTemps.at(iCombinedIntensitiesColumn);
		vecPepIntensitiesTemp.clear();
		begin = 0;
		end = strPepIntensitiesTemp.find(";", begin);
		while (end != strPepIntensitiesTemp.npos)
		{   //many peptides
			strTemp = strPepIntensitiesTemp.substr(begin, end - begin);
			vecPepIntensitiesTemp.push_back(atof(strTemp.c_str()));
			begin = end + 1;
			end = strPepIntensitiesTemp.find(";", begin);
		}
		assert(vecPeptideSequencesTemp.size() == vecPepIntensitiesTemp.size());

		mapSequenceAndFlankingRegions.clear();
		mapPeptideAndIntensitys.clear();
		for (int i = 0; i < vecPeptideSequencesTemp.size(); i++)
		{
			if (vecPepIntensitiesTemp[i] == 0.0)
			{
				continue;
			}
			strWithFlankingRegionTemp = ProteinTemp.mf_GetPeptidesAdjacentSequence(vecPeptideSequencesTemp[i]);
			if (strWithFlankingRegionTemp == "NULL")
			{
				cout << "\tWarning:\tCannot find peptide " << vecPeptideSequencesTemp[i] << "'s adjacent Sequence!\n";
				flog.mf_Input("\tWarning:\tCannot find peptide " + vecPeptideSequencesTemp[i] + "'s adjacent Sequence!\n");
				continue;
			}
			mapSequenceAndFlankingRegions.insert(pair<string, string>(vecPeptideSequencesTemp[i], strWithFlankingRegionTemp));
			mapPeptideAndIntensitys.insert(pair<string, double>(vecPeptideSequencesTemp[i], vecPepIntensitiesTemp[i]));
		}
		if (mapPeptideAndIntensitys.size() == 0){
			fgets(Buffer, BUFFERLENGTH, pFile);
			strLine = Buffer;
			continue;
		}

		m_iPeptidesNumberOfProteins += mapPeptideAndIntensitys.size();;
		mapPeptideAndIntensityIter = mapPeptideAndIntensitys.begin();
		for (; mapPeptideAndIntensityIter != mapPeptideAndIntensitys.end(); mapPeptideAndIntensityIter++)
		{
			ProteinTemp.m_vPeptidesSequences.push_back(mapPeptideAndIntensityIter->first);
			mapSequenceAndFlankingRegionsIter = mapSequenceAndFlankingRegions.find(mapPeptideAndIntensityIter->first);
			assert(mapSequenceAndFlankingRegionsIter != mapSequenceAndFlankingRegions.end());
			ProteinTemp.m_vPeptidesSequencesWithFlankingRegion.push_back(mapSequenceAndFlankingRegionsIter->second);
			peptideFeaturesTemp.Clear();
			peptideFeaturesTemp.m_vecFeatures = PeptideAttributeWorkerTemp.mf_GetAttributeFromSequence(mapPeptideAndIntensityIter->first, \
				mapSequenceAndFlankingRegionsIter->second);
			ProteinTemp.m_vPeptidesFeatures.push_back(peptideFeaturesTemp);
			ProteinTemp.m_vdPeptidesMW.push_back(peptideFeaturesTemp.m_vecFeatures.at(2));  //The third attribute is MW
			ProteinTemp.m_vPeptidesNativeIntensity.push_back(mapPeptideAndIntensityIter->second);
			ProteinTemp.m_vPeptidesCorrectedIntensity.push_back(mapPeptideAndIntensityIter->second); //assign the same value with m_vPeptidesNativeIntensity for iterration
			ProteinTemp.m_vPeptideExperimentalQFactors.push_back(0.0);
		}

		ProteinTemp.m_iPeptidesNumber = ProteinTemp.m_vPeptidesSequences.size();

		m_vTrainProteins.push_back(ProteinTemp);

		fgets(Buffer, BUFFERLENGTH, pFile);
		strLine = Buffer;
	}

	if (m_vTrainProteins.size() == 0)
	{
		cout << "Error:\tPlease check the fasta file or the input file!" << endl;
		flog.mf_Input("Error:\tPlease check the fasta file or the input file!\n");
		flog.mf_Destroy();
		exit(1);
	}
	cout << "\t" << m_vTrainProteins.size() << " proteins with " << m_iPeptidesNumberOfProteins << " peptides were loaded" << endl;
	flog.mf_Input("\t" + fSize_t2String(m_vTrainProteins.size()) + " proteins with " + fInt2String(m_iPeptidesNumberOfProteins) + " peptides were loaded\n");
	flog.mf_Input("Done\n");
	return 1;

}

bool CProteinWorker::mf_GetProteinInfoAttributeColumns(const map<string, int> &mapAttrtibuteAndcolumns, int &iProteinIDColumnNum, int &iProteinFullNameColumnNum,
	int &iPeptideSequenceColumnNum,	vector<int> &vecIntensitiesColumn, int &iiBAQExistColumn, vector<int> &veciBAQColumns)
{
	
	map<string, int>::const_iterator mapAttrtibuteAndcolumnsIter;
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("ProteinID");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinIDColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"ProteinID\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"ProteinID\" in the ProteinsInfo.txt file.\n");
		return false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Majority protein IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinFullNameColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"Majority protein IDs\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"Majority protein IDs\" in the ProteinsInfo.txt file.\n");
		return false;
	}
	
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("vPeptidesSequence");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptideSequenceColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"vPeptidesSequence\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"vPeptidesSequence\" in the ProteinsInfo.txt file.\n");
		return false;
	}

	vecIntensitiesColumn.clear();
	for (int i = 0; i < m_vecStrExperiments.size(); i++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Intensity "+m_vecStrExperiments[i]);
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			vecIntensitiesColumn.push_back(mapAttrtibuteAndcolumnsIter->second);
		}
		else
		{
			cout << "Error:\tCannot find column \"Intensity "+m_vecStrExperiments[i]+"\" in the ProteinsInfo.txt file." << endl;
			flog.mf_Input("Error:\tCannot find column \"Intensity " + m_vecStrExperiments[i] + "\" in the ProteinsInfo.txt file.\n");
			return false;
		}
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("bIfiBAQExist");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iiBAQExistColumn = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"bIfiBAQExist\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"bIfiBAQExist\" in the ProteinsInfo.txt file.\n");
		return false;
	}

	veciBAQColumns.clear();
	for (int i = 0; i < m_vecStrExperiments.size(); i++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("iBAQ " + m_vecStrExperiments[i]);
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			veciBAQColumns.push_back(mapAttrtibuteAndcolumnsIter->second);
		}
	}

	return true;

}

void CProteinWorker::mf_DetermineTrainSet()
{
	vector<CProtein>::iterator  ProteinIter;
	int iTrainPeptidesNum = 0;
	for (ProteinIter = m_vTrainProteins.begin(); ProteinIter != m_vTrainProteins.end(); ProteinIter++)
	{
		ProteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(ProteinIter->m_vPeptidesNativeIntensity);
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			if (ProteinIter->m_vPeptidesNativeIntensity.at(ProteinIter->m_iPeptidesMaxIntensityIndex) != 0.0)
			{
				ProteinIter->m_bIfInTrainSet = true;
				iTrainPeptidesNum += ProteinIter->m_iPeptidesNumber;
			}
		}
	}
	cout << "There are " << iTrainPeptidesNum << " peptides in the train set.\n";
	flog.mf_Input("There are " + fInt2String(iTrainPeptidesNum) + " peptides in the train set.\n");
}
bool CProteinWorker::mf_GetTrainProteinInfoAttributeColumns(const map<string, int> &mapAttrtibuteAndcolumns, int &iProteinIDColumnNum, \
	int &iProteinFullNameColumnNum, int &iPeptideSequenceColumnNum, int &iCombinedIntensityColumn)
{
	map<string, int>::const_iterator mapAttrtibuteAndcolumnsIter;
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("ProteinID");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinIDColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"ProteinID\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"ProteinID\" in the ProteinsInfo.txt file.\n");
		return false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Majority protein IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinFullNameColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"Majority protein IDs\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"Majority protein IDs\" in the ProteinsInfo.txt file.\n");
		return false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("vPeptidesSequence");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptideSequenceColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"vPeptidesSequence\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"vPeptidesSequence\" in the ProteinsInfo.txt file.\n");
		return false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Intensity");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iCombinedIntensityColumn = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"Intensity\" in the ProteinsInfo.txt file." << endl;
		flog.mf_Input("Error:\tCannot find column \"Intensity\" in the ProteinsInfo.txt file.\n");
		return false;
	}
	return true;

}

// Method One
void  CProteinWorker::mf_CalculateExperimentalQfactor()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	vector<double> vPeptideExperitalQfactors;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex);
				vPeptideExperitalQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));					
			}
		}
	}


#ifdef Debug_Mode
	string strAllPeptideExperitalQfactors = m_Params.m_strQuantiResultPath + "\\AllPeptideExperitalQfactors.txt";
	SaveVector(strAllPeptideExperitalQfactors, vPeptideExperitalQfactors);
#endif
}

//Method Two
void CProteinWorker::CalExperQfactorWithMontoKarlo()
{
	vector<CProtein>::iterator proteinIter;
	int iIterations = 1000;
	int iDistributionLength;
	vector<double> vecdDistribution;
	vector<double> vecdMax;
	double dMaxInOneDist;
	double dTemp;
	vector<double> vecdQfactorOfMaxpeptide;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			if (m_Params.strDistributionLenghType == "IdentifiedNum")
			{
				iDistributionLength = 100 * proteinIter->m_vPeptidesSequences.size();
			}
			else if (m_Params.strDistributionLenghType == "IntheoryNum")
			{
				iDistributionLength = 100 * proteinIter->m_iNumberOfTheoreticEnzyme;
			}
			else
			{
				cout << "Error\t Cannot parse the parameter DistributionLenghType\n";
			}
			vecdMax.clear();
			for (int i = 0; i < iIterations; i++)
			{
				int seed = std::chrono::system_clock::now().time_since_epoch().count();
				std::default_random_engine generator(seed);
				boost::exponential_distribution<double> expdis(1.0 / m_Params.dMeanQfactor); // exponential_distribution(RealType lambda = 1);
				dMaxInOneDist = 0.0;
				for (int i = 0; i < iDistributionLength; i++)
				{
					dTemp = expdis(generator);
					if (dTemp > dMaxInOneDist)
						dMaxInOneDist = dTemp;
				}
				vecdMax.push_back((dMaxInOneDist));
			}
			vecdQfactorOfMaxpeptide.push_back(Average(vecdMax));
		}
	}
	dTemp = GetMaxValue(vecdQfactorOfMaxpeptide);
	for (int i = 0; i < vecdQfactorOfMaxpeptide.size(); i++)
	{
		vecdQfactorOfMaxpeptide.at(i) = vecdQfactorOfMaxpeptide.at(i) / dTemp;
	}
	int iIndex = 0;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{

			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = proteinIter->m_vPeptidesNativeIntensity.at(k)\
					/ proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex);
				proteinIter->m_vPeptideExperimentalQFactors.at(k) =\
					proteinIter->m_vPeptideExperimentalQFactors.at(k)*vecdQfactorOfMaxpeptide.at(iIndex);
			}
			iIndex++;
		}
	}
}

void CProteinWorker::CalExperQfactorWithEExponentDist()
{
	map<int, double> mapCumulSum;
	mapCumulSum.insert(pair<int, double>(1, 1.0));
	vector<CProtein>::iterator proteinIter;
	int iDistributionLength;
	vector<double> vecdQfactorOfMaxpeptide;
	double dTempSum = 0.0;
	for (proteinIter = m_vTrainProteins.begin(); proteinIter != m_vTrainProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			iDistributionLength = proteinIter->m_strProteinSequence.size();
			iDistributionLength = iDistributionLength / m_Params.dDistributionCorrectFactor;

			if (mapCumulSum.size() < iDistributionLength)
			{
				dTempSum = mapCumulSum[mapCumulSum.size()];
				for (int i = mapCumulSum.size() + 1; i <= iDistributionLength; i++)
				{
					dTempSum += (double)(1.0 / i);
					mapCumulSum.insert(pair<int, double>(i, dTempSum));
				}
			}
			vecdQfactorOfMaxpeptide.push_back(mapCumulSum[iDistributionLength]);
		}
	}
	double dTemp;
	dTemp = GetMaxValue(vecdQfactorOfMaxpeptide);
	for (int i = 0; i < vecdQfactorOfMaxpeptide.size(); i++)
	{
		vecdQfactorOfMaxpeptide.at(i) = vecdQfactorOfMaxpeptide.at(i) / dTemp;
	}
	int iIndex = 0;
	for (proteinIter = m_vTrainProteins.begin(); proteinIter != m_vTrainProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at\
					(proteinIter->m_iPeptidesMaxIntensityIndex);
				proteinIter->m_vPeptideExperimentalQFactors.at(k) =\
					proteinIter->m_vPeptideExperimentalQFactors.at(k)*vecdQfactorOfMaxpeptide.at(iIndex);
			}
			iIndex++;
		}
	}
}

void CProteinWorker::CalExperQfactorWithExtremDist()
{
	vector<CProtein>::iterator proteinIter;
	int iIterations = 1000;
	int iDistributionLength;
	vector<double> vecdDistribution;
	vector<double> vecdMax;
	double dMaxInOneDist;
	double dTemp;
	vector<double> vecdQfactorOfMaxpeptide;
	int seed = 0;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			if (m_Params.strDistributionLenghType == "IdentifiedNum")
			{
				iDistributionLength = 100 * proteinIter->m_vPeptidesSequences.size();
			}
			else if (m_Params.strDistributionLenghType == "IntheoryNum")
			{
				iDistributionLength = 100 * proteinIter->m_iNumberOfTheoreticEnzyme;
			}
			else
			{
				cout << "Error\t Cannot parse the parameter DistributionLenghType\n";
			}

			vecdMax.clear();
			seed++;
			for (int i = 0; i < iIterations; i++)
			{
				seed = seed + i;
				std::default_random_engine generator(seed);
				extreme_value_distribution_MinCase<double> ExtremeDist(m_Params.dMeanQfactor, m_Params.dStdQfactor);
				dMaxInOneDist = 0.0;
				for (int i = 0; i < iDistributionLength; i++)
				{
					dTemp = ExtremeDist(generator);
					if (dTemp > dMaxInOneDist)
						dMaxInOneDist = dTemp;
				}
				vecdMax.push_back(exp(dMaxInOneDist));
			}
			seed = seed - iIterations;
			vecdQfactorOfMaxpeptide.push_back(Average(vecdMax));
		}
	}
	dTemp = GetMaxValue(vecdQfactorOfMaxpeptide);
	for (int i = 0; i < vecdQfactorOfMaxpeptide.size(); i++)
	{
		vecdQfactorOfMaxpeptide.at(i) = vecdQfactorOfMaxpeptide.at(i) / dTemp;
	}
	int iIndex = 0;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex);
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k)*vecdQfactorOfMaxpeptide.at(iIndex);

			}
			iIndex++;
		}
	}
}
// Method Three
void CProteinWorker::mf_CalcExperiQfactorWithMaxIdentiNum()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dProteinIntensities;
	vector<double> vecdPeptideQfactors;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			dProteinIntensities = 0.0;
			for (int i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
			{
				dProteinIntensities +=proteinIter->m_vPeptidesNativeIntensity.at(i);
			}

			proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex) = \
				(proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex) / dProteinIntensities)*proteinIter->m_vPeptidesSequences.size();
			vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex));
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				if (k != proteinIter->m_iPeptidesMaxIntensityIndex)
				{
					proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
						proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)\
						*proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex);
					vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));
				}
			}

		}
	}

	// normalization by divide by max Qfactor
	double dMaxPeptideQfactor = 0.0;
	dMaxPeptideQfactor = GetMaxValue(vecdPeptideQfactors);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k) / dMaxPeptideQfactor;
			}
		}
	}

}
void CProteinWorker::mf_CalcExperiQfactorWithMaxTheoryNum()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dProteinIntensities;
	vector<double> vecdPeptideQfactors;

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			dProteinIntensities = 0.0;
			for (int i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
			{
				dProteinIntensities += proteinIter->m_vPeptidesNativeIntensity.at(i);
			}

			if (proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex) != 0.0)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex) = \
					(proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex) / dProteinIntensities)*proteinIter->m_iNumberOfTheoreticEnzyme;
				vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex));

				for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
				{
					if (k != proteinIter->m_iPeptidesMaxIntensityIndex)
					{
						proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
							proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)\
							*proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex);
						vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));

					}
				}
			}
		}
	}

	// normalization by divide by max Qfactor
	double dMaxPeptideQfactor = 0.0;
	dMaxPeptideQfactor = GetMaxValue(vecdPeptideQfactors);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k) / dMaxPeptideQfactor;
			}
		}
	}
	
}
void CProteinWorker::mf_CalcExperiQfactorWithMaxSequenceLenth()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dProteinIntensities;
	vector<double> vecdPeptideQfactors;

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			dProteinIntensities = 0.0;
			for (int i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
			{
				dProteinIntensities += proteinIter->m_vPeptidesNativeIntensity.at(i);
			}

			proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex) = \
				(proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex) / dProteinIntensities)*proteinIter->m_strProteinSequence.size();
			vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex));

			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				if (k != proteinIter->m_iPeptidesMaxIntensityIndex)
				{
					proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
						proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)\
						*proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex);
					vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));

				}
			}
		}
	}

	// normalization by divide by max Qfactor
	double dMaxPeptideQfactor = 0.0;
	dMaxPeptideQfactor = GetMaxValue(vecdPeptideQfactors);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k) / dMaxPeptideQfactor;
			}
		}
	}
}
void CProteinWorker::mf_CalcExperiQfactorWithlogMaxIdentiNum()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dProteinIntensities;
	vector<double> vecdPeptideQfactors;

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			dProteinIntensities = 0.0;
			for (int i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
			{
				dProteinIntensities += proteinIter->m_vPeptidesNativeIntensity.at(i);
			}

			if (proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex) != 0.0)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex) = \
					log10(proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)) / log10(dProteinIntensities/proteinIter->m_vPeptidesSequences.size());
				vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex));

				for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
				{
					if (k != proteinIter->m_iPeptidesMaxIntensityIndex)
					{
						proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
							proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)\
							*proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex);
						vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));

					}
				}
			}
		}
	}

	// normalization by divide by max Qfactor
	double dMaxPeptideQfactor = 0.0;
	dMaxPeptideQfactor = GetMaxValue(vecdPeptideQfactors);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k) / dMaxPeptideQfactor;
			}
		}
	}
}
void CProteinWorker::mf_CalcExperiQfactorWithlogMaxTheoryNum()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dProteinIntensities;	
	vector<double> vecdPeptideQfactors;

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{
			dProteinIntensities = 0.0;
			for (int i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
			{
				dProteinIntensities += proteinIter->m_vPeptidesNativeIntensity.at(i);
			}

			proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex) = \
				log10(proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)) / log10(dProteinIntensities / proteinIter->m_iNumberOfTheoreticEnzyme);
			vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex));

			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				if (k != proteinIter->m_iPeptidesMaxIntensityIndex)
				{
					proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
						proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)\
						*proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex);
					vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));

				}
			}
		}
	}

	// normalization by divide by max Qfactor
	double dMaxPeptideQfactor = 0.0;
	dMaxPeptideQfactor = GetMaxValue(vecdPeptideQfactors);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k) / dMaxPeptideQfactor;
			}
		}
	}
}
void CProteinWorker::mf_CalcExperiQfactorWithlogMaxSequenceLenth()
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dProteinIntensities;
	vector<double> vecdPeptideQfactors;

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		proteinIter->m_iPeptidesMaxIntensityIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
		if (proteinIter->m_bIfInTrainSet == true)
		{

			dProteinIntensities = 0.0;
			for (int i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
			{
				dProteinIntensities += proteinIter->m_vPeptidesNativeIntensity.at(i);
			}

			proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex) = \
				log(proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)) / log(dProteinIntensities / proteinIter->m_strProteinSequence.size());
			vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex));

			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				if (k != proteinIter->m_iPeptidesMaxIntensityIndex)
				{
					proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
						proteinIter->m_vPeptidesNativeIntensity.at(k) / proteinIter->m_vPeptidesNativeIntensity.at(proteinIter->m_iPeptidesMaxIntensityIndex)\
						*proteinIter->m_vPeptideExperimentalQFactors.at(proteinIter->m_iPeptidesMaxIntensityIndex);
					vecdPeptideQfactors.push_back(proteinIter->m_vPeptideExperimentalQFactors.at(k));

				}
			}
		}
	}
	// normalization by divide by max Qfactor
	double dMaxPeptideQfactor = 0.0;
	dMaxPeptideQfactor = GetMaxValue(vecdPeptideQfactors);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet == true)
		{
			for (size_t k = 0; k < proteinIter->m_vPeptideExperimentalQFactors.size(); k++)
			{
				proteinIter->m_vPeptideExperimentalQFactors.at(k) = \
					proteinIter->m_vPeptideExperimentalQFactors.at(k) / dMaxPeptideQfactor;
			}
		}
	}
}

//load protein fasta 
bool CProteinWorker::mf_LoadProteinFasta( string strProteinFasterFilePath, FastaType fastatype)
{
	char Buffer[BUFFERLENGTH];
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, strProteinFasterFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << strProteinFasterFilePath << endl;
		flog.mf_Input("Error:\tCannot open " +strProteinFasterFilePath+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	char *pstr;
	string strIDTemp;
	string strSequenceTemp;
	string strTemp;
	fgets(Buffer, BUFFERLENGTH, pFile);
	while ((Buffer[0] == '\n') && (!feof(pFile)))  //allowing empty lines
	{
		fgets(Buffer, BUFFERLENGTH, pFile);
	}


	std::match_results<std::string::const_iterator> Matchresult;
	bool valid;
	string strFirstRow;
	pstr = Buffer;
	if (Buffer[0] == '>')  //get the protein id according to the regular expression fastatype;
	{
		strFirstRow = Buffer;
		valid = std::regex_search(strFirstRow, Matchresult, fastatype);
		if (valid == true)
		{
			strIDTemp = Matchresult[1];
		}
		else
		{
			cout << "Error:\tCannot parse the fasta file by the regular expression.\n";
			flog.mf_Input("Error:\tCannot parse the fasta file by the regular expression.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	while (fgets(Buffer, BUFFERLENGTH, pFile))
	{
		if (Buffer[0] == '\0')   //allow the blank lines
			continue;
		pstr = Buffer;

		if (Buffer[0] == '>')  //get the protein id according to the regular expression fastatype;
		{
			m_mapProteinsIdAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));
			strSequenceTemp.clear();
			strFirstRow = Buffer;
			valid = std::regex_search(strFirstRow, Matchresult, fastatype);
			if (valid == true)
			{
				strIDTemp = Matchresult[1];
			}
			else
			{
				cout << "Error:\tCannot parse the fasta file by the regular expression.\n";
				flog.mf_Input("Error:\tCannot parse the fasta file by the regular expression.\n");
				flog.mf_Destroy();
				exit(1);
			}
		}
		else // get the sequence of the protein id
		{
			strTemp.clear();
			strTemp = pstr;
			strSequenceTemp = strSequenceTemp + strTemp.substr(0, strTemp.size() - 1);
		}
	}
	m_mapProteinsIdAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));

	fclose(pFile);
	return 1;
}

//calculate the iBAQ and LFAQpep by calling mf_CalculateProteiniBAQ and mf_CalculateProteinLFAQpep
void CProteinWorker::mf_CalculateProteinIntensity(int iExperimentIndex)
{
	if (m_Params.m_bIfCalculateTop3)
	{
		mf_CalculateTop3(iExperimentIndex);
	}
	mf_CalculateProteiniBAQ();
	mf_CalculateProteinLFAQpep( iExperimentIndex);
}

// recalculate the iBAQ according to the protein' peptides intensity;
void  CProteinWorker::mf_CalculateProteiniBAQ()
{
	if (m_Params.m_bIfCalculateiBAQ)
	{
		cout << "Calculate proteins' abundance using iBAQ method\n";
		flog.mf_Input("Calculate proteins' abundance using iBAQ method\n");
	}
	CProteinDigestion proteinDigestion(m_Params);
	for (int i = 0; i < m_vecProteins.size(); i++)
	{
		for (size_t j = 0; j < m_vecProteins.at(i).m_vPeptidesNativeIntensity.size(); j++)
		{
			m_vecProteins.at(i).m_dPeptidesIntensitySum += m_vecProteins.at(i).m_vPeptidesNativeIntensity.at(j);
		}
		m_vecProteins.at(i).m_iPeptidesNumber = m_vecProteins.at(i).m_vPeptidesSequences.size();
		// enzyme number
		m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme = proteinDigestion.mf_CalculateOptDigestionNumber(m_vecProteins.at(i).m_strProteinSequence);

		if (m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme <= 0)
			m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme=1;
		m_vecProteins.at(i).m_dReCalculateiBAQ = m_vecProteins.at(i).m_dPeptidesIntensitySum / m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme;

	}
	
}

// recalculate the LFAQpep according to the protein' peptides intensity;
void CProteinWorker::mf_CalculateProteinLFAQpep(int iExperimentIndex)
{
	cout << "Calculate proteins' abundance using LFAQ method\n";
	flog.mf_Input("Calculate proteins' abundance using LFAQ method\n");

	mf_DetermineTrainSet();
	if (m_Params.strMaxPeptideQfactorType == "MaxQfactorOne")
	{
		mf_CalculateExperimentalQfactor();
	}
	else if (m_Params.strMaxPeptideQfactorType == "ExponDistMonteCarlo")
	{
		CalExperQfactorWithMontoKarlo();
	}
	else if (m_Params.strMaxPeptideQfactorType == "ExponDistExpect")
	{
		CalExperQfactorWithEExponentDist();
	}
	else if (m_Params.strMaxPeptideQfactorType == "ExtremDistMonteCarlo")
	{
		CalExperQfactorWithExtremDist();
	}
	else if (m_Params.strMaxPeptideQfactorType == "MaxAverageAndIdentiNum")
	{
		mf_CalcExperiQfactorWithMaxIdentiNum();
	}
	else if (m_Params.strMaxPeptideQfactorType == "MaxAverageAndTheoryNum")
	{
		mf_CalcExperiQfactorWithMaxTheoryNum();
	}
	else if (m_Params.strMaxPeptideQfactorType == "MaxAverageAndSequenceLength")
	{
		mf_CalcExperiQfactorWithMaxSequenceLenth();
	}
	else if (m_Params.strMaxPeptideQfactorType == "logMaxAverageAndIdentiNum")
	{
		mf_CalcExperiQfactorWithlogMaxIdentiNum();
	}
	else if (m_Params.strMaxPeptideQfactorType == "logMaxAverageAndTheoryNum")
	{
		mf_CalcExperiQfactorWithlogMaxTheoryNum();
	}
	else if (m_Params.strMaxPeptideQfactorType == "logMaxAverageAndSequenceLength")
	{
		mf_CalcExperiQfactorWithlogMaxSequenceLenth();
	}
	else
	{
		cout << "Error:\tCannot parse parameter MaxPeptideQfactorType \n";
		flog.mf_Input("Error:\tCannot parse parameter MaxPeptideQfactorType \n");
		flog.mf_Destroy();
		exit(-1);
	}

	CRegression regression;
	regression.mf_RegressionRun(m_vTrainProteins, m_vecProteins, m_Params, m_vecStrExperiments.at(iExperimentIndex));

	cout << "Correct the peptides' native intensities by Qfactor\n";
	flog.mf_Input("Correct the peptides' native intensities by Qfactor\n");
	
	mf_CorrectProteinPeptides(m_vecProteins);
}
void CProteinWorker::mf_CorrectProteinPeptides(vector<CProtein>& proteins)
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	vector<double>::iterator PeptideQfactorsIter;
	vector<double> CorrectedPeptidesIntensity;
	vector<double> NativePeptidesIntensity;
	double dPeptideIntensityTemp;
	int i = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		i = 0;
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			proteinIter->m_bIfCalculateLFAQpep = true;
			proteinIter->m_bIfCalculateLFAQpro = true;
			if (proteinIter->m_iPeptidesNumber > UniquePeptidesCorrectThreshold)
			{
				proteinIter->m_dProteinLFAQpep = 0.0;
				proteinIter->m_dProteinLFAQpro = 0.0;
				proteinIter->m_vPeptidesCorrectedIntensity.clear();
				for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++)
				{
					dPeptideIntensityTemp = *peptidesIter;
					proteinIter->m_dProteinLFAQpro += dPeptideIntensityTemp;
					NativePeptidesIntensity.push_back(dPeptideIntensityTemp/proteinIter->m_vdPeptidesMW[i]);

					if (proteinIter->m_vPeptideQfactors.at(i) != 0.0)
						dPeptideIntensityTemp = dPeptideIntensityTemp / proteinIter->m_vPeptideQfactors.at(i);
					else
						dPeptideIntensityTemp = 0.0;
					i++;
					CorrectedPeptidesIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_vPeptidesCorrectedIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_dProteinLFAQpep += dPeptideIntensityTemp;
				}
			}
			else
			{ //not enough peptides for correction
				proteinIter->m_dProteinLFAQpep = 0.0;
				proteinIter->m_dProteinLFAQpro = 0.0;
				proteinIter->m_vPeptidesCorrectedIntensity.clear();
				for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++)
				{
					dPeptideIntensityTemp = *peptidesIter;
					proteinIter->m_dProteinLFAQpro += dPeptideIntensityTemp;
					NativePeptidesIntensity.push_back(dPeptideIntensityTemp / proteinIter->m_vdPeptidesMW[i]);
					CorrectedPeptidesIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_vPeptidesCorrectedIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_dProteinLFAQpep += dPeptideIntensityTemp;
					i++;
				}
			}
			proteinIter->m_dCVNative = CalculatelogCV(NativePeptidesIntensity);
			proteinIter->m_dCVAfterCorrected = CalculatelogCV(CorrectedPeptidesIntensity);
			NativePeptidesIntensity.clear();
			CorrectedPeptidesIntensity.clear();
			if (proteinIter->m_iPeptidesNumber != 0)
			{
				proteinIter->m_dProteinLFAQpep = proteinIter->m_dProteinLFAQpep / proteinIter->m_iPeptidesNumber;

				PeptideQfactorsIter = proteinIter->m_vPeptideQfactors.begin();
				proteinIter->m_dProteinQfactor = 0.0;
				for (; PeptideQfactorsIter != proteinIter->m_vPeptideQfactors.end(); PeptideQfactorsIter++)
				{
					proteinIter->m_dProteinQfactor += *PeptideQfactorsIter;
				}

				proteinIter->m_dProteinLFAQpro = proteinIter->m_dProteinLFAQpro / proteinIter->m_dProteinQfactor;
			}
			else
			{
				proteinIter->m_dProteinLFAQpep = 0.0;
				proteinIter->m_dProteinLFAQpro = 0.0;
				cout << "\tWarning:\tProtein " << proteinIter->m_strProteinID << " do not have unique peptides.\n";
				cout << "\tSo we make its LFAQ 0.0.\n";
				flog.mf_Input("\tWarning:\tProtein " + proteinIter->m_strProteinID + " do not have unique peptides. \n");
				flog.mf_Input("So we make its LFAQ 0.0.\n");

			}
		}
	}
}

void  CProteinWorker::mf_PredictMolsByStandAndLFAQ(vector<CProtein>& proteins, double& REsq)
{
	map<string, double> mapUPS2Id2Mols;
	map<string, double> mapUPS2iBAQs;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	map<string, double>::iterator UPS2iBAQIter;
	vector<CProtein>::iterator proteinIter;
	string strTemp;

	CProteinWorker proteinworker;
	cout << "Predict proteins' abundance by the standard proteins\n";
	flog.mf_Input("Predict proteins' abundance by the standard proteins\n");
	proteinworker.mf_LoadUPS2Mols(m_Params.m_strStandardProteinsFilePath, mapUPS2Id2Mols);

	double dPearsonCorr = 0.0;
	vector<double> veclogLFAQpepOfStandProteins;//the LFAQpep of identified standard proteins
	vector<double> vecUPS2logMolsIdentified; // the spiked-in amount of identified standard proteins;
	vector<double> vecPredictedLogMolsOfStandProteinsByLFAQpep;

	int iNumberofPeptidesOfStandProteins = 0;
	string::size_type stStart, stPosition, stEnd;
	vector<string> vecProteinsTemp;
	bool bIfOthersOnlyContam = true;
	map<string, int>::iterator mapLowSpikedInUPS2Iter;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateLFAQpep == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stPosition = strTemp.find(m_Params.m_strIdentifierOfStandPro);
			vecProteinsTemp = split(proteinIter->m_strProteinFullName, ";");
			if (stPosition != strTemp.npos)
			{
				bIfOthersOnlyContam = true;
				if (m_Params.m_bIfExistContanProtein == true && vecProteinsTemp.size() > 1)
				{
					for (size_t i = 1; i < vecProteinsTemp.size(); i++)
					{
						if (vecProteinsTemp[i].find(m_Params.m_strContaminantfix) == vecProteinsTemp[i].npos)
						{
							bIfOthersOnlyContam = false;
						}
					}
				}
				if ((proteinIter->m_strProteinFullName.find(";") == \
					proteinIter->m_strProteinFullName.npos) || (bIfOthersOnlyContam == true))
				{ //
					Ups2MolsIter = mapUPS2Id2Mols.begin();
					for (; Ups2MolsIter != mapUPS2Id2Mols.end(); Ups2MolsIter++)
					{
						if (strTemp.find(Ups2MolsIter->first) != strTemp.npos)
						{
							if (proteinIter->m_dProteinLFAQpep != 0.0)
							{
								iNumberofPeptidesOfStandProteins += proteinIter->m_vPeptidesSequences.size();
								if (Ups2MolsIter->second == 0.0)
								{
									vecUPS2logMolsIdentified.push_back(0.0);
								}
								else
								{
									vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
								}
								if (proteinIter->m_dProteinLFAQpep == 0.0)
								{
									veclogLFAQpepOfStandProteins.push_back(0.0);
								}
								else
								{
									veclogLFAQpepOfStandProteins.push_back(log10(proteinIter->m_dProteinLFAQpep));
								}
								
							} // end if 0.0
							break;
						}
					} // for mapUPS2Id2Mols
					if (Ups2MolsIter == mapUPS2Id2Mols.end())
					{
						mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
						if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
						{
							cout << "Warning:\tCannot find protein " << strTemp << " in StandardProteins.txt\n";
							flog.mf_Input("Warning:\tCannot find protein " + strTemp + " in StandardProteins.txt\n");
						}
					}
				}
			}
		}
	}

	cout << "\tLoaded " << vecUPS2logMolsIdentified.size() << " standard proteins " << " with " << \
		iNumberofPeptidesOfStandProteins << " peptides\n";
	flog.mf_Input("\tLoaded " + fSize_t2String(vecUPS2logMolsIdentified.size()) + " standard proteins with " + \
		fInt2String(iNumberofPeptidesOfStandProteins) + " peptides\n");

	BootLineaRegress bootLinear;
	double a, b;
	double DevSQ;
	double sd;
	double RSS;
	double Maxd;
	double Mind;
	double Meand;
	double SdOfa;
	double SdOfb;

	cout << "\tPredict protein amount by bootstrap linear regression\n";
	flog.mf_Input("\tPredict protein amount by bootstrap linear regression\n");

	// Predict protein amount by LFAQpep
	cout << "\t\tPredict protein amount by LFAQ" << endl;
	flog.mf_Input("\t\tPredict protein amount by LFAQ\n");

	// No normalized
	bootLinear.mf_BootstrapLinearRegression(veclogLFAQpepOfStandProteins, vecUPS2logMolsIdentified, 10000, a, b, SdOfa, SdOfb, DevSQ, sd, RSS, Maxd, Mind, Meand);

	double R2;
	bootLinear.mf_Predict(veclogLFAQpepOfStandProteins, a, b, vecPredictedLogMolsOfStandProteinsByLFAQpep);
	R2 = spearsonCorrelation(vecUPS2logMolsIdentified, vecPredictedLogMolsOfStandProteinsByLFAQpep);
	if (vecUPS2logMolsIdentified.size() != vecPredictedLogMolsOfStandProteinsByLFAQpep.size())
	{
		cout << "Error:\tThe number of standard proteins does not equal to the number of the predicted proteins.\n";
		flog.mf_Input("Error:\tThe number of standard proteins does not equal to the number of the predicted proteins.\n");
		exit(0);
	}
	REsq = 0.0;
	for (int i = 0; i < vecUPS2logMolsIdentified.size(); i++)
	{
		REsq += (vecUPS2logMolsIdentified.at(i) - vecPredictedLogMolsOfStandProteinsByLFAQpep.at(i))\
			*(vecUPS2logMolsIdentified.at(i) - vecPredictedLogMolsOfStandProteinsByLFAQpep.at(i));
	}
	cout << "\t\t\ta= " << a << endl;
	cout << "\t\t\tb= " << b << endl;
	cout << "\t\t\tSd of a= " << SdOfa << endl;
	cout << "\t\t\tSd of b= " << SdOfb << endl;
	cout << "\t\t\tMean Standard Error= " << sd << endl;
	cout << "\t\t\tAutocorrelation R2= " << R2*R2 << endl;
	flog.mf_Input("\t\t\ta= " + fDouble2String(a) + "\n");
	flog.mf_Input("\t\t\tb= " + fDouble2String(b) + "\n");
	flog.mf_Input("\t\t\tSd of a= " + fDouble2String(SdOfa) + "\n");
	flog.mf_Input("\t\t\tSd of b= " + fDouble2String(SdOfb) + "\n");
	flog.mf_Input("\t\t\tMean Standard Error= " + fDouble2String(sd) + "\n");
	flog.mf_Input("\t\t\tAutocorrelation R2= " + fDouble2String(R2*R2) + "\n");

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateLFAQpep == true)
		{
			if (proteinIter->m_dProteinLFAQpep == 0.0)
			{
				proteinIter->m_dPredictedMolOfLFAQpep =0.0;
			}
			else
			{
				proteinIter->m_dPredictedMolOfLFAQpep = pow(10, bootLinear.mf_Predict(log10(proteinIter->m_dProteinLFAQpep), a, b));
			}
			
		}
	}


}

void CProteinWorker::mf_CalculateTop3(int iExperimentIndex)
{
	cout << "Calculate proteins' abundance using Top3 method\n";
	flog.mf_Input("Calculate proteins' abundance using Top3 method\n");
	vector<double> vPeptideIntensitiesTemp;
	vector<size_t> vPeptideIndexesTemp;
	for (int i = 0; i < m_vecProteins.size(); i++)
	{
		if (m_vecProteins.at(i).m_vPeptidesNativeIntensity.size() == 0)
		{
			m_vecProteins.at(i).m_dProteinIntensityTop3 = 0.0;
			continue;
		}
		if (m_vecProteins.at(i).m_vPeptidesNativeIntensity.size() == 1)
		{
			m_vecProteins.at(i).m_dProteinIntensityTop3 = m_vecProteins.at(i).m_vPeptidesNativeIntensity.at(0);
			continue;
		}
		else if (m_vecProteins.at(i).m_vPeptidesNativeIntensity.size() == 2)
		{
			m_vecProteins.at(i).m_dProteinIntensityTop3 += m_vecProteins.at(i).m_vPeptidesNativeIntensity.at(0);
			m_vecProteins.at(i).m_dProteinIntensityTop3 += m_vecProteins.at(i).m_vPeptidesNativeIntensity.at(1);
			m_vecProteins.at(i).m_dProteinIntensityTop3 = m_vecProteins.at(i).m_dProteinIntensityTop3 / 2.0;
			continue;
		}
		else
		{
			vPeptideIntensitiesTemp.clear();
			vPeptideIndexesTemp.clear();
			for (int j = 0; j < m_vecProteins.at(i).m_vPeptidesNativeIntensity.size(); j++)
			{
				vPeptideIntensitiesTemp.push_back(m_vecProteins.at(i).m_vPeptidesNativeIntensity.at(j));
				vPeptideIndexesTemp.push_back(j);
			}
			quickSort(vPeptideIntensitiesTemp, vPeptideIndexesTemp, 0, vPeptideIntensitiesTemp.size() - 1);
			m_vecProteins.at(i).m_dProteinIntensityTop3 += vPeptideIntensitiesTemp.at(0);
			m_vecProteins.at(i).m_dProteinIntensityTop3 += vPeptideIntensitiesTemp.at(1);
			m_vecProteins.at(i).m_dProteinIntensityTop3 += vPeptideIntensitiesTemp.at(2);
			m_vecProteins.at(i).m_dProteinIntensityTop3 = m_vecProteins.at(i).m_dProteinIntensityTop3 / 3;
		}

	}
}

void CProteinWorker::mf_SaveProteinsNativePeptides(string path)
{
	ofstream ofile(path.c_str());
	if (!ofile)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " + path+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "Peptides' features saved to " << path << endl;
		flog.mf_Input("Peptides' features saved to " + path+"\n");
	}
		
	ofile << "PeptidesSequence\tProteinID\tPeptideQFactor\tpeptideOriginalIntensity\tPeptideFeatures\n";
	vector<CProtein>::iterator  ProteinIter;
	int i = 0;
	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
	{
		i = 0;

		for (i = 0; i < ProteinIter->m_iPeptidesNumber; i++)
		{		
			ofile << ProteinIter->m_vPeptidesSequences.at(i) << "\t" ;
			ofile << ProteinIter->m_strProteinID << "\t";
			ofile << ProteinIter->m_vPeptideQfactors.at(i) << "\t";
			ofile << ProteinIter->m_vPeptidesNativeIntensity.at(i) << "\t";
			ofile << ProteinIter->m_vPeptidesFeatures.at(i) << endl;
		}

	}
	ofile.close();
}

void CProteinWorker::mf_saveProteinLFAQpep(string ProteinLFAQpepResultPath)
{
	ofstream ofile;
	ofile.open(ProteinLFAQpepResultPath.c_str());
	if (!ofile)
	{
		cout << "Error:\tCannot open " << ProteinLFAQpepResultPath << endl;
		flog.mf_Input("Error:\tCannot open " + ProteinLFAQpepResultPath+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "\tSaved " << m_vecProteins.size() << " proteins' LFAQ into " << ProteinLFAQpepResultPath << endl;
		flog.mf_Input("\tSaved "+ fSize_t2String(m_vecProteins.size()) + " proteins' LFAQ into " + ProteinLFAQpepResultPath +"\n");
	}

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	ofile << "Protein IDs\tMajority protein IDs\tNumberOfUniquePeptides\t";
	ofile << "LFAQ\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		ofile << "iBAQ\t";
	}

	if (m_Params.m_bIfCalculateTop3)
	{
		ofile << "Top3\t";
	}
	ofile << "CVOriginal\tCVAfterCorrected\t";
	ofile << "NumberOfTheoreticEnzyme\tPredictedMol(LFAQ)\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		ofile << "PredictedMol(iBAQ)\t";
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		ofile << "PredictedMol(Top3)\t";
	}
	ofile << "Coverage\tPeptideSequences\tPeptideOriginalIntensities\tPeptideQfactors\t";
	ofile << "PeptideCorrectedIntensities\tPeptideMWs\n";

	// Calculate the sequence coverage 
	vector<int> vecAppearedSequence;
	int iBegin, iEnd;
	int i, j;
	double dCoverageTemp;

    ofile << setprecision(6);
    ofile.setf(ios::fixed);
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		ofile << proteinIter->m_strProteinFullName << "\t" << proteinIter->m_strProteinID << "\t";
		ofile << proteinIter->m_iPeptidesNumber << "\t";
		ofile << proteinIter->m_dProteinLFAQpep << "\t";
		if (m_Params.m_bIfCalculateiBAQ)
		{
			if (m_bIfiBAQIntensityExist)
			{
				ofile << proteinIter->m_dMaxQuantiBAQ << "\t";
			}
			else
			{
				ofile << proteinIter->m_dReCalculateiBAQ << "\t";
			}
		}
		if (m_Params.m_bIfCalculateTop3)
		{
			ofile << proteinIter->m_dProteinIntensityTop3 << "\t";
		}
		ofile<< proteinIter->m_dCVNative << "\t" << proteinIter->m_dCVAfterCorrected << "\t";
		ofile << proteinIter->m_iNumberOfTheoreticEnzyme << "\t";		
		ofile << proteinIter->m_dPredictedMolOfLFAQpep << "\t";
		if (m_Params.m_bIfCalculateiBAQ)
		{
			ofile << proteinIter->m_dPredictedMolOfiBAQ << "\t";
		}
		if (m_Params.m_bIfCalculateTop3)
		{
			ofile << proteinIter->m_dPredictedMolOfTop3 << "\t";
		}
		// Calculate the sequence coverage 
		vecAppearedSequence.clear();
		for (i = 0; i < proteinIter->m_strProteinSequence.length(); i++)
			vecAppearedSequence.push_back(0);
		for (i = 0; i < proteinIter->m_vPeptidesSequences.size(); i++)
		{
			iBegin = proteinIter->m_strProteinSequence.find(proteinIter->m_vPeptidesSequences.at(i));
			if (iBegin == string::npos)
			{
				cout << "Warning:\tCannot find peptide " << proteinIter->m_vPeptidesSequences.at(i) << \
					" in protein " << proteinIter->m_strProteinID << "'s sequence\n";
				flog.mf_Input("Warning:\tCannot find peptide " + proteinIter->m_vPeptidesSequences.at(i)+\
					+ " in protein " + proteinIter->m_strProteinID + "'s sequence\n");
				break;
			}
			iEnd = iBegin + proteinIter->m_vPeptidesSequences.at(i).size() - 1;
			for (j = iBegin; j <= iEnd; j++)
			{
				vecAppearedSequence.at(j) = 1;
			}
		}
		dCoverageTemp = 0.0;
		for (i = 0; i < vecAppearedSequence.size(); i++)
		{
			dCoverageTemp += vecAppearedSequence.at(i);
		}
		dCoverageTemp = dCoverageTemp / proteinIter->m_strProteinSequence.size();
		ofile << dCoverageTemp << "\t";
		if (proteinIter->m_vPeptidesSequences.size() == 0)
		{
			ofile << "\t";
			ofile << "\t";
			if (proteinIter->m_bIfCalculateLFAQpep)
			{
				ofile << "\t";
			}
			ofile << "\n";
		}
		else
		{
			for (i = 0; i < proteinIter->m_vPeptidesSequences.size() - 1; i++)
			{
				ofile << proteinIter->m_vPeptidesSequences.at(i) << ";";
			}
			ofile << proteinIter->m_vPeptidesSequences.at(i) << "\t";

			for (i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size() - 1; i++)
			{
				ofile << proteinIter->m_vPeptidesNativeIntensity.at(i) << ";";
			}
			ofile << proteinIter->m_vPeptidesNativeIntensity.at(i) << "\t";
			if (proteinIter->m_bIfCalculateLFAQpep)
			{
				for (i = 0; i < proteinIter->m_vPeptideQfactors.size() - 1; i++)
				{
					ofile << proteinIter->m_vPeptideQfactors.at(i) << ";";
				}
				ofile << proteinIter->m_vPeptideQfactors.at(i) << "\t";

				for (i = 0; i < proteinIter->m_vPeptidesCorrectedIntensity.size() - 1; i++)
				{
					ofile << proteinIter->m_vPeptidesCorrectedIntensity.at(i) << ";";
				}
				ofile << proteinIter->m_vPeptidesCorrectedIntensity.at(i) << "\t";
			}
			else
			{
				ofile << "\t";
			}
		}
		if (proteinIter->m_vdPeptidesMW.size()>0)
		{
			for (i = 0; i < proteinIter->m_vdPeptidesMW.size() - 1; i++)
			{
				ofile << proteinIter->m_vdPeptidesMW.at(i) << ";";
			}
			ofile << proteinIter->m_vdPeptidesMW.at(i) << "\n";
		}
		else
		{
			ofile << "\n";
		}
	}
	ofile.close();

}


void CProteinWorker::mf_SavePeptidesQfactors(string path){
    ofstream ofile;
    ofile.open(path.c_str());
    if (!ofile)
    {
        cout << "Error:\tCannot open " << path << endl;
        flog.mf_Input("Error:\tCannot open " + path + "\n");
        flog.mf_Destroy();
        exit(1);
    }

    vector<CProtein>::iterator proteinIter;
    int iPeptideNumber = 0;

    for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++){
        for (int i = 0; i < proteinIter->m_vPeptidesSequences.size(); i++)
        {
            ofile << proteinIter->m_vPeptidesSequences[i] << "\t";
            ofile << proteinIter->m_strProteinID << "\t";
            ofile << proteinIter->m_vPeptideQfactors[i] << endl;
            iPeptideNumber++;
        }
    }


    cout << "\tSaved " << iPeptideNumber << " peptides' Qfactors into " << path << endl;
    flog.mf_Input("\tSaved " + fSize_t2String(iPeptideNumber) + \
        " peptides' Qfactors  into " + path + "\n");
    ofile.close();
}
void CProteinWorker::mf_SaveQfactorsOfUPSPeptides(string path)
{
	map<string, double> mapUPS2Id2Mols;
	mf_LoadUPS2Mols(m_Params.m_strStandardProteinsFilePath, mapUPS2Id2Mols);

	ofstream ofile;
	ofile.open(path.c_str());
	if (!ofile)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " + path + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "\tSaved UPS2 peptides' Q-factor into " << path << endl;
		flog.mf_Input("\tSaved UPS2 peptides' Qfactor into " + path + "\n");
	}
	ofile << "Peptide Sequence\tProtein ID\tQ-factor\n";

	int iBegin, iEnd = 0;
	string strUPS2ProteinIDTemp;
	map<string, double>::iterator Ups2MolsIter;
	for (int i = 0; i < m_vecProteins.size(); i++)
	{
		if (m_vecProteins[i].m_strProteinID.find(m_Params.m_strIdentifierOfStandPro) != \
			m_vecProteins[i].m_strProteinID.npos)
		{
			iBegin = m_vecProteins[i].m_strProteinID.find("|");
			iEnd = m_vecProteins[i].m_strProteinID.find("|", iBegin + 1);
			strUPS2ProteinIDTemp = m_vecProteins[i].m_strProteinID.substr(iBegin + 1, iEnd - iBegin - 1);

			Ups2MolsIter = mapUPS2Id2Mols.find(strUPS2ProteinIDTemp);
			if (Ups2MolsIter == mapUPS2Id2Mols.end())
			{
				cout << "Can not find UPS2 protein " << strUPS2ProteinIDTemp << endl;
				continue;
			}
			assert(m_vecProteins[i].m_vPeptidesSequences.size() == m_vecProteins[i].m_vPeptidesNativeIntensity.size());
			for (int j = 0; j < m_vecProteins[i].m_vPeptidesNativeIntensity.size(); j++)
			{
				ofile << m_vecProteins[i].m_vPeptidesSequences[j] << "\t";
				ofile << m_vecProteins[i].m_strProteinID << "\t";
				ofile << m_vecProteins[i].m_vPeptidesNativeIntensity[j] / Ups2MolsIter->second << endl;;
			}

		}
		
	}
	ofile.close();

}
bool CProteinWorker::mf_LoadProteinLFAQpep( string path)
{

	//open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, path.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " +path+"\n");
		flog.mf_Destroy();
		exit(1);
		return false;
	}

	//temporary variable
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	int icolumns = 0;
	int count = 0;

	int iMajorProteinIDcolumnNum = 0;
	int iProteinIDColumnNum = 0;
	int iiBAQColumnNum = 0;
	int iLFAQpepColumnNum = 0;
	int iTop3ColumnNum = 0;
	int iPredictedMolOfLFAQpepColumnNum = 0;
	int iPredictedMolOfiBAQNum = 0;
	int iPredictedMolOfTop3Num = 0;
	int iCorrectedPepIntensitiesColumnNum = 0;

	bool bMajorProteinIDcolumn = false;
	bool bProteinIDColumnNum = false;
	bool biBAQColumn = false;
	bool bLFAQpepColumn = false;
	bool bTop3Column = false;
	bool bPredictedMolOfLFAQpepColumn = false;
	bool bPredictedMolOfiBAQColumn = false;
	bool bPredictedMolOfTop3Column = false;
	bool bCorrectedPepIntensitiesColumn = false;
	int iAttributeNumbers = 9;

	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Protein IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iMajorProteinIDcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bMajorProteinIDcolumn = true;
		count++;
	}
	else
	{
		bMajorProteinIDcolumn = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Majority protein IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinIDColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bProteinIDColumnNum = true;
		count++;
	}
	else
	{
		bProteinIDColumnNum = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("iBAQ");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iiBAQColumnNum = mapAttrtibuteAndcolumnsIter->second;
		biBAQColumn = true;
		count++;
	}
	else
	{
		biBAQColumn = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("LFAQ");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iLFAQpepColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bLFAQpepColumn = true;
		count++;
	}
	else
	{
		bLFAQpepColumn = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Top3");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iTop3ColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bTop3Column = true;
		count++;
	}
	else
	{
		bTop3Column = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("PredictedMol(LFAQ)");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPredictedMolOfLFAQpepColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bPredictedMolOfLFAQpepColumn = true;
		count++;
	}
	else
	{
		bPredictedMolOfLFAQpepColumn = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("PredictedMol(iBAQ)");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPredictedMolOfiBAQNum = mapAttrtibuteAndcolumnsIter->second;
		bPredictedMolOfiBAQColumn = true;
		count++;
	}
	else
	{
		bPredictedMolOfiBAQColumn = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("PredictedMol(Top3)");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPredictedMolOfTop3Num = mapAttrtibuteAndcolumnsIter->second;
		bPredictedMolOfTop3Column = true;
		count++;
	}
	else
	{
		bPredictedMolOfTop3Column = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("PeptideCorrectedIntensities");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iCorrectedPepIntensitiesColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bCorrectedPepIntensitiesColumn = true;
		count++;
	}
	else
	{
		bCorrectedPepIntensitiesColumn = false;
	}

	if (!m_Params.m_bIfCalculateiBAQ)
	{
		iAttributeNumbers-=2;
	}
	if (!m_Params.m_bIfCalculateTop3)
	{
		iAttributeNumbers-=2;
	}
	if (count < iAttributeNumbers)
	{
		if (!bMajorProteinIDcolumn)
		{
			cout << "Error:\tCannot find the column \"Protein IDs\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Protein IDs\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bProteinIDColumnNum)
		{
			cout << "Error:\tCannot find the column \"Majority protein IDs\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Majority protein IDs\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}

		if (m_Params.m_bIfCalculateiBAQ&&!biBAQColumn)
		{
			cout << "Error:\tCannot find the column \"iBAQ\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"iBAQ\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bLFAQpepColumn)
		{
			cout << "Error:\tCannot find the column \"LFAQpep\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"LFAQpep\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (m_Params.m_bIfCalculateTop3&&!bTop3Column)
		{
			cout << "Error:\tCannot find the column \"Top3\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Top3\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bPredictedMolOfLFAQpepColumn)
		{
			cout << "Error:\tCannot find the column \"PredictedMol(LFAQ)\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"PredictedMol(LFAQ)\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (m_Params.m_bIfCalculateTop3&&!bPredictedMolOfTop3Column)
		{
			cout << "Error:\tCannot find the column \"PredictedMol(Top3)\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"PredictedMol(Top3)\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (m_Params.m_bIfCalculateiBAQ&&!bPredictedMolOfiBAQColumn)
		{
			cout << "Error:\tCannot find the column \"PredictedMol(iBAQ)\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"PredictedMol(iBAQ)\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bCorrectedPepIntensitiesColumn)
		{
			cout << "Error:\tCannot find the column \"PeptidesCorrectedIntensity\" in the ProteinResults.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"PeptidesCorrectedIntensity\" in the ProteinResults.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	icolumns = 0;
	count = 0;
	CProtein proteinTemp;
	string strPeptidesIntensitiesTemp;
	int iBegin, iEnd;
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;

	while (!feof(pFile))
	{
		proteinTemp.Clear();
		count = 0;
		icolumns = 0;
		while (count < iAttributeNumbers)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				pstr1 = strstr(pstr, "\n");
				if (pstr1 != NULL)
				{
					*pstr1 = '\0';
				}
				else
				{
					cout << "Error:\tThe format of proteinResults*.txt  is wrong.\n";
					flog.mf_Input("Error:\tThe format of proteinResults*.txt  is wrong.\n");
					flog.mf_Destroy();
					exit(1);
				}
			}
			if (icolumns == iMajorProteinIDcolumnNum)
			{
				proteinTemp.m_strProteinFullName = pstr;  //Majority protein IDs
				count++;
			}
			if (icolumns == iProteinIDColumnNum)
			{
				proteinTemp.m_strProteinID = pstr;  //protein IDs
				count++;
			}
			if (m_Params.m_bIfCalculateiBAQ&&icolumns == iiBAQColumnNum)
			{
				proteinTemp.m_dMaxQuantiBAQ = atof(pstr);  //iBAQ
				count++;
			}
			if (icolumns == iLFAQpepColumnNum)
			{
				proteinTemp.m_dProteinLFAQpep = atof(pstr);  //LFAQpep
				count++;
			}
			if (m_Params.m_bIfCalculateTop3&&icolumns == iTop3ColumnNum)
			{
				proteinTemp.m_dProteinIntensityTop3 = atof(pstr);  //Top3
				count++;
			}
			if (icolumns == iPredictedMolOfLFAQpepColumnNum)
			{
				proteinTemp.m_dPredictedMolOfLFAQpep = atof(pstr);  //PredictedMolOfLFAQpep
				count++;
			}
			if (m_Params.m_bIfCalculateTop3&&icolumns == iPredictedMolOfTop3Num)
			{
				proteinTemp.m_dPredictedMolOfTop3 = atof(pstr);  //PredictedMolOfTop3
				count++;
			}
			if (m_Params.m_bIfCalculateiBAQ&&icolumns == iPredictedMolOfiBAQNum)
			{
				proteinTemp.m_dPredictedMolOfiBAQ = atof(pstr);  //PredictedMolOfiBAQ
				count++;
			}
			if (icolumns == iCorrectedPepIntensitiesColumnNum)
			{
				strPeptidesIntensitiesTemp = pstr; //PeptidesCorrectedIntensity
				iBegin = 0;
				iEnd = strPeptidesIntensitiesTemp.find(";", iBegin);
				while (iEnd != strPeptidesIntensitiesTemp.npos)
				{
					proteinTemp.m_vPeptidesCorrectedIntensity.push_back(atof((strPeptidesIntensitiesTemp.substr(iBegin, iEnd - iBegin)).c_str()));
					iBegin = iEnd + 1;
					iEnd = strPeptidesIntensitiesTemp.find(";", iBegin);
				}
				proteinTemp.m_vPeptidesCorrectedIntensity.push_back(atof((strPeptidesIntensitiesTemp.substr(iBegin, strPeptidesIntensitiesTemp.size() - iBegin)).c_str()));

				count++;
			}
			pstr = pstr1 + 1;
			icolumns++;
			
		}
		m_vecProteins.push_back(proteinTemp);
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
}

// merge all ProteinLFAQpep in different experiment;
void CProteinWorker::mf_MergeProteinLFAQpep(string MergedProteinResultPath)
{
	cout << "Merging all proteins' LFAQ from different experiments\n";
	flog.mf_Input("Merging all proteins' LFAQ from different experiments\n");
	string strProteinLFAQpepPath;

	// temp variables
	string strIDTemp;
	string strTemp;
	CMergedProteins MergedProteins;
	CMergedProtein mergedProteinTemp;
	vector<double> dvecPeptidesIntensitiesTemp;
	map<string, int> mapProteinNamesAndIndex; // to memorize if the protein appear or not, if yes, the value is the index in MergedProteins
	map<string, int>::iterator mapProteinNamesAndIndexIter;

	dvecPeptidesIntensitiesTemp.clear();
	mapProteinNamesAndIndex.clear();
	// assuming that every LFAQpep result is sorted, our result is OK;
	for (int i = 0; i < m_vecStrExperiments.size(); i++)
	{
		m_vecProteins.clear();
		strProteinLFAQpepPath = m_Params.m_strProteinLFAQpepPath + m_vecStrExperiments.at(i) + ".txt";
		mf_LoadProteinLFAQpep(strProteinLFAQpepPath);
		if (MergedProteins.m_vecMergedProteins.size() == 0) // first experiment
		{
			for (int j = 0; j < m_vecProteins.size(); j++)
			{
				mergedProteinTemp.Clear();
				mergedProteinTemp.m_strPrtoteinName = m_vecProteins[j].m_strProteinID;
				mergedProteinTemp.m_strSubsetProteins = m_vecProteins[j].m_strProteinFullName;
				mergedProteinTemp.m_vecExperiments.push_back(m_vecStrExperiments.at(i));
				mergedProteinTemp.m_vecLFAQpepOfExperiments.push_back(m_vecProteins[j].m_dProteinLFAQpep);
				mergedProteinTemp.m_vecPredictedMolsByLFAQpepOfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfLFAQpep);
				if (m_Params.m_bIfCalculateiBAQ)
				{
					mergedProteinTemp.m_vecMaxQuantiBAQOfExperiments.push_back(m_vecProteins[j].m_dMaxQuantiBAQ);
					mergedProteinTemp.m_vecPredictedMolsByiBAQOfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfiBAQ);
				}
				if (m_Params.m_bIfCalculateTop3)
				{
					mergedProteinTemp.m_vecTop3OfExperiments.push_back(m_vecProteins[j].m_dProteinIntensityTop3);
					mergedProteinTemp.m_vecPredictedMolsByTop3OfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfTop3);
				}
				dvecPeptidesIntensitiesTemp.clear();
				for (int k = 0; k < m_vecProteins[j].m_vPeptidesCorrectedIntensity.size(); k++)
				{
					dvecPeptidesIntensitiesTemp.push_back(m_vecProteins[j].m_vPeptidesCorrectedIntensity[k]);
				}
				mergedProteinTemp.m_mapExperimentsAndPeptidesIntensity.insert(pair<string, vector<double>>(m_vecStrExperiments[i], dvecPeptidesIntensitiesTemp));
				MergedProteins.m_vecMergedProteins.push_back(mergedProteinTemp);
				mapProteinNamesAndIndex.insert(pair<string, int>(m_vecProteins[j].m_strProteinID, j));
			}

		}
		else
		{  // not first experiment
			for (int j = 0; j < m_vecProteins.size(); j++)
			{
				mapProteinNamesAndIndexIter = mapProteinNamesAndIndex.find(m_vecProteins.at(j).m_strProteinID);
				if (mapProteinNamesAndIndexIter != mapProteinNamesAndIndex.end())
				{ // not new protein
					// make up missing experiment
					int iBegin = 0;
					for (; MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecExperiments.back() != m_vecStrExperiments[iBegin]; iBegin++){};
					for (; iBegin < i - 1; iBegin++)
					{
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecLFAQpepOfExperiments.push_back(0.0);
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecPredictedMolsByLFAQpepOfExperiments.push_back(0.0);
						if (m_Params.m_bIfCalculateiBAQ)
						{
							MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecMaxQuantiBAQOfExperiments.push_back(0.0);
							MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecPredictedMolsByiBAQOfExperiments.push_back(0.0);
						}
						if (m_Params.m_bIfCalculateTop3)
						{
							MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecTop3OfExperiments.push_back(0.0);
							MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecPredictedMolsByTop3OfExperiments.push_back(0.0);
						}
						dvecPeptidesIntensitiesTemp.clear();
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_mapExperimentsAndPeptidesIntensity.\
							insert(pair<string, vector<double>>(m_vecStrExperiments[iBegin+1], dvecPeptidesIntensitiesTemp));
					}

					MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecExperiments.push_back(m_vecStrExperiments.at(i));
					MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecLFAQpepOfExperiments.push_back(m_vecProteins[j].m_dProteinLFAQpep);
					MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecPredictedMolsByLFAQpepOfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfLFAQpep);
					if (m_Params.m_bIfCalculateiBAQ)
					{
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecMaxQuantiBAQOfExperiments.push_back(m_vecProteins[j].m_dMaxQuantiBAQ);
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecPredictedMolsByiBAQOfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfiBAQ);
					}
					if (m_Params.m_bIfCalculateTop3)
					{
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecTop3OfExperiments.push_back(m_vecProteins[j].m_dProteinIntensityTop3);
						MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_vecPredictedMolsByTop3OfExperiments.\
							push_back(m_vecProteins[j].m_dPredictedMolOfTop3);
					}
					dvecPeptidesIntensitiesTemp.clear();
					for (int k = 0; k < m_vecProteins[j].m_vPeptidesCorrectedIntensity.size(); k++)
					{
						dvecPeptidesIntensitiesTemp.push_back(m_vecProteins[j].m_vPeptidesCorrectedIntensity[k]);
					}
					MergedProteins.m_vecMergedProteins[mapProteinNamesAndIndexIter->second].m_mapExperimentsAndPeptidesIntensity.\
						insert(pair<string, vector<double>>(m_vecStrExperiments[i], dvecPeptidesIntensitiesTemp));
				}
				else
				{ // new protein
					mergedProteinTemp.Clear();
					mergedProteinTemp.m_strPrtoteinName = m_vecProteins[j].m_strProteinID;
					mergedProteinTemp.m_strSubsetProteins = m_vecProteins[j].m_strProteinFullName;
					dvecPeptidesIntensitiesTemp.clear();
					//make up missing experiment
					for (int t = 0; t < i; t++)
					{
						mergedProteinTemp.m_vecLFAQpepOfExperiments.push_back(0.0);
						mergedProteinTemp.m_vecPredictedMolsByLFAQpepOfExperiments.push_back(0.0);
						if (m_Params.m_bIfCalculateiBAQ)
						{
							mergedProteinTemp.m_vecMaxQuantiBAQOfExperiments.push_back(0.0);
							mergedProteinTemp.m_vecPredictedMolsByiBAQOfExperiments.push_back(0.0);
						}
						if (m_Params.m_bIfCalculateTop3)
						{
							mergedProteinTemp.m_vecTop3OfExperiments.push_back(0.0);
							mergedProteinTemp.m_vecPredictedMolsByTop3OfExperiments.push_back(0.0);
						}
						mergedProteinTemp.m_mapExperimentsAndPeptidesIntensity.insert(pair < string, vector<double>>(m_vecStrExperiments[t], dvecPeptidesIntensitiesTemp));
					}
					mergedProteinTemp.m_vecExperiments.push_back(m_vecStrExperiments.at(i));
					mergedProteinTemp.m_vecLFAQpepOfExperiments.push_back(m_vecProteins[j].m_dProteinLFAQpep);
					mergedProteinTemp.m_vecPredictedMolsByLFAQpepOfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfLFAQpep);
					if (m_Params.m_bIfCalculateiBAQ)
					{
						mergedProteinTemp.m_vecMaxQuantiBAQOfExperiments.push_back(m_vecProteins[j].m_dMaxQuantiBAQ);
						mergedProteinTemp.m_vecPredictedMolsByiBAQOfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfiBAQ);
					}
					if (m_Params.m_bIfCalculateTop3)
					{
						mergedProteinTemp.m_vecTop3OfExperiments.push_back(m_vecProteins[j].m_dProteinIntensityTop3);
						mergedProteinTemp.m_vecPredictedMolsByTop3OfExperiments.push_back(m_vecProteins[j].m_dPredictedMolOfTop3);
					}
					dvecPeptidesIntensitiesTemp.clear();
					for (int k = 0; k < m_vecProteins[j].m_vPeptidesCorrectedIntensity.size(); k++)
					{
						dvecPeptidesIntensitiesTemp.push_back(m_vecProteins[j].m_vPeptidesCorrectedIntensity[k]);
					}
					mergedProteinTemp.m_mapExperimentsAndPeptidesIntensity.insert(\
						pair<string, vector<double>>(m_vecStrExperiments[i], dvecPeptidesIntensitiesTemp));

					MergedProteins.m_vecMergedProteins.push_back(mergedProteinTemp);
					mapProteinNamesAndIndex.insert(pair<string, int>(m_vecProteins[j].m_strProteinID, MergedProteins.m_vecMergedProteins.size() - 1));
				}
			} // end for m_vecProteins

		} //end not first experiment
	} // end for m_vecStrExperiments

	vector<CMergedProtein>::iterator mergeProteinIter;
	mergeProteinIter = MergedProteins.m_vecMergedProteins.begin();
	int iExperimentNum = m_vecStrExperiments.size();
	for (; mergeProteinIter != MergedProteins.m_vecMergedProteins.end(); mergeProteinIter++)
	{
		int iBegin = 0;
		for (; mergeProteinIter->m_vecExperiments.back() != m_vecStrExperiments[iBegin]; iBegin++){};
		if (mergeProteinIter->m_vecExperiments.back() != m_vecStrExperiments.back())
		{
			for (; iBegin < iExperimentNum - 1; iBegin++)
			{
				mergeProteinIter->m_vecLFAQpepOfExperiments.push_back(0.0);
				mergeProteinIter->m_vecPredictedMolsByLFAQpepOfExperiments.push_back(0.0);
				if (m_Params.m_bIfCalculateiBAQ)
				{
					mergeProteinIter->m_vecMaxQuantiBAQOfExperiments.push_back(0.0);
					mergeProteinIter->m_vecPredictedMolsByiBAQOfExperiments.push_back(0.0);
				}
				if (m_Params.m_bIfCalculateTop3)
				{
					mergeProteinIter->m_vecTop3OfExperiments.push_back(0.0);
					mergeProteinIter->m_vecPredictedMolsByTop3OfExperiments.push_back(0.0);
				}
				dvecPeptidesIntensitiesTemp.clear();
				mergeProteinIter->m_mapExperimentsAndPeptidesIntensity.\
					insert(pair<string, vector<double>>(m_vecStrExperiments[iBegin+1], dvecPeptidesIntensitiesTemp));
			}
		}

		assert(m_vecStrExperiments.size() == mergeProteinIter->m_vecLFAQpepOfExperiments.size());
		assert(m_vecStrExperiments.size() == mergeProteinIter->m_vecPredictedMolsByLFAQpepOfExperiments.size());
		if (m_Params.m_bIfCalculateiBAQ)
		{
			assert(m_vecStrExperiments.size() == mergeProteinIter->m_vecMaxQuantiBAQOfExperiments.size());
			assert(m_vecStrExperiments.size() == mergeProteinIter->m_vecPredictedMolsByiBAQOfExperiments.size());
		}
		if (m_Params.m_bIfCalculateTop3)
		{
			assert(m_vecStrExperiments.size() == mergeProteinIter->m_vecTop3OfExperiments.size());
			assert(m_vecStrExperiments.size() == mergeProteinIter->m_vecPredictedMolsByTop3OfExperiments.size());
		}

	}

	MergedProteins.mf_ExperimentsAnalysis();

	mf_SaveMergedProtein(MergedProteinResultPath, MergedProteins);
}

void CProteinWorker::mf_SaveMergedProtein(const string& MergedProteinResultPath, const CMergedProteins& MergedProteins)
{
	ofstream oMergefile;
	oMergefile.open(MergedProteinResultPath.c_str());
	if (!oMergefile)
	{
		cout << "Error:\tCannot open " << MergedProteinResultPath << endl;
		flog.mf_Input("Error:\tCannot open " + MergedProteinResultPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "\tMerged " << MergedProteins.m_vecMergedProteins.size() << " proteins' LFAQ into " << MergedProteinResultPath << endl;
		flog.mf_Input("\tMerged " + fSize_t2String(MergedProteins.m_vecMergedProteins.size()) + " proteins' LFAQ into " + MergedProteinResultPath + "\n");
	}

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	oMergefile << "Protein IDs\tMajority protein IDs\tExperiments\t";
	oMergefile << "LFAQ\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		oMergefile << "iBAQ\t";
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		oMergefile << "Top3\t";
	}
	oMergefile << "PredictedMol(LFAQ)\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		oMergefile << "PredictedMol(iBAQ)\t";
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		oMergefile << "PredictedMol(Top3)\t";
	}
	oMergefile << "LFAQCV\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		oMergefile << "iBAQCV\t";
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		oMergefile << "Top3CV\t";
	}
	oMergefile << "PeptidesIntensities\n";

	int indexTemp = 0;
	int j = 0;
	map<string, vector<double>>::const_iterator mapExperimentAndPeptidesIntensitiesIter;
	vector<double>::const_iterator vecPeptidesintensitiesIter;
	int iCountTemp;
	for (indexTemp = 0; indexTemp < MergedProteins.m_vecMergedProteins.size(); indexTemp++)
	{
		oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_strSubsetProteins << "\t";
		oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_strPrtoteinName << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecExperiments.size() - 1; j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecExperiments.at(j) << ";";
		}
		oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecExperiments.at(j) << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecLFAQpepOfExperiments.size() - 1; j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecLFAQpepOfExperiments.at(j) << ";";
		}
		oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecLFAQpepOfExperiments.at(j) << "\t";
		if (m_Params.m_bIfCalculateiBAQ)
		{
			for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecMaxQuantiBAQOfExperiments.size() - 1; j++)
			{
				oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecMaxQuantiBAQOfExperiments.at(j) << ";";
			}
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecMaxQuantiBAQOfExperiments.at(j) << "\t";
		}
		if (m_Params.m_bIfCalculateTop3)
		{
			for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecTop3OfExperiments.size() - 1; j++)
			{
				oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecTop3OfExperiments.at(j) << ";";
			}
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecTop3OfExperiments.at(j) << "\t";
		}
		for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByLFAQpepOfExperiments.size() - 1; j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByLFAQpepOfExperiments.at(j) << ";";
		}
		oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByLFAQpepOfExperiments.at(j) << "\t";
		if (m_Params.m_bIfCalculateiBAQ)
		{
			for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByiBAQOfExperiments.size() - 1; j++)
			{
				oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByiBAQOfExperiments.at(j) << ";";
			}
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByiBAQOfExperiments.at(j) << "\t";
		}
		if (m_Params.m_bIfCalculateTop3)
		{
			for (j = 0; j < MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByTop3OfExperiments.size() - 1; j++)
			{
				oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByTop3OfExperiments.at(j) << ";";
			}
			oMergefile << MergedProteins.m_vecMergedProteins[indexTemp].m_vecPredictedMolsByTop3OfExperiments.at(j) << "\t";
		}

		oMergefile << MergedProteins.m_vecLFAQpepCV.at(indexTemp) << "\t";
		if (m_Params.m_bIfCalculateiBAQ)
		{
			oMergefile << MergedProteins.m_vecMaxQuantiBAQCV.at(indexTemp) << "\t";
		}
		if (m_Params.m_bIfCalculateTop3)
		{
			oMergefile << MergedProteins.m_vecTop3CV.at(indexTemp) << "\t";
		}
		mapExperimentAndPeptidesIntensitiesIter = MergedProteins.m_vecMergedProteins[indexTemp].m_mapExperimentsAndPeptidesIntensity.begin();
		for (iCountTemp = 0; mapExperimentAndPeptidesIntensitiesIter != MergedProteins.m_vecMergedProteins[indexTemp].m_mapExperimentsAndPeptidesIntensity.end(); \
			mapExperimentAndPeptidesIntensitiesIter++, iCountTemp++)
		{
			vecPeptidesintensitiesIter = mapExperimentAndPeptidesIntensitiesIter->second.begin();

			if (mapExperimentAndPeptidesIntensitiesIter->second.size() > 0)
			{
				for (; vecPeptidesintensitiesIter != mapExperimentAndPeptidesIntensitiesIter->second.end() - 1; vecPeptidesintensitiesIter++)
				{
					oMergefile << *vecPeptidesintensitiesIter << ",";
				}
				if (iCountTemp != MergedProteins.m_vecMergedProteins[indexTemp].m_mapExperimentsAndPeptidesIntensity.size() - 1)
				{
					oMergefile << *vecPeptidesintensitiesIter << ";";
				}
				else
				{
					oMergefile << *vecPeptidesintensitiesIter;
				}
			}
			else
			{
				if (iCountTemp != MergedProteins.m_vecMergedProteins[indexTemp].m_mapExperimentsAndPeptidesIntensity.size() - 1)
				{
					oMergefile << 0.0 << ";";
				}
				else
				{
					oMergefile << 0.0;
				}
			}
		}
		oMergefile << endl;
	}
	flog.mf_Input("Done\n");
}
void CProteinWorker::mf_MergeRegressionResult( string RegressionResultPath)
{
	string strRegressionResultOfOneExperimentPath;
	FILE * Ofile; 
	char Buffer[BUFFERLENGTH];
	char *pstr;
	string strIDTemp;
	string strTemp;

	ofstream oMergefile;
	oMergefile.open(RegressionResultPath.c_str());
	if (!oMergefile)
	{
		cout << "Error:\tCannot open " << RegressionResultPath << endl;
		flog.mf_Input("Error:\tCannot open " + RegressionResultPath+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		flog.mf_Input("Merging regression results from different experiments into " + RegressionResultPath);
	}
	for (int i = 0; i < m_vecStrExperiments.size(); i++)
	{
		strRegressionResultOfOneExperimentPath = m_Params.m_strRegressionResult + m_vecStrExperiments.at(i) + ".txt";
		errno_t err;
		err = fopen_s(&Ofile, strRegressionResultOfOneExperimentPath.c_str(), "r");
		if (err != 0)
		{
			flog.mf_Input("Error:\tCannot open " + strRegressionResultOfOneExperimentPath+"\n");
			flog.mf_Destroy();
			exit(1);
		}

		fgets(Buffer, BUFFERLENGTH, Ofile);
		pstr = Buffer;
		oMergefile << m_vecStrExperiments.at(i) << endl;
		while (!feof(Ofile))
		{
			oMergefile << pstr;
			fgets(Buffer, BUFFERLENGTH, Ofile);
			pstr = Buffer;
		}

		fclose(Ofile);
	}
	flog.mf_Input("...Done\n");

}
//to validation the correctness of the program, calculate the pearson correlation 
// between  the UPS2 protein's quantifications and UPS2 protein's real moles.
void CProteinWorker::mf_ShowUPS2AnalysisResult( int iExperimentIndex)
{
	map<string, double> mapUPS2Id2Mols;
	map<string, double> mapUPS2iBAQs;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	map<string, double>::iterator UPS2iBAQIter;
	map<string, int>::iterator mapLowSpikedInUPS2Iter;
	vector<CProtein>::iterator proteinIter;
	string strTemp;

	mf_LoadUPS2Mols(m_Params.m_strStandardProteinsFilePath, mapUPS2Id2Mols);

	double dPearsonCorr = 0.0;
	vector<string> vecUPS2IDIdentified;
	vector<double> vecUPS2logLFAQpep;
	vector<double> vecUPS2logLFAQpro;
	vector<double> vecUPS2logMolsIdentified;
	vector<double> vecUPS2logiBAQIdentified;
	vector<double> vecUPS2logTop3;

	string::size_type  stPosition;
	ofstream ofile;
	int iPosTemp = m_Params.m_strProteinLFAQpepPath.find_last_of("\\");
	string CorrelationPath = m_Params.m_strProteinLFAQpepPath.substr(0, iPosTemp + 1) + "LFAQResultsForStandardProteins" + m_vecStrExperiments.at(iExperimentIndex) + ".txt";

	ofile.open(CorrelationPath);
	if (!ofile)
	{
		cout << "Error:\tCannot open " << CorrelationPath << endl;
		flog.mf_Input("Error:\tCannot open " + CorrelationPath+"\n");
		flog.mf_Destroy();
		exit(-1);
	}
	// Calculate the sequence coverage 
	vector<int> vecAppearedSequence;
	int iBegin, iEnd;
	vector<double> vecdCoverageOfUPS2;
	double dCoverageTemp;
	vector<string> vecProteinsTemp;
	bool bIfOthersOnlyContam;
	ofile << "Protein ID\tNumberOfUniquePeptides\tSpikedInMols\tLFAQ\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		ofile << "iBAQ\t";
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		ofile << "Top3\t";
	}
	ofile << "PredictedMol(LFAQ)\t";
	if (m_Params.m_bIfCalculateiBAQ)
	{
		ofile << "PredictedMol(iBAQ)\t";
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		ofile << "PredictedMol(Top3)\t";
	}
	ofile << "Coverage\tProteinSequenceLength\t\n";
	
    ofile << setprecision(6);
    ofile.setf(ios::fixed);

	double dPeptideIntSum;
	double dMaxintensityTemp;
	int iMaxIndex;
	double IdentifiedQfactorMax;
	double TheoryQfactorMax;
	double SequenceQfactorMax;
	double dlogIdentifiedQfactorMax;
	double dlogTheoryQfactorMax;
	double dlogSequenceQfactorMax;
	size_t i,j;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateLFAQpep == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stPosition = strTemp.find(m_Params.m_strIdentifierOfStandPro);
			vecProteinsTemp = split(proteinIter->m_strProteinFullName, ";");
			if (stPosition != strTemp.npos)
			{
				bIfOthersOnlyContam = true;
				if (m_Params.m_bIfExistContanProtein == true && vecProteinsTemp.size() > 1)
				{
					for (i = 1; i < vecProteinsTemp.size(); i++)
					{
						if (vecProteinsTemp[i].find(m_Params.m_strContaminantfix) == vecProteinsTemp[i].npos)
						{
							bIfOthersOnlyContam = false;
						}
					}
				}
				if (proteinIter->m_strProteinFullName.find(";") == proteinIter->m_strProteinFullName.npos || bIfOthersOnlyContam==true)
				{ //
					Ups2MolsIter = mapUPS2Id2Mols.begin();
					for (; Ups2MolsIter != mapUPS2Id2Mols.end(); Ups2MolsIter++)
					{
						if (strTemp.find(Ups2MolsIter->first) != strTemp.npos)
						{
							dPeptideIntSum = 0.0;
							for (i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
							{
								dPeptideIntSum += proteinIter->m_vPeptidesNativeIntensity[i];
							}
							dMaxintensityTemp = GetMaxValue(proteinIter->m_vPeptidesNativeIntensity);
							iMaxIndex = GetMaxIndex(proteinIter->m_vPeptidesNativeIntensity);
							IdentifiedQfactorMax = (proteinIter->m_vPeptidesNativeIntensity[iMaxIndex] / dPeptideIntSum)*proteinIter->m_iPeptidesNumber;
							TheoryQfactorMax = (proteinIter->m_vPeptidesNativeIntensity[iMaxIndex] / dPeptideIntSum)*proteinIter->m_iNumberOfTheoreticEnzyme;
							SequenceQfactorMax = (proteinIter->m_vPeptidesNativeIntensity[iMaxIndex] / dPeptideIntSum)*proteinIter->m_strProteinSequence.size();
							dlogIdentifiedQfactorMax = log(proteinIter->m_vPeptidesNativeIntensity[iMaxIndex]) / log(dPeptideIntSum/proteinIter->m_iPeptidesNumber);
							dlogTheoryQfactorMax = log(proteinIter->m_vPeptidesNativeIntensity[iMaxIndex]) / log(dPeptideIntSum / proteinIter->m_iNumberOfTheoreticEnzyme);
							dlogSequenceQfactorMax = log(proteinIter->m_vPeptidesNativeIntensity[iMaxIndex]) / log(dPeptideIntSum / proteinIter->m_strProteinSequence.size());
							// Calculate the sequence coverage 
							vecAppearedSequence.clear();
							for (i = 0; i < proteinIter->m_strProteinSequence.length(); i++)
								vecAppearedSequence.push_back(0);
							for (i = 0; i < proteinIter->m_vPeptidesSequences.size(); i++)
							{
								iBegin = proteinIter->m_strProteinSequence.find(proteinIter->m_vPeptidesSequences.at(i));
								if (iBegin == string::npos)
								{
									cout << "Cannot find " << proteinIter->m_vPeptidesSequences.at(i) << " in " << proteinIter->m_strProteinID << "'s sequence\n";
									break;
								}
								iEnd = iBegin + proteinIter->m_vPeptidesSequences.at(i).size() - 1;
								for (j = iBegin; j <= iEnd; j++)
								{
									vecAppearedSequence.at(j) = 1;
								}
							}
							dCoverageTemp = 0.0;
							for (i = 0; i < vecAppearedSequence.size(); i++)
							{
								dCoverageTemp += vecAppearedSequence.at(i);
							}
							dCoverageTemp = dCoverageTemp / proteinIter->m_strProteinSequence.size();

							vecUPS2IDIdentified.push_back(Ups2MolsIter->first);
							ofile << Ups2MolsIter->first << "\t";
							ofile << proteinIter->m_iPeptidesNumber << "\t";
							ofile << Ups2MolsIter->second << "\t";

							if (proteinIter->m_dProteinLFAQpep != 0.0 && (proteinIter->m_dMaxQuantiBAQ != 0.0 || (!m_bIfiBAQIntensityExist)))
							{
								if (Ups2MolsIter->second == 0.0)
								{
									vecUPS2logMolsIdentified.push_back(0.0);
								}
								else
								{
									vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
								}
								if (proteinIter->m_dProteinLFAQpep == 0.0)
								{
									vecUPS2logLFAQpep.push_back(0.0);
								}
								else
								{
									vecUPS2logLFAQpep.push_back(log10(proteinIter->m_dProteinLFAQpep));
								}
								if (proteinIter->m_dProteinIntensityTop3 == 0.0)
								{
									vecUPS2logTop3.push_back(0.0);
								}
								else
								{
									vecUPS2logTop3.push_back(log10(proteinIter->m_dProteinIntensityTop3));
								}
								
								ofile << proteinIter->m_dProteinLFAQpep << "\t";

								if (m_Params.m_bIfCalculateiBAQ)
								{
									if (m_bIfiBAQIntensityExist)
									{
										if (proteinIter->m_dMaxQuantiBAQ == 0.0)
										{
											vecUPS2logiBAQIdentified.push_back(0.0);
										}
										else
										{
											vecUPS2logiBAQIdentified.push_back(log10(proteinIter->m_dMaxQuantiBAQ));
										}
										
										ofile << proteinIter->m_dMaxQuantiBAQ << "\t";
									}
									else
									{
										if (proteinIter->m_dReCalculateiBAQ == 0.0)
										{
											vecUPS2logiBAQIdentified.push_back(0.0);
										}
										else
										{
											vecUPS2logiBAQIdentified.push_back(log10(proteinIter->m_dReCalculateiBAQ));
										}
										
										ofile << proteinIter->m_dReCalculateiBAQ << "\t";
									}
								}
								if (m_Params.m_bIfCalculateTop3)
								{
									ofile << proteinIter->m_dProteinIntensityTop3 << "\t";
								}

								ofile << proteinIter->m_dPredictedMolOfLFAQpep << "\t";
								if (m_Params.m_bIfCalculateiBAQ)
								{
									ofile << proteinIter->m_dPredictedMolOfiBAQ << "\t";
								}
								if (m_Params.m_bIfCalculateTop3)
								{
									ofile << proteinIter->m_dPredictedMolOfTop3 << "\t";
								}
								ofile << dCoverageTemp << "\t";
								
							} // end if 0.0
							else
							{
								ofile << proteinIter->m_dProteinLFAQpep<< "\t";
								if (m_Params.m_bIfCalculateTop3)
								{
									ofile << proteinIter->m_dProteinIntensityTop3 << "\t";
								}
								if (m_Params.m_bIfCalculateiBAQ)
								{
									ofile << proteinIter->m_dMaxQuantiBAQ << "\t";
								}
								ofile << proteinIter->m_dPredictedMolOfLFAQpep << "\t";
								if (m_Params.m_bIfCalculateTop3)
								{
									ofile << proteinIter->m_dPredictedMolOfTop3 << "\t";
								}
								if (m_Params.m_bIfCalculateiBAQ)
								{
									ofile << proteinIter->m_dPredictedMolOfiBAQ << "\t";
								}
								ofile << dCoverageTemp << "\t";
							}
							ofile << proteinIter->m_strProteinSequence.length() << "\t";
							ofile << endl;
							break;
						}
					} // for mapUPS2Id2Mols
					if (Ups2MolsIter == mapUPS2Id2Mols.end())
					{
						mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
						if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
						{
							cout << "Warning:\tCannot find protein " << strTemp << " in StandardProteins.txt\n";
							flog.mf_Input("Warning:\tCannot find protein " + strTemp + " in StandardProteins.txt\n");

						}
					}
				}
			}
						
		}
	
	}
	ofile.close();
	
	if (vecUPS2IDIdentified.size() == 0)
	{
		cout << "Error:\tCannot find any standard protein in the input file using the standard protein identifier " << m_Params.m_strIdentifierOfStandPro << endl;
		flog.mf_Input("Error:\tCannot find any standard protein in the input file using the standard protein identifier " +m_Params.m_strIdentifierOfStandPro+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	
	cout << "\t\tR2 between predicted and real amounts of " << vecUPS2logMolsIdentified.size() << "  standard proteins:" << endl;
	flog.mf_Input("\t\tR2 between predicted and real amounts of " + fSize_t2String(vecUPS2logMolsIdentified.size()) + "  standard proteins:\n");
	cout << "\t\t\tLFAQ\t";
	flog.mf_Input("\t\t\tLFAQ\t");
	if (m_Params.m_bIfCalculateiBAQ)
	{
		cout << "iBAQ\t";
		flog.mf_Input("iBAQ\t");
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		cout << "Top3\n";
		flog.mf_Input("Top3\n");
	}
	if (!m_Params.m_bIfCalculateiBAQ&&!m_Params.m_bIfCalculateTop3)
	{
		cout << "\n";
		flog.mf_Input("\n");
	}
	dPearsonCorr = spearsonCorrelation(vecUPS2logLFAQpep, vecUPS2logMolsIdentified);
	cout<<"\t\t\t"<<round(dPearsonCorr*dPearsonCorr,3) <<"\t";
	flog.mf_Input("\t\t\t" + fDouble2String(round(dPearsonCorr*dPearsonCorr,3)) + "\t");
	
	if (m_Params.m_bIfCalculateiBAQ)
	{
		dPearsonCorr = spearsonCorrelation(vecUPS2logiBAQIdentified, vecUPS2logMolsIdentified);
		cout << round(dPearsonCorr*dPearsonCorr,3) << "\t";
		flog.mf_Input(fDouble2String(round(dPearsonCorr*dPearsonCorr,3)) + "\t");
	}
	if (m_Params.m_bIfCalculateTop3)
	{
		dPearsonCorr = spearsonCorrelation(vecUPS2logTop3, vecUPS2logMolsIdentified);
		cout << round(dPearsonCorr*dPearsonCorr,3) << "\n";
		flog.mf_Input(fDouble2String(round(dPearsonCorr*dPearsonCorr,3)) + "\n");
	}
	if (!m_Params.m_bIfCalculateiBAQ&&!m_Params.m_bIfCalculateTop3)
	{
		cout << "\n";
		flog.mf_Input("\n");
	}
}


void CProteinWorker::mf_PredictMolsByStandProteinsAndLFAQ(int iExperimentIndex)
{

	map<string, double> mapUPS2Id2Mols;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	vector<CProtein>::iterator proteinIter;
	string strTemp;

	cout << "Predict proteins' abundance by the standard proteins' LFAQ\n";
	flog.mf_Input("Predict proteins' abundance by the standard proteins' LFAQ\n");
	mf_LoadUPS2Mols(m_Params.m_strStandardProteinsFilePath, mapUPS2Id2Mols);

	double dPearsonCorr = 0.0;
	vector<double> veclogLFAQpepOfStandProteins;//the LFAQpep of identified standard proteins
	vector<double> veclogLFAQproOfStandProteins; //the LFAQpro of identified standard proteins
	vector<double> vecUPS2logMolsIdentified; // the spiked-in amount of identified standard proteins;
	vector<double> vecPredictedLogMolsOfStandProteinsByLFAQpep;
	vector<double> vecPredictedLogMolsOfStandProteinsByLFAQpro;

	int iNumberofPeptidesOfStandProteins = 0;
	string::size_type stStart, stPosition, stEnd;
	vector<string> vecProteinsTemp;
	bool bIfOthersOnlyContam = true;
	map<string, int>::iterator mapLowSpikedInUPS2Iter;

	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateLFAQpep == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stPosition = strTemp.find(m_Params.m_strIdentifierOfStandPro);
			vecProteinsTemp = split(proteinIter->m_strProteinFullName, ";");
			if (stPosition != strTemp.npos)
			{
				bIfOthersOnlyContam = true;
				if (m_Params.m_bIfExistContanProtein == true && vecProteinsTemp.size() > 1)
				{
					for (size_t i = 1; i < vecProteinsTemp.size(); i++)
					{
						if (vecProteinsTemp[i].find(m_Params.m_strContaminantfix) == vecProteinsTemp[i].npos)
						{
							bIfOthersOnlyContam = false;
						}
					}
				}
				if ((proteinIter->m_strProteinFullName.find(";") == \
					proteinIter->m_strProteinFullName.npos) || (bIfOthersOnlyContam == true))
				{ //
					Ups2MolsIter = mapUPS2Id2Mols.begin();
					for (; Ups2MolsIter != mapUPS2Id2Mols.end(); Ups2MolsIter++)
					{
						if (strTemp.find(Ups2MolsIter->first) != strTemp.npos)
						{
							if (proteinIter->m_dProteinLFAQpep != 0.0)
							{
								iNumberofPeptidesOfStandProteins += proteinIter->m_vPeptidesSequences.size();
								if (Ups2MolsIter->second == 0.0)
								{
									vecUPS2logMolsIdentified.push_back(0.0);
									cout << "Warning:\tThe mols of standard protein " << Ups2MolsIter->first << " is 0\n";
									flog.mf_Input("Warning:\tThe mols of standard protein " + Ups2MolsIter->first + " is 0\n");
								}
								else
								{
									vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
								}

								veclogLFAQpepOfStandProteins.push_back(log10(proteinIter->m_dProteinLFAQpep));

								if (proteinIter->m_dProteinLFAQpro == 0.0)
								{
									veclogLFAQproOfStandProteins.push_back(0.0);
								}
								else
								{
									veclogLFAQproOfStandProteins.push_back(log10(proteinIter->m_dProteinLFAQpro));
								}

							} // end if 0.0
							break;
						} // end if strTemp
					} // for mapUPS2Id2Mols

					if (Ups2MolsIter == mapUPS2Id2Mols.end())
					{
						mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
						if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
						{
							cout << "Warning:\tCannot find protein " << strTemp << " in StandardProteins.txt\n";
							flog.mf_Input("Warning:\tCannot find protein " + strTemp + " in StandardProteins.txt\n");
						}
					}

				}
			} // end if (stPosition != strTemp.npos)

		}// end if m_bIfCalculateLFAQpep
	}// end for m_vecProteins

	cout << "\tLoaded " << vecUPS2logMolsIdentified.size() << " standard proteins " << " with " << \
		iNumberofPeptidesOfStandProteins << " peptides\n";
	flog.mf_Input("\tLoaded " + fSize_t2String(vecUPS2logMolsIdentified.size()) + " standard proteins with " + \
		fInt2String(iNumberofPeptidesOfStandProteins) + " peptides\n");


	BootLineaRegress bootLinear;
	double a, b;
	double DevSQ;
	double sd;
	double RSS;
	double Maxd;
	double Mind;
	double Meand;
	double SdOfa;
	double SdOfb;

	cout << "\tPredict protein amount by bootstrap linear regression\n";
	flog.mf_Input("\tPredict protein amount by bootstrap linear regression\n");

	double R2;

	// Predict protein amount by LFAQpep
	cout << "\t\tPredict protein amount by LFAQ" << endl;
	flog.mf_Input("\t\tPredict protein amount by LFAQ\n");

	// Normalized
	//vector<double> vecNormalizedX;
	//vector<double> vecNormalizedY;
	//bootLinear.mf_Normalization(veclogLFAQpepOfStandProteins, vecNormalizedX, sdx); //or veclogLFAQproOfStandProteins
	//bootLinear.mf_Normalization(vecUPS2logMolsIdentified, vecNormalizedY,sdy);
	//bootLinear.mf_BootstrapLinearRegression(vecNormalizedX, vecNormalizedY, 10000, a, b,SdOfa,SdOfb, DevSQ, sd, RSS, Maxd, Mind, Meand);
	//a = (sdy / sdx)*a;
	//b = (sdy / sdx)*b;
	//cout << "sdx= " << sdx << endl;
	//cout << "sdy= " << sdy <<endl;

	// No normalized
	bootLinear.mf_BootstrapLinearRegression(veclogLFAQpepOfStandProteins, vecUPS2logMolsIdentified, 10000, a, b, SdOfa, SdOfb, DevSQ, sd, RSS, Maxd, Mind, Meand);


	bootLinear.mf_Predict(veclogLFAQpepOfStandProteins, a, b, vecPredictedLogMolsOfStandProteinsByLFAQpep);
	R2 = spearsonCorrelation(vecUPS2logMolsIdentified, vecPredictedLogMolsOfStandProteinsByLFAQpep);
	cout << "\t\t\ta= " << a << endl;
	cout << "\t\t\tb= " << b << endl;
	cout << "\t\t\tSd of a= " << SdOfa << endl;
	cout << "\t\t\tSd of b= " << SdOfb << endl;
	cout << "\t\t\tMean Standard Error= " << sd << endl;
	cout << "\t\t\tAutocorrelation R2= " << R2*R2 << endl;
	flog.mf_Input("\t\t\ta= " + fDouble2String(a) + "\n");
	flog.mf_Input("\t\t\tb= " + fDouble2String(b) + "\n");
	flog.mf_Input("\t\t\tSd of a= " + fDouble2String(SdOfa) + "\n");
	flog.mf_Input("\t\t\tSd of b= " + fDouble2String(SdOfb) + "\n");
	flog.mf_Input("\t\t\tMean Standard Error= " + fDouble2String(sd) + "\n");
	flog.mf_Input("\t\t\tAutocorrelation R2= " + fDouble2String(R2*R2) + "\n");
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateLFAQpep == true)
		{
			if (proteinIter->m_dProteinLFAQpep == 0.0)
			{
				proteinIter->m_dPredictedMolOfLFAQpep = 0.0;
			}
			else
			{
				proteinIter->m_dPredictedMolOfLFAQpep = pow(10, bootLinear.mf_Predict(log10(proteinIter->m_dProteinLFAQpep), a, b));
			}
		}
	}

}
void CProteinWorker::mf_PredictMolsByStandProteinsAndiBAQ(int iExperimentIndex)
{
	map<string, double> mapUPS2Id2Mols;
	map<string, double> mapUPS2iBAQs;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	map<string, double>::iterator UPS2iBAQIter;
	vector<CProtein>::iterator proteinIter;
	string strTemp;

	cout << "Predict proteins' abundance by the standard proteins' iBAQ\n";
	flog.mf_Input("Predict proteins' abundance by the standard proteins'iBAQ\n");
	mf_LoadUPS2Mols(m_Params.m_strStandardProteinsFilePath, mapUPS2Id2Mols);


	double dPearsonCorr = 0.0;
	vector<double> veclogiBAQOfStandProteins; //the iBAQ of identified standard proteins
	vector<double> vecUPS2logMolsIdentified; // the spiked-in amount of identified standard proteins;
	vector<double> vecPredictedLogMolsOfStandProteinsByiBAQ;

	int iNumberofPeptidesOfStandProteins = 0;
	string::size_type stStart, stPosition, stEnd;
	vector<string> vecProteinsTemp;
	bool bIfOthersOnlyContam = true;
	map<string, int>::iterator mapLowSpikedInUPS2Iter;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		strTemp = proteinIter->m_strProteinID;
		stPosition = strTemp.find(m_Params.m_strIdentifierOfStandPro);
		vecProteinsTemp = split(proteinIter->m_strProteinFullName, ";");
		if (stPosition != strTemp.npos)
		{
			bIfOthersOnlyContam = true;
			if (m_Params.m_bIfExistContanProtein == true && vecProteinsTemp.size() > 1)
			{
				for (size_t i = 1; i < vecProteinsTemp.size(); i++)
				{
					if (vecProteinsTemp[i].find(m_Params.m_strContaminantfix) == vecProteinsTemp[i].npos)
					{
						bIfOthersOnlyContam = false;
					}
				}
			}
			if ((proteinIter->m_strProteinFullName.find(";") == \
				proteinIter->m_strProteinFullName.npos) || (bIfOthersOnlyContam == true))
			{ //
				Ups2MolsIter = mapUPS2Id2Mols.begin();
				for (; Ups2MolsIter != mapUPS2Id2Mols.end(); Ups2MolsIter++)
				{
					if (strTemp.find(Ups2MolsIter->first) != strTemp.npos)
					{
						if (proteinIter->m_dMaxQuantiBAQ != 0.0 || (!m_bIfiBAQIntensityExist))
						{
							iNumberofPeptidesOfStandProteins += proteinIter->m_vPeptidesSequences.size();
							if (Ups2MolsIter->second == 0.0)
							{
								vecUPS2logMolsIdentified.push_back(0.0);
							}
							else
							{
								vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
							}


							if (m_bIfiBAQIntensityExist)
							{
								if (proteinIter->m_dMaxQuantiBAQ == 0.0)
								{
									veclogiBAQOfStandProteins.push_back(0.0);
								}
								else
								{
									veclogiBAQOfStandProteins.push_back(log10(proteinIter->m_dMaxQuantiBAQ));
								}

							}
							else
							{
								if (proteinIter->m_dReCalculateiBAQ == 0.0)
								{
									veclogiBAQOfStandProteins.push_back(0.0);
								}
								else
								{
									veclogiBAQOfStandProteins.push_back(log10(proteinIter->m_dReCalculateiBAQ));
								}
							}

						} // end if 0.0
						break;
					}
				} // for mapUPS2Id2Mols

				if (Ups2MolsIter == mapUPS2Id2Mols.end())
				{
					mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
					if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
					{
						cout << "Warning:\tCannot find protein " << strTemp << " in StandardProteins.txt\n";
						flog.mf_Input("Warning:\tCannot find protein " + strTemp + " in StandardProteins.txt\n");
					}
				}

			}
		} // end if (stPosition != strTemp.npos)

	}// end for m_vecProteins

	cout << "\tLoaded " << vecUPS2logMolsIdentified.size() << " standard proteins " << " with " << \
		iNumberofPeptidesOfStandProteins << " peptides\n";
	flog.mf_Input("\tLoaded " + fSize_t2String(vecUPS2logMolsIdentified.size()) + " standard proteins with " + \
		fInt2String(iNumberofPeptidesOfStandProteins) + " peptides\n");

	BootLineaRegress bootLinear;
	double a, b;
	double DevSQ;
	double sd;
	double RSS;
	double Maxd;
	double Mind;
	double Meand;
	double SdOfa;
	double SdOfb;
	double R2;

	cout << "\tPredict protein amount by bootstrap linear regression\n";
	flog.mf_Input("\tPredict protein amount by bootstrap linear regression\n");


	// Predict protein amount by iBAQ
	if (m_Params.m_bIfCalculateiBAQ)
	{
		cout << "\t\tPredict protein amount by iBAQ" << endl;
		flog.mf_Input("\t\tPredict protein amount by iBAQ\n");
		// Normalized
		//vector<double> vecNormalizedX;
		//vector<double> vecNormalizedY;
		//bootLinear.mf_Normalization(veclogiBAQOfStandProteins, vecNormalizedX, sdx); 
		//bootLinear.mf_Normalization(vecUPS2logMolsIdentified, vecNormalizedY,sdy);
		//bootLinear.mf_BootstrapLinearRegression(vecNormalizedX, vecNormalizedY, 10000, a, b,SdOfa,SdOfb, DevSQ, sd, RSS, Maxd, Mind, Meand);

		// No normalized
		bootLinear.mf_BootstrapLinearRegression(veclogiBAQOfStandProteins, vecUPS2logMolsIdentified, 10000, a, b, SdOfa, SdOfb, DevSQ, sd, RSS, Maxd, Mind, Meand);

		bootLinear.mf_Predict(veclogiBAQOfStandProteins, a, b, vecPredictedLogMolsOfStandProteinsByiBAQ);
		R2 = spearsonCorrelation(vecUPS2logMolsIdentified, vecPredictedLogMolsOfStandProteinsByiBAQ);
		cout << "\t\t\ta= " << a << endl;
		cout << "\t\t\tb= " << b << endl;
		cout << "\t\t\tSd of a= " << SdOfa << endl;
		cout << "\t\t\tSd of b= " << SdOfb << endl;
		cout << "\t\t\tMean Standard Error= " << sd << endl;
		cout << "\t\t\tAutocorrelation R2= " << R2*R2 << endl;
		flog.mf_Input("\t\t\ta= " + fDouble2String(a) + "\n");
		flog.mf_Input("\t\t\tb= " + fDouble2String(b) + "\n");
		flog.mf_Input("\t\t\tSd of a= " + fDouble2String(SdOfa) + "\n");
		flog.mf_Input("\t\t\tSd of b= " + fDouble2String(SdOfb) + "\n");
		flog.mf_Input("\t\t\tMean Standard Error= " + fDouble2String(sd) + "\n");
		flog.mf_Input("\t\t\tAutocorrelation R2= " + fDouble2String(R2*R2) + "\n");
		for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
		{
			if (m_bIfiBAQIntensityExist)
			{
				if (proteinIter->m_dMaxQuantiBAQ == 0.0)
				{
					proteinIter->m_dPredictedMolOfiBAQ = 0.0;
				}
				else
				{
					proteinIter->m_dPredictedMolOfiBAQ = pow(10, bootLinear.mf_Predict(log10(proteinIter->m_dMaxQuantiBAQ), a, b));
				}

			}
			else
			{
				if (proteinIter->m_dReCalculateiBAQ == 0.0)
				{
					proteinIter->m_dPredictedMolOfiBAQ = 0.0;
				}
				else
				{
					proteinIter->m_dPredictedMolOfiBAQ = pow(10, bootLinear.mf_Predict(log10(proteinIter->m_dReCalculateiBAQ), a, b));
				}

			}
		}


	}

}

void CProteinWorker::mf_PredictMolsByStandProteinsAndTopN(int iExperimentIndex)
{
	map<string, double> mapUPS2Id2Mols;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	vector<CProtein>::iterator proteinIter;
	string strTemp;

	cout << "Predict proteins' abundance by the standard proteins' TopN\n";
	flog.mf_Input("Predict proteins' abundance by the standard proteins' TopN\n");
	mf_LoadUPS2Mols(m_Params.m_strStandardProteinsFilePath, mapUPS2Id2Mols);


	double dPearsonCorr = 0.0;
	vector<double> veclogTop3OfStandProteins;
	vector<double> vecUPS2logMolsIdentified; // the spiked-in amount of identified standard proteins;
	vector<double> vecPredictedLogMolsOfStandProteinsByTop3;

	int iNumberofPeptidesOfStandProteins = 0;
	string::size_type stStart, stPosition, stEnd;
	vector<string> vecProteinsTemp;
	bool bIfOthersOnlyContam = true;
	map<string, int>::iterator mapLowSpikedInUPS2Iter;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		strTemp = proteinIter->m_strProteinID;
		stPosition = strTemp.find(m_Params.m_strIdentifierOfStandPro);
		vecProteinsTemp = split(proteinIter->m_strProteinFullName, ";");
		if (stPosition != strTemp.npos)
		{
			bIfOthersOnlyContam = true;
			if (m_Params.m_bIfExistContanProtein == true && vecProteinsTemp.size() > 1)
			{
				for (size_t i = 1; i < vecProteinsTemp.size(); i++)
				{
					if (vecProteinsTemp[i].find(m_Params.m_strContaminantfix) == vecProteinsTemp[i].npos)
					{
						bIfOthersOnlyContam = false;
					}
				}
			}
			if ((proteinIter->m_strProteinFullName.find(";") == \
				proteinIter->m_strProteinFullName.npos) || (bIfOthersOnlyContam == true))
			{ //
				Ups2MolsIter = mapUPS2Id2Mols.begin();
				for (; Ups2MolsIter != mapUPS2Id2Mols.end(); Ups2MolsIter++)
				{
					if (strTemp.find(Ups2MolsIter->first) != strTemp.npos)
					{
						iNumberofPeptidesOfStandProteins += proteinIter->m_vPeptidesSequences.size();
						if (Ups2MolsIter->second == 0.0)
						{
							vecUPS2logMolsIdentified.push_back(0.0);
						}
						else
						{
							vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
						}
						if (proteinIter->m_dProteinIntensityTop3 == 0.0)
						{
							veclogTop3OfStandProteins.push_back(0.0);
						}
						else
						{
							veclogTop3OfStandProteins.push_back(log10(proteinIter->m_dProteinIntensityTop3));
						}

						break;
					}
				} // for mapUPS2Id2Mols

				if (Ups2MolsIter == mapUPS2Id2Mols.end())
				{
					mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
					if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
					{
						cout << "Warning:\tCannot find protein " << strTemp << " in StandardProteins.txt\n";
						flog.mf_Input("Warning:\tCannot find protein " + strTemp + " in StandardProteins.txt\n");
					}
				}

			}
		} // end if (stPosition != strTemp.npos)
	}// end for m_vecProteins


	cout << "\tLoaded " << vecUPS2logMolsIdentified.size() << " standard proteins " << " with " << \
		iNumberofPeptidesOfStandProteins << " peptides\n";
	flog.mf_Input("\tLoaded " + fSize_t2String(vecUPS2logMolsIdentified.size()) + " standard proteins with " + \
		fInt2String(iNumberofPeptidesOfStandProteins) + " peptides\n");

	BootLineaRegress bootLinear;
	double a, b;
	double DevSQ;
	double sd;
	double RSS;
	double Maxd;
	double Mind;
	double Meand;
	double SdOfa;
	double SdOfb;
	double R2;

	cout << "\tPredict protein amount by bootstrap linear regression\n";
	flog.mf_Input("\tPredict protein amount by bootstrap linear regression\n");



	if (m_Params.m_bIfCalculateTop3)
	{

		cout << "\t\tPredict protein amount by Top3" << endl;
		flog.mf_Input("\t\tPredict protein amount by Top3\n");

		// No normalized
		bootLinear.mf_BootstrapLinearRegression(veclogTop3OfStandProteins, vecUPS2logMolsIdentified, 10000, a, b, SdOfa, SdOfb, DevSQ, sd, RSS, Maxd, Mind, Meand);

		bootLinear.mf_Predict(veclogTop3OfStandProteins, a, b, vecPredictedLogMolsOfStandProteinsByTop3);
		R2 = spearsonCorrelation(vecUPS2logMolsIdentified, vecPredictedLogMolsOfStandProteinsByTop3);
		cout << "\t\t\ta= " << a << endl;
		cout << "\t\t\tb= " << b << endl;
		cout << "\t\t\tSd of a= " << SdOfa << endl;
		cout << "\t\t\tSd of b= " << SdOfb << endl;
		cout << "\t\t\tMean Standard Error= " << sd << endl;
		cout << "\t\t\tAutocorrelation R2= " << R2*R2 << endl;
		flog.mf_Input("\t\t\ta= " + fDouble2String(a) + "\n");
		flog.mf_Input("\t\t\tb= " + fDouble2String(b) + "\n");
		flog.mf_Input("\t\t\tSd of a= " + fDouble2String(SdOfa) + "\n");
		flog.mf_Input("\t\t\tSd of b= " + fDouble2String(SdOfb) + "\n");
		flog.mf_Input("\t\t\tMean Standard Error= " + fDouble2String(sd) + "\n");
		flog.mf_Input("\t\t\tAutocorrelation R2= " + fDouble2String(R2*R2) + "\n");
		for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
		{
			if (proteinIter->m_dProteinIntensityTop3 == 0.0)
			{
				proteinIter->m_dPredictedMolOfTop3 = 0.0;
			}
			else
			{
				proteinIter->m_dPredictedMolOfTop3 = pow(10, bootLinear.mf_Predict(log10(proteinIter->m_dProteinIntensityTop3), a, b));
			}

		}

	}

}



int CProteinWorker::mf_GetExperimentNames(CQuantificationParam trainparam)
{
	ifstream fin(trainparam.m_strExperimentDesignPath.c_str());
	if (!fin)
	{
		cout << "Error:\tCannot open " << trainparam.m_strExperimentDesignPath << endl;
		flog.mf_Input("Error:\tCannot open " + trainparam.m_strExperimentDesignPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	if (trainparam.m_strIdentifySoftwareType == "maxquant")
	{
		map<string, int> mapExperiments;
		map<string, int>::iterator mapExperimentsIter;
		string strTemp1, strTemp2;
		int iTemp1 = 0, iTemp2 = 0;
		m_iNumberOfExperiments = 0;
		getline(fin, strTemp1);//jump the first row
		while (getline(fin, strTemp1))
		{
			iTemp1 = strTemp1.find("\t");
			iTemp1 = strTemp1.find("\t", iTemp1 + 1);
			iTemp2 = strTemp1.find("\n", iTemp1 + 1);
			strTemp2 = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			mapExperiments.insert(pair<string, int>(strTemp2, 0));
		}

		m_iNumberOfExperiments = mapExperiments.size();
		mapExperimentsIter = mapExperiments.begin();
		for (; mapExperimentsIter != mapExperiments.end(); mapExperimentsIter++)
		{
			m_vecStrExperiments.push_back(mapExperimentsIter->first);
		}
	}
	else if (trainparam.m_strIdentifySoftwareType == "mzQuantML")
	{
		map<string, string> mapAttributeAndvalues;
		map<string, string>::iterator mapAttributeAndvaluesIter;
		string strTemp1;
		string strAttributetemp, strValueTemp;
		int iTemp1 = 0, iTemp2 = 0;
		m_iNumberOfExperiments = 0;
		getline(fin, strTemp1);
		iTemp1 = strTemp1.find("\t");
		strAttributetemp = strTemp1.substr(0, iTemp1);
		iTemp2 = strTemp1.find("\n", iTemp1 + 1);
		strValueTemp = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		if (strAttributetemp == "Experiment Number")
			m_iNumberOfExperiments = atoi(strValueTemp.c_str());

		for (int i = 1; i <= m_iNumberOfExperiments; i++)
		{
			getline(fin, strTemp1);
			iTemp1 = strTemp1.find("\t");
			strAttributetemp = strTemp1.substr(0, iTemp1);
			iTemp2 = strTemp1.find("\n", iTemp1 + 1);
			strValueTemp = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			m_vecStrExperiments.push_back(strValueTemp);
		}
		fin.close();
		if (remove(trainparam.m_strExperimentDesignPath.c_str()) != 0)
		{
			cout << "Cannot remove file " << trainparam.m_strExperimentDesignPath << endl;
			flog.mf_Input("Cannot remove file " + trainparam.m_strExperimentDesignPath);
			exit(0);
		}
		
	}
	else if (trainparam.m_strIdentifySoftwareType == "PeakView")
	{
		map<string, string> mapAttributeAndvalues;
		map<string, string>::iterator mapAttributeAndvaluesIter;
		string strTemp1;
		string strAttributetemp, strValueTemp;
		int iTemp1 = 0, iTemp2 = 0;
		m_iNumberOfExperiments = 0;
		getline(fin, strTemp1);
		iTemp1 = strTemp1.find("\t");
		strAttributetemp = strTemp1.substr(0, iTemp1);
		iTemp2 = strTemp1.find("\n", iTemp1 + 1);
		strValueTemp = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		if (strAttributetemp == "Experiment Number")
			m_iNumberOfExperiments = atoi(strValueTemp.c_str());

		for (int i = 1; i <= m_iNumberOfExperiments; i++)
		{
			getline(fin, strTemp1);
			iTemp1 = strTemp1.find("\t");
			strAttributetemp = strTemp1.substr(0, iTemp1);
			iTemp2 = strTemp1.find("\n", iTemp1 + 1);
			strValueTemp = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			m_vecStrExperiments.push_back(strValueTemp);
		}
		fin.close();
		if (remove(trainparam.m_strExperimentDesignPath.c_str()) != 0)
		{
			cout << "Cannot remove file " << trainparam.m_strExperimentDesignPath << endl;
			flog.mf_Input("Cannot remove file " + trainparam.m_strExperimentDesignPath);
			exit(0);
		}
	}
	
	return m_iNumberOfExperiments;
}


bool CProteinWorker::mf_LoadUPS2Mols( string path, map<string, double>&mapUPS2Id2Mols)
{
	map<string, int>::iterator mapLowSpikedInUPS2Iter;
	m_mapLowSpikedInUPS2.clear();
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P08758ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P02741ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P05413ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P10145ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P02788ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P10636ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P00441ups", 1));
	m_mapLowSpikedInUPS2.insert(pair<string, int>("P01375ups", 1));

	mapUPS2Id2Mols.clear();

	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	string strTemp;
	FILE * pFile;
	double dMolsTemp;
	errno_t err;

	err = fopen_s(&pFile, path.c_str(),"r");
	string DefaultPathOfStandProteins = "StandardProteins.txt";
	if (err !=0 )
	{
		err = fopen_s(&pFile, DefaultPathOfStandProteins.c_str(), "r");
		if (err != 0)
		{
			cout << "Error:\tCannot open " << path << endl;
			flog.mf_Input("Error:\tCannot open " + path + "\n");
			flog.mf_Destroy();
			exit(1);
		}		
	}

	fgets(Buffer, BUFFERLENGTH, pFile);
	fgets(Buffer, BUFFERLENGTH, pFile);

	while (!feof(pFile))
	{
		if (Buffer[0] == '\n')
		{
			fgets(Buffer, BUFFERLENGTH, pFile);
			continue;
		}
		pstr = Buffer;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			flog.mf_Input("Error:\tThe format of \"StandardProteins.txt\" is wrong.\n");
			flog.mf_Destroy();
			exit(1);
		}
		strTemp = pstr;

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 == NULL)
		{
			pstr1 = strstr(pstr, "\n");
		}
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			flog.mf_Input("Error:\tThe format of \"StandardProteins.txt\" is wrong.\n");
			flog.mf_Destroy();
			exit(1);
		}
		dMolsTemp = atof(pstr);

		mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
		if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
		{
			mapUPS2Id2Mols.insert(pair<string, double>(strTemp, dMolsTemp));
		}
		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	fclose(pFile);
	return true;
}

bool CProteinWorker::mf_AnalysisPeptidesCVForSameUPS2Mols(string path, vector<CProtein>& vecProteins, map<string, double>&mapUPS2Id2Mols)
{
	ofstream ofileCVForSameMols(path);
	if (!ofileCVForSameMols)
	{
		cout << "Error:\tCannot open " <<path<< endl;
		flog.mf_Input("Error:\tCannot open " + path +"\n");
		flog.mf_Destroy();
		exit(1);
	}

	map<double, vector<double>> mapdUPS2Mols2PeptidesNativeIntensitys;
	map<double, vector<double>> mapdUPS2Mols2PeptidesCorrectedIntensitys;
	map<double, vector<double>>::iterator mapdUPS2Mols2PeptidesNativeIntensitysIter;
	map<double, vector<double>>::iterator mapdUPS2Mols2PeptidesCorrectedIntensitysIter;
	vector<double> vecdPeptideIntensityTemp;
	vector<CProtein>::iterator proteinIter;
	string::size_type stPosition;
	string strTemp;
	map<string, double>::iterator Ups2MolsIter;
	map<string, int>::iterator mapLowSpikedInUPS2Iter;

	int i = 0; // just for index
	for (proteinIter = vecProteins.begin(); proteinIter != vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateLFAQpep == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stPosition = strTemp.find(m_Params.m_strIdentifierOfStandPro);
			if (stPosition != strTemp.npos&&proteinIter->m_strProteinFullName.find(";") == proteinIter->m_strProteinFullName.npos)
			{ //
				Ups2MolsIter = mapUPS2Id2Mols.begin();
				for (; Ups2MolsIter != mapUPS2Id2Mols.end(); Ups2MolsIter++)
				{
					if (strTemp.find(Ups2MolsIter->first) != strTemp.npos)
					{
						mapdUPS2Mols2PeptidesNativeIntensitysIter = mapdUPS2Mols2PeptidesNativeIntensitys.find(Ups2MolsIter->second);
						mapdUPS2Mols2PeptidesCorrectedIntensitysIter = mapdUPS2Mols2PeptidesCorrectedIntensitys.find(Ups2MolsIter->second);
						if (mapdUPS2Mols2PeptidesNativeIntensitysIter == mapdUPS2Mols2PeptidesNativeIntensitys.end())
						{//new mols number group;
							vecdPeptideIntensityTemp.clear();
							vecdPeptideIntensityTemp = proteinIter->m_vPeptidesNativeIntensity;
							mapdUPS2Mols2PeptidesNativeIntensitys.insert(pair<double, vector<double>>(Ups2MolsIter->second, vecdPeptideIntensityTemp));

							vecdPeptideIntensityTemp.clear();
							vecdPeptideIntensityTemp = proteinIter->m_vPeptidesCorrectedIntensity;
							mapdUPS2Mols2PeptidesCorrectedIntensitys.insert(pair<double, vector<double>>(Ups2MolsIter->second, vecdPeptideIntensityTemp));

						}
						else
						{
							for (i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
							{
								mapdUPS2Mols2PeptidesNativeIntensitysIter->second.push_back(proteinIter->m_vPeptidesNativeIntensity.at(i));
							}
							for (i = 0; i < proteinIter->m_vPeptidesCorrectedIntensity.size(); i++)
							{
								mapdUPS2Mols2PeptidesCorrectedIntensitysIter->second.push_back(proteinIter->m_vPeptidesCorrectedIntensity.at(i));
							}
						}//end if new mols
						break;
					}
				} // for mapUPS2Id2Mols
				if (Ups2MolsIter == mapUPS2Id2Mols.end())
				{
					mapLowSpikedInUPS2Iter = m_mapLowSpikedInUPS2.find(strTemp);
					if (mapLowSpikedInUPS2Iter == m_mapLowSpikedInUPS2.end())
					{
						cout << "\tWarning:\tCannot find protein " << strTemp << " in StandardProteins.txt\n";
						flog.mf_Input("\tWarning:\t Cannot find protein " + strTemp + " in StandardProteins.txt\n");
					}
				}				
			}//end if UPS2 exist
		}// end if (proteinIter->m_bIfCalculateLFAQpep == true)
	} //end for proteins

	mapdUPS2Mols2PeptidesNativeIntensitysIter = mapdUPS2Mols2PeptidesNativeIntensitys.begin();
	mapdUPS2Mols2PeptidesCorrectedIntensitysIter = mapdUPS2Mols2PeptidesCorrectedIntensitys.begin();

	for (; mapdUPS2Mols2PeptidesNativeIntensitysIter != mapdUPS2Mols2PeptidesNativeIntensitys.end(); mapdUPS2Mols2PeptidesNativeIntensitysIter++, mapdUPS2Mols2PeptidesCorrectedIntensitysIter++)
	{
		ofileCVForSameMols << mapdUPS2Mols2PeptidesNativeIntensitysIter->first << "\t";
		ofileCVForSameMols << CalculatelogCV(mapdUPS2Mols2PeptidesNativeIntensitysIter->second) << "\t";
		ofileCVForSameMols << CalculatelogCV(mapdUPS2Mols2PeptidesCorrectedIntensitysIter->second) << endl;
	}

	ofileCVForSameMols.close();
	return true;
}