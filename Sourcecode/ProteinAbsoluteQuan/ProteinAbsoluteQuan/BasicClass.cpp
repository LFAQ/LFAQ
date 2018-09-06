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
#pragma once
#include"stdafx.h"
#include"BasicClass.h"


int UniquePeptidesTrainThreshold=5; // Default value
CProtein::CProtein()
{
	m_dPeptidesIntensityMax = 0.0;
	m_dPeptidesIntensityMedian = 0.0;
	m_dProteinLFAQpep = 0.0;
	m_dProteinLFAQpro = 0.0;
	m_dProteinQfactor = 0.0;
	m_dMaxQuantiBAQ = 0.0;
	m_dReCalculateiBAQ = 0.0;
	m_iPeptidesNumber = 0;
	m_dPeptidesIntensitySum = 0.0;
	m_dCVNative = 0.0;
	m_dCVAfterCorrected = 0.0;
	m_bIfCalculateLFAQpep = false;
	m_bIfCalculateLFAQpro = false;
	m_dPredictedMolOfLFAQpep = 0.0;
	m_dPredictedMolOfLFAQpro=0.0;
	m_dPredictedMolOfiBAQ=0.0;
	m_dProteinIntensityTop3 = 0.0;
	m_dPredictedMolOfTop3=0.0;
	m_bIfInTrainSet = false;
}
void CProtein::Clear()
{
	m_vPeptidesIDs.clear();
	m_vPeptidesSequencesWithFlankingRegion.clear();
	m_vPeptidesSequences.clear();
	m_vPeptidesFeatures.clear();
	m_vPeptideExperimentalQFactors.clear();
	m_vPeptidesNativeIntensity.clear();
	m_vPeptidesCombinedIntensity.clear();
	m_vPeptidesCorrectedIntensity.clear();
	m_dPeptidesIntensityMedian = 0.0;
	m_dPeptidesIntensityMax = 0.0;
	m_mapPeptidesSequenceAndCorrespondRawFiles.clear();
	m_mapExperimentAndPeptidesIntensity.clear();
	m_vPeptideQfactors.clear();
	m_vecdMaxQuantiBAQ.clear();
	m_dMaxQuantiBAQ = 0.0;
	m_vdPeptidesMW.clear();
	m_dReCalculateiBAQ = 0.0;
	m_bIfInTrainSet = false;
}

string CProtein::mf_GetPeptidesAdjacentSequence(string PeptideSequence)
{
	string::size_type peptideLocation = 0, begin = 0, Getwidth = 0, end;
	end = PeptideSequence.find("{");
	PeptideSequence = PeptideSequence.substr(0, end);
	peptideLocation = m_strProteinSequence.find(PeptideSequence);
	if (peptideLocation != m_strProteinSequence.npos)
	{
		if ((peptideLocation >= FlankingRegionWidth) && (peptideLocation <= m_strProteinSequence.size() - PeptideSequence.size() - FlankingRegionWidth))
		{
			begin = peptideLocation - FlankingRegionWidth;
			Getwidth = PeptideSequence.size() + 2 * FlankingRegionWidth;
		}
		else if (peptideLocation < FlankingRegionWidth)
		{
			begin = 0;
			Getwidth = peptideLocation + PeptideSequence.size() + FlankingRegionWidth;
		}
		else if (peptideLocation > m_strProteinSequence.size() - PeptideSequence.size() - FlankingRegionWidth)
		{
			begin = peptideLocation - FlankingRegionWidth;
			Getwidth = m_strProteinSequence.size() - peptideLocation - 1 + PeptideSequence.size() + FlankingRegionWidth;
		}
		return m_strProteinSequence.substr(begin, Getwidth);
	}
	else
	{
		cout << "\tWarning:\tCannot find peptide " << PeptideSequence << " in protein " << m_strProteinID << endl;
		cout << "The protein sequence is " << m_strProteinSequence << endl; 
		flog.mf_Input("\tWarning:\tCannot find peptide " + PeptideSequence + " in protein " + m_strProteinID +"\n");
		flog.mf_Input("The protein sequence is " + m_strProteinSequence + "\n");
		return "NULL";
	}
}
// judge if RelatedProtein A is the subset of RelatedProtein B;
bool CRelatedProtein::mf_ifSubset(CRelatedProtein& A)
{
	vector<string> vecA;
	vector<string> vecB;
	for (size_t i = 0; i < A.vecPeptides.size(); i++)
	{
		vecA.push_back(A.vecPeptides[i]);
	}
	for (size_t j = 0; j < this->vecPeptides.size(); j++)
	{
		vecB.push_back(this->vecPeptides[j]);
	}

	if (ifSubset(vecA,vecB))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void CRelatedProtein::Clear()
{
	bIfSubSet = false;
	strProteinID="";
	vecPeptides.clear();
	vecSubsetProteinIDs.clear(); 
	vecBIfPeptideGroupUnique.clear();
}
void Cpeptides::mf_Clear()
{
	m_mapPeptideIDAndSequence.clear();
	m_mapPeptideIDAndIntensity.clear();
	m_mapPeptideSequenceAndIntensity.clear();
	m_mapPeptideSequencuAndIfUse.clear();
	m_mapPeptideSequenceAndAttribute.clear();
	m_mapPeptideIDAndIntensitys.clear();
}

void CMergedProtein::Clear()
{
	m_vecExperiments.clear();
	m_vecMaxQuantiBAQOfExperiments.clear();
	m_vecReCaliBAQOfExperiments.clear();
	m_vecLFAQproOfExperiments.clear();
	m_vecLFAQpepOfExperiments.clear();
	m_vecTop3OfExperiments.clear();
	m_vecPredictedMolsByLFAQpepOfExperiments.clear();
	m_vecPredictedMolsByLFAQproOfExperiments.clear();
	m_vecPredictedMolsByiBAQOfExperiments.clear();
	m_vecPredictedMolsByTop3OfExperiments.clear();
	m_mapExperimentsAndPeptidesIntensity.clear();
}

void CMergedProteins::mf_ExperimentsAnalysis()
{
	double dCVTemp;
	for (size_t i = 0; i < m_vecMergedProteins.size(); i++)
	{		
		dCVTemp = CalculatelogCV(m_vecMergedProteins[i].m_vecMaxQuantiBAQOfExperiments);
		m_vecMaxQuantiBAQCV.push_back(dCVTemp);
		dCVTemp = CalculatelogCV(m_vecMergedProteins[i].m_vecReCaliBAQOfExperiments);
		m_vecReCaliBAQCV.push_back(dCVTemp);
		dCVTemp = CalculatelogCV(m_vecMergedProteins[i].m_vecLFAQpepOfExperiments);
		m_vecLFAQpepCV.push_back(dCVTemp);
		dCVTemp = CalculatelogCV(m_vecMergedProteins[i].m_vecLFAQproOfExperiments);
		m_vecLFAQproCV.push_back(dCVTemp);
		dCVTemp = CalculatelogCV(m_vecMergedProteins[i].m_vecTop3OfExperiments);
		m_vecTop3CV.push_back(dCVTemp);
	}
}

string AddDecorateStar(const string& strProjectName)
{
	string strTemp;
	for (int i = 0; i < 80; i++)
	{
		strTemp += "*";
	}
	strTemp += "\n";
	for (int i = 0; i < 25; i++)
	{
		strTemp += " ";
	}
	strTemp += strProjectName;
	strTemp += "\n";
	for (int i = 0; i < 80; i++)
	{
		strTemp += "*";
	}
	strTemp += "\n";
	return strTemp;
}
string Unicode2Multibyte(_TCHAR * chars)
{
	DWORD dwNum = WideCharToMultiByte(CP_OEMCP, NULL, chars, -1, NULL, 0, NULL, FALSE);
	char *psText;
	psText = new char[dwNum];
	if (!psText)
	{
		delete[]psText;
	}
	WideCharToMultiByte(CP_OEMCP, NULL, chars, -1, psText, dwNum, NULL, FALSE);
	string strTemp = psText;
	delete[]psText;
	psText = NULL;
	return strTemp;
}
map<string, string> GetParametersFromFile(const string &ParamFilePath)
{
	map<string,string> mapAttributeNameAndValues;

	ifstream fin(ParamFilePath.c_str());
	if (!fin)
	{
		cout << "Error:\tCannot open parameter file: " << ParamFilePath << endl;
		if (flog.m_bIfReady)
		{
			flog.mf_Input("Error:\tCannot open parameter file: " + ParamFilePath + ".\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	string strAttributeName, strAtributeValue;
	string strTemp = "";
	int iTemp1 = 0, iTemp2 = 0;
	while (getline(fin, strTemp))
	{
		if (strTemp == "")
		{
			continue;
		}
		iTemp2 = strTemp.find("=", 0);
		strAttributeName = strTemp.substr(0, iTemp2);
		iTemp1 = strTemp.find("\"");
		if (iTemp1 == strTemp.npos)
		{
			cout << "Error:\tCannot read parameter from " << strTemp << endl;
			if (flog.m_bIfReady)
			{
				flog.mf_Input("Error:\tCannot read parameter from " + strTemp + ".\n");
				flog.mf_Destroy();
				exit(1);
			}
			return mapAttributeNameAndValues;
		}
		iTemp2 = strTemp.find("\"", iTemp1 + 1);
		if (iTemp2 == strTemp.npos)
		{
			cout << "Error:\tCannot read parameter from " << strTemp << endl;
			if (flog.m_bIfReady)
			{
				flog.mf_Input("Error:\tCannot read parameter from " + strTemp + ".\n");
				flog.mf_Destroy();
				exit(1);
			}
			return mapAttributeNameAndValues;
		}
		strAtributeValue = strTemp.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		transform(strAttributeName.begin(), strAttributeName.end(), strAttributeName.begin(), ::tolower);//change AttributeNameto lowcase
		mapAttributeNameAndValues.insert(pair<string, string>(strAttributeName, strAtributeValue));
	}
	return mapAttributeNameAndValues;
}
string GetResultPath(string ParamFilePath)
{
	map<string, string> mapAttributeNameAndValues;
	map<string, string>::iterator mapAttributeNameAndValuesIter;
	mapAttributeNameAndValues = GetParametersFromFile(ParamFilePath);	
	mapAttributeNameAndValuesIter = mapAttributeNameAndValues.find("resultpath");
	string strAtributeValue;
	if (mapAttributeNameAndValuesIter != mapAttributeNameAndValues.end())
	{
		return mapAttributeNameAndValuesIter->second;
	}
	else
	{
		return "NULL";
	}
}

// check if the file exist and if it contain double quotation marks
bool CheckFilePath(string& filepath)
{
	if (filepath.at(0) == '\"'&&filepath.at(filepath.size() - 1) == '\"')
	{
		filepath = filepath.substr(1, filepath.size() - 2);
	}
	fstream fc(filepath);
	if (!fc)
	{
		return false;
	}
	return true;
}

bool bPrefix(string str, string prefix)
{
	if (prefix.size() > str.size())
	{
		return false;
	}
	if (prefix.empty())
	{
		return true;
	}
	return str.compare(0, prefix.size(), prefix) == 0;
}

string stringTrim(string str)
{
	//search for the begin of truncated string
	string::iterator begin = str.begin();
	while (begin != str.end() && (*begin == ' ' || *begin == '\t' || *begin == '\n' || *begin == '\r'))
	{
		++begin;
	}

	//all characters are whitespaces
	if (begin == str.end())
	{
		str.clear();
		return str;
	}

	//search for the end of truncated string
	string::iterator end = str.end();
	end--;
	while (end != begin && (*end == ' ' || *end == '\n' || *end == '\t' || *end == '\r'))
	{
		--end;
	}
	++end;

	//no characters are whitespaces
	if (begin == str.begin() && end == str.end())
	{
		return str;
	}

	str=string(begin, end);

	return str;
}

string removestringWhitespaces(string str)
{
	bool contains_ws = false;
	for (string::const_iterator it = str.begin(); it != str.end(); ++it)
	{
		char c = *it;
		if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
		{
			contains_ws = true;
			break;
		}
	}

	if (contains_ws)
	{
		string tmp;
		tmp.reserve(str.size());
		for (string::const_iterator it = str.begin(); it != str.end(); ++it)
		{
			char c = *it;
			if (c != ' ' && c != '\t' && c != '\n' && c != '\r')
				tmp += c;
		}
		str.swap(tmp);
	}

	return str;
}

string stringtoLower(string str)
{
	std::transform(str.begin(), str.end(), str.begin(), (int(*)(int))tolower);
	return str;

}
string PrefixOfstring(string str, char delim)
{
	size_t pos = str.find(delim);
	if (pos == string::npos) //char not found
	{
		cout << "Error:\tCannot find delimiter " << delim << endl;
		flog.mf_Input("Error:\tCannot find delimiter " + delim );
		flog.mf_Input("\n");
		flog.mf_Destroy();
	}
	return str.substr(0, pos);
}
string SuffixOfstring(string str, char delim)
{
	size_t pos = str.rfind(delim);
	if (pos == string::npos) //char not found
	{
		cout << "Error:\tCannot find delimiter " << delim << endl;
		flog.mf_Input("Error:\tCannot find delimiter " + delim);
		flog.mf_Input("\n");
		flog.mf_Destroy();
	}
	return str.substr(++pos);
}
bool stringSplit(string str, const string splitter, vector<string>& substrings)
{
	substrings.clear();
	if (str.empty())
		return false;

	if (splitter.empty()) // split after every character:
	{
		substrings.resize(str.size());
		for (size_t i = 0; i < str.size(); ++i)
			substrings[i] = str[i];
		return true;
	}

	size_t len = splitter.size(), start = 0, pos = str.find(splitter);
	if (len == 0)
		len = 1;
	while (pos != string::npos)
	{
		substrings.push_back(str.substr(start, pos - start));
		start = pos + len;
		pos = str.find(splitter, start);
	}
	substrings.push_back(str.substr(start, str.size() - start));
	return substrings.size() > 1;
}
double Average(double *Array, int Number)
{
	double ave = 0.0;
	for (int i = 0; i < Number; i++)
	{
		ave += Array[i];
	}
	return ave / Number;

}
double Average(vector<double> Array)
{
	if (Array.size() == 0)
	{
		return 0.0;
	}
	double ave = 0.0;
	vector<double>::iterator doubleIter;
	for (doubleIter = Array.begin(); doubleIter != Array.end(); doubleIter++)
	{
		ave += *doubleIter;
	}
	double len = Array.size();
	return ave / len;
}

void quickSort(double v[], int indexTemp[], int left, int right)
{ // A problem: Stack overflow may happen in the sort function when the dataset is too big.
	if (left < right)
	{
		double key = v[left];
		int indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		quickSort(v, indexTemp, left, low - 1);
		quickSort(v, indexTemp, low + 1, right);
	}
}

void quickSort(vector<size_t> &v, vector<size_t> &indexTemp, int left, int right)
{
	if (left < right)
	{
		size_t key = v[left];
		size_t indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		quickSort(v, indexTemp, left, low - 1);
		quickSort(v, indexTemp, low + 1, right);
	}
}

void quickSort(vector<double> &v, vector<size_t> &indexTemp, int left, int right)
{
	if (left < right)
	{
		double key = v[left];
		size_t indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		quickSort(v, indexTemp, left, low - 1);
		quickSort(v, indexTemp, low + 1, right);
	}
}


double spearsonCorrelation(double s1[], double s2[], int NumberOfSamples)
{   
	double dCorrelation = 0.0;
	double means2 = 0.0;
	double means1 = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	means1 = Average(s1, NumberOfSamples);
	for (int i = 0; i < NumberOfSamples; i++)
		means2 += s2[i];
	means2 = means2 / NumberOfSamples;
	for (int i = 0; i < NumberOfSamples; i++)
		dCorrelation += (s2[i] - means2)*(s1[i] - means1);
	for (int i = 0; i < NumberOfSamples; i++)
	{
		sd1 += (s1[i] - means1)*(s1[i] - means1);
		sd2 += (s2[i] - means2)*(s2[i] - means2);
	}
	dCorrelation = dCorrelation / sqrt(sd1*sd2);
	if (sd1*sd2 <= 0) 
		return -1;
	else	return dCorrelation;
}

double spearsonCorrelation(vector<double> s1, double s2[])
{
	double dCorrelation = 0.0;
	double means2 = 0.0;
	double means1 = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	int NumberOfSamples = s1.size();
	means1 = Average(s1);
	means2 = Average(s2, NumberOfSamples);
	for (int i = 0; i < NumberOfSamples; i++)
		dCorrelation += (s2[i] - means2)*(s1[i] - means1);
	for (int i = 0; i < NumberOfSamples; i++)
	{
		sd1 += (s1[i] - means1)*(s1[i] - means1);
		sd2 += (s2[i] - means2)*(s2[i] - means2);
	}
	dCorrelation = dCorrelation / sqrt(sd1*sd2);
	if (sd1*sd2 <= 0) 
		return -1;
	else	return dCorrelation;
}

double spearsonCorrelation(vector<double> s1, vector<double> s2)
{
	double dCorrelation = 0.0;
	double means2 = 0.0;
	double means1 = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	int NumberOfSamples = s1.size();
	means1 = Average(s1);
	means2 = Average(s2);
	for (int i = 0; i < NumberOfSamples; i++)
		dCorrelation += (s2[i] - means2)*(s1[i] - means1);
	for (int i = 0; i < NumberOfSamples; i++)
	{
		sd1 += (s1[i] - means1)*(s1[i] - means1);
		sd2 += (s2[i] - means2)*(s2[i] - means2);
	}
	dCorrelation = dCorrelation / sqrt(sd1*sd2);
	if (sd1*sd2 <= 0) 
		return 0;
	else	return dCorrelation;
}

double CalculatelogCV(const vector<double>&s)
{  // add log 20170607
	vector<double> nonZeroS;
	nonZeroS.clear();
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] != 0.0)
		{
			nonZeroS.push_back(s[i]);
		}
	}

	double mean = 0.0;
	double sd = 0.0, sd2 = 0.0;
	int len = nonZeroS.size();


	for (int i = 0; i < len; i++)
	{
		nonZeroS[i] = log10(nonZeroS[i]);
		mean += (nonZeroS[i]);
	}

	if (len <= 1)
		return 0.0;
	else
	{
		mean = mean / len;
		for (int i = 0; i < len; i++)
		{
			sd += ((nonZeroS[i]) - mean)*((nonZeroS[i]) - mean);
		}
	}

	sd = sqrt(sd / len);

	if (mean <= 0)
		return -1;
	else	return sd / mean;
}

void Discretization(double *TestY, int m_iTestSampleNumber, int iBinNum)
{
	vector <pair<double, int>> vctemp;
	int bin_volume = m_iTestSampleNumber / iBinNum;
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		vctemp.push_back(make_pair(TestY[i], i));
	}
	sort(vctemp.begin(), vctemp.end());
	double dBinGap = 1 / (double)iBinNum;
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		TestY[vctemp[i].second] = int(i / bin_volume + 1)*dBinGap;
		if (TestY[vctemp[i].second]>1) 
			TestY[vctemp[i].second] = 1;
	}
}

void Discretization(vector<double> &TestY,int iBinNum)
{
	vector <pair<double, int>> vctemp;
	int m_iTestSampleNumber = TestY.size();
	int bin_volume = m_iTestSampleNumber / iBinNum;
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		vctemp.push_back(make_pair(TestY[i], i));
	}
	sort(vctemp.begin(), vctemp.end());
	double dBinGap = 1 / (double)iBinNum;
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		TestY[vctemp[i].second] = int(i / bin_volume + 1)*dBinGap;
		if (TestY[vctemp[i].second]>1) 
			TestY[vctemp[i].second] = 1;
	}
}
double GetMaxValue(const vector<double> PeptidesIntensity)
{
	double dMaxIntensitytemp = 0.0;
	for (size_t i = 0; i < PeptidesIntensity.size(); i++)
	{
		if (dMaxIntensitytemp < PeptidesIntensity.at(i))
		{
			dMaxIntensitytemp = PeptidesIntensity.at(i);
		}

	}
	return dMaxIntensitytemp;
}
double GetMaxValue(double *PeptidesIntensity, int iNum)
{
	double dMaxIntensitytemp = 0.0;
	for (size_t i = 0; i < iNum; i++)
	{
		if (dMaxIntensitytemp < PeptidesIntensity[i])
		{
			dMaxIntensitytemp = PeptidesIntensity[i];
		}
	}
	return dMaxIntensitytemp;
}

double GetMinValue(vector<double> PeptidesIntensity)
{
	double dMinTemp = 1e15;
	for (size_t i = 0; i < PeptidesIntensity.size(); i++)
	{
		if (dMinTemp > PeptidesIntensity.at(i))
		{
			dMinTemp = PeptidesIntensity.at(i);
		}
	}
	return dMinTemp;
}
double GetMinValue(double *PeptidesIntensity, int iNum)
{
	double dMinTemp = 1e15;
	for (size_t i = 0; i < iNum; i++)
	{
		if (dMinTemp > PeptidesIntensity[i])
		{
			dMinTemp = PeptidesIntensity[i];
		}
	}
	return dMinTemp;
}

double GetSum(vector<double> vec)
{
	double dSumTemp = 0.0;
	for (size_t i = 0; i < vec.size(); i++)
	{
		dSumTemp += vec[i];
	}
	return dSumTemp;
}
int GetMaxIndex(vector<double> vec)
{
	double dMaxIntensitytemp = 0.0;
	int iMaxIndex = 0;
	for (size_t i = 0; i < vec.size(); i++)
	{
		if (dMaxIntensitytemp < vec.at(i))
		{
			dMaxIntensitytemp = vec.at(i);
			iMaxIndex = i;
		}
	}
	return iMaxIndex;
}

double GetMedianIntensity(vector<double> PeptidesIntensity)
{
	sort(PeptidesIntensity.begin(), PeptidesIntensity.end());
	size_t len = PeptidesIntensity.size();
	if (len == 0)             
		return 0.0;
	if (len % 2 == 0)
		return (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at((len - 2) / 2)) / 2;
	else
		return PeptidesIntensity.at(floor(len / 2)); 
}

//Ensure that there is not negative and zeros;
void Assertion(vector<double> &TestY)
{
	for (size_t i = 0; i < TestY.size(); i++)
	{
		if (TestY.at(i) <= 0.1)
			TestY.at(i) = 0.1;
		if (TestY.at(i) > 1.0)
			TestY.at(i) = 1.0;
	}
}

//Ensure that there is not negative and zeros;
void Assertion(double *TestY, int m_iTestSampleNumber)
{
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		if (TestY[i] <= 0.1)
		{
			TestY[i] = 0.1;
		}
		if (TestY[i] > 1.0)
			TestY[i] = 1.0;
	}
}
void Scaling(vector<double> &TestY)
{
	double dMax = GetMaxValue(TestY);
	double dMin = GetMinValue(TestY);
	for (size_t i = 0; i < TestY.size(); i++)
	{
		TestY[i] = ((TestY[i] - dMin)*0.9) / dMax + 0.1;
	}
}
void Scaling(double *TestY, int m_iTestSampleNumber)
{
	double dMax = GetMaxValue(TestY,m_iTestSampleNumber);
	double dMin = GetMinValue(TestY, m_iTestSampleNumber);
	double t = 0.05;
	for (size_t i = 0; i < m_iTestSampleNumber; i++)
	{
		TestY[i] = ((TestY[i] - dMin)*(1-t)) / (dMax-dMin) + t;
	}
}
//integral using trapzoid formular
double IntegralFromVector(vector<double> Values, vector<double> times)
{
	double result = 0.0;
	if (Values.size() != times.size())
	{
		cout << "Error:\tThe size of Values does not correspond to times'" << endl;
		flog.mf_Input("Error:\tThe size of Values does not correspond to times'\n");
		flog.mf_Destroy();
		exit(1);
	}
	for (size_t i = 1; i < times.size(); i++)
	{
		result = result + Values.at(i)*(times.at(i) - times.at(i - 1)) / 2;
	}
	return result;
}

void GetQuantiles(vector<double> PeptidesIntensity, double &FirstQ, double &median, double &ThirdQ)
{
	sort(PeptidesIntensity.begin(), PeptidesIntensity.end());
	int len = PeptidesIntensity.size();
	if (len ==0)            
	{
		FirstQ = 0.0;
		median = 0.0;
		ThirdQ = 0.0;
	}
	else if (len == 1)
	{
		FirstQ = PeptidesIntensity.front();
		median = PeptidesIntensity.front();
		ThirdQ = PeptidesIntensity.front();
	}
	else if (len == 2)
	{
		FirstQ = PeptidesIntensity.at(0);
		median = (PeptidesIntensity.at(0) + PeptidesIntensity.at(1)) / 2;
		ThirdQ = PeptidesIntensity.at(1);
	}
	else if (len == 3)
	{
		FirstQ = PeptidesIntensity.at(0);
		median = PeptidesIntensity.at(1);
		ThirdQ = PeptidesIntensity.at(2);
	}
	else if (len % 4 == 0)
	{
		FirstQ = (PeptidesIntensity.at(len/4)+PeptidesIntensity.at(len/4-1))/2;
		median = (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at(len /2-1))/ 2;
		ThirdQ = (PeptidesIntensity.at(3 * (len / 4)) + PeptidesIntensity.at(3 * (len / 4) - 1)) / 2;
	} 
	else if (len % 4 == 1)
	{
		FirstQ = (PeptidesIntensity.at(floor(len / 4)) + PeptidesIntensity.at(ceil(len / 4))) / 2;
		median = PeptidesIntensity.at(floor(len / 2));
		ThirdQ = (PeptidesIntensity.at(floor(3 * (len / 4))) + PeptidesIntensity.at(ceil(3 * (len / 4) ))) / 2;
	}
	else if (len % 4 == 2)
	{
		FirstQ = PeptidesIntensity.at(floor(len / 4)) ;
		median = (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at(len / 2 - 1)) / 2;
		ThirdQ = PeptidesIntensity.at(floor(3 * (len / 4))) ;

	}
	else if (len % 4 == 3)
	{
		FirstQ = PeptidesIntensity.at(floor(len / 4));
		median = PeptidesIntensity.at(ceil(len / 2));
		ThirdQ = PeptidesIntensity.at(floor(3 * (len / 4)));
	}
}

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns, string sep)
{
    char *pstr1;
    int icolumns = 0;

    pstr1 = strstr(pstr, sep.c_str());
    while (pstr1 != NULL)
    {
        if (pstr1 != NULL)
        {
            *pstr1 = '\0';
        }
        else
        {
            cout << "Error:\tThe format is wrong.\n";
            exit(1);
        }
        mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
        icolumns++;
        pstr = pstr1 + 1;
        pstr1 = strstr(pstr, sep.c_str());
    }
    pstr1 = strstr(pstr, "\n");
    if (pstr1 != NULL)
    {
        *pstr1 = '\0';
    }
    else
    {
        cout << "Error:\tThe format is wrong.\n";
        exit(1);
    }
    mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
}
void GetAttributesFromFirstRow(string str, map<string, int>& mapAttributesAndColumns, string sep)
{
	mapAttributesAndColumns.clear();
	if (str == "")
	{
		return;
	}
	int iColumn = 0;
	int iBegin = 0, iEnd = 0;
	iEnd = str.find(sep, iBegin);
	string strsub;
	while (iEnd!=str.npos && iBegin < str.size() - 1)
	{
		strsub = str.substr(iBegin, iEnd - iBegin);
		mapAttributesAndColumns.insert(pair<string, int>(strsub, iColumn));
		iColumn++;
		iBegin = iEnd + 1;
		iEnd = str.find(sep, iBegin);
	}
	if (iBegin < str.size() - 1)
	{
		iEnd = str.size();
		strsub = str.substr(iBegin, iEnd - iBegin);
		mapAttributesAndColumns.insert(pair<string, int>(strsub, iColumn));
	}
}


bool fStringToBool(string str, bool &bl)
{
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	if (str == "true")
		bl = true;
	else if (str == "false")
		bl= false;
	else
	{
		return false;
	}
	return true;
}
string fInt2String(int i)
{
	char buf[10];
	sprintf(buf, "%d", i);
	string str = buf;
	return str;
}
string fSize_t2String(size_t si)
{
	int iTemp = si;
	string str = fInt2String(iTemp);
	return str;
}
string fDouble2String(double d)
{
	ostringstream os;
	os << d;
	return os.str();
}
bool is_dir(const char* path) {
	struct _stat buf = { 0 };
	_stat(path, &buf);
	return buf.st_mode & _S_IFDIR;
}
bool is_file(const char* FILENAME)
{
	fstream _file;
	_file.open(FILENAME, ios::in);
	if (!_file)
	{
		return false;
	}
	else
	{
		return true;
	}

}

//remove the whitespace at the end of the string or at the begin of the string
string strim(string &str)
{
	int ilocation = str.find_last_of(' ');
	if (ilocation != -1)
	{
		//remove the whitespace at the end of the string
		str.erase(str.find_last_not_of(' ') + 1);    
	}
	ilocation = str.find_first_not_of(' ');
	//remove the whitespace at the begin of the string
	str.erase(0, ilocation);    
	return str;
}

std::vector<string> split(const string &strLine, const string& splitter)
{
	std::vector<string> substrings;
	substrings.clear();
	string str = strLine;
	//delete the newline character if it has
	if (str.at(str.size() - 1) == '\n')
	{
		str = str.substr(0, str.size() - 1);
	}
	// split after every character:
	if (splitter.empty()) 
	{
		substrings.resize(str.size());
		for (size_t i = 0; i < str.size(); ++i)
			substrings[i] = str[i];
		return substrings;
	}

	size_t len = splitter.size(), start = 0, pos = str.find(splitter);
	if (len == 0)
		len = 1;
	while (pos != string::npos)
	{
		substrings.push_back(str.substr(start, pos - start));
		start = pos + len;
		pos = str.find(splitter, start);
	}
	substrings.push_back(str.substr(start, str.size() - start));
	return substrings;
}

void SaveVector(string path, vector<double> vec)
{
	ofstream oFile(path);
	if (!oFile)
	{
		cout << "Cannot open " << path << endl;
		exit(-1);
	}

	for (int i = 0; i < vec.size(); i++)
	{
		oFile << vec.at(i) << endl;
	}
	oFile.close();
}
void SaveVector(string path, vector<string> vec)
{
	ofstream oFile(path);
	if (!oFile)
	{
		cout << "Cannot open " << path << endl;
		exit(-1);
	}

	for (int i = 0; i < vec.size(); i++)
	{
		oFile << vec.at(i) << endl;
	}
	oFile.close();
}

void SaveArray(string path, double * arr, int len)
{
	ofstream oFile(path);
	if (!oFile)
	{
		cout << "Cannot open " << path << endl;
		exit(-1);
	}

	for (int i = 0; i < len; i++)
	{
		oFile << arr[i] << endl;
	}
	oFile.close();
}
// Judge if vector A is the subset of vector B;
bool ifSubset(vector<string>&vecA, vector<string>& vecB)
{

	int iNumA = vecA.size();
	int iNumB = vecB.size();
	if (iNumA > iNumB)
	{
		return false;
	}
	vector<string>A;
	vector<string>B;
	for (int i = 0; i < iNumA; i++)
	{
		A.push_back(vecA[i]);
	}
	for (int j = 0; j <iNumB; j++)
	{
		B.push_back(vecB[j]);
	}

	sort(A.begin(), A.end());
	sort(B.begin(), B.end());

	int iA = 0, iB = 0;
	while (iA < iNumA&&iB < iNumB)
	{
		if (B[iB] < A[iA])
		{
			iB++;
		}
		else if (A[iA] == B[iB])
		{
			iA++;
			iB++;
		}
		else if (B[iB]>A[iA])
		{
			return false;
		}
	}
	if (iA < iNumA)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool MoreThan(double a, double b)
{
	return a > b;
}

double round(double number, unsigned int bits)
{
	double integerPart = floor(number);
	number -= integerPart;
	for (unsigned int i = 0; i < bits; i++)
	{
		number *= 10;
	}
	number = floor(number + 0.5);
	for (unsigned int i = 0; i < bits; i++)
	{
		number /= 10;
	}
	return number + integerPart;
}