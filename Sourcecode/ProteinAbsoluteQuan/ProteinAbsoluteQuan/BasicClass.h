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
#include"PepAttributeWorker.h"
#include<sys/stat.h>
#include<windows.h>
#include<regex>

using namespace std;


typedef std::regex FastaType;

extern int UniquePeptidesTrainThreshold;
const int FlankingRegionWidth = 15;
const int NumberOfFeature = 587;

class CProtein
{
public:
	CProtein();
	void Clear();
	string m_strProteinID;
	string m_strProteinFullName;
	vector<int> m_vPeptidesIDs;
	vector<string> m_vPeptidesSequencesWithFlankingRegion;
	vector<string> m_vPeptidesSequences;
	vector<double> m_vdPeptidesMW;
	map<string, vector<string>> m_mapPeptidesSequenceAndCorrespondRawFiles;
	map < string, vector<double>> m_mapExperimentAndPeptidesIntensity;
	vector<CPeptideAttribute> m_vPeptidesFeatures;
	vector<double> m_vPeptideExperimentalQFactors; 
	vector<double> m_vPeptidesNativeIntensity; 
	vector<double> m_vPeptidesCorrectedIntensity;
	vector<double> m_vPeptideQfactors;
	vector<double> m_vecdMaxQuantiBAQ;
	vector<double> m_vPeptidesCombinedIntensity;

	double m_dPeptidesIntensityMedian;
	double m_dPeptidesIntensityMax;
	double m_dProteinLFAQpep;
	double m_dProteinLFAQpro;
	double m_dProteinQfactor;
	double m_dMaxQuantiBAQ;
	double m_dReCalculateiBAQ;
	bool m_bIfInTrainSet;  // if use it in the training set

	double m_dPredictedMolOfLFAQpep;
	double m_dPredictedMolOfLFAQpro;
	double m_dPredictedMolOfiBAQ;
	double m_dPredictedMolOfTop3;

	double m_dCVNative;
	double m_dCVAfterCorrected;
	double m_dPeptidesIntensitySum;
	int m_iPeptidesNumber;
	int m_iNumberOfTheoreticEnzyme;
	bool m_bIfCalculateLFAQpep;
	bool m_bIfCalculateLFAQpro;
	double m_dProteinIntensityTop3;
	int m_iPeptidesMaxIntensityIndex;
	
	string m_strProteinSequence;
	string mf_GetPeptidesAdjacentSequence(string PeptideSequence);
};
class Cpeptides
{
public:
	void mf_Clear();

	map<int, vector<double>> m_mapPeptideIDAndIntensitys;
	map<int, double> m_mapPeptideIDAndCombinedIntensity;

	map<int, string> m_mapPeptideIDAndSequence;
	map<int, double> m_mapPeptideIDAndIntensity;
	map<string, double> m_mapPeptideSequenceAndIntensity;
	map<string, bool> m_mapPeptideSequencuAndIfUse;
	map<string, CPeptideAttribute> m_mapPeptideSequenceAndAttribute;
	string m_strPeptideIntensityName;
};

class CRelatedProtein
{
public:
	CRelatedProtein()
	{
		bIfSubSet = false;
		strProteinID = "";
		vecPeptides.clear();
		vecSubsetProteinIDs.clear();
		vecBIfPeptideGroupUnique.clear();
	}
	string strProteinID;
	bool bIfSubSet;  // If true, this protein is the subset of some protein
	vector<string> vecPeptides;
	vector<string> vecSubsetProteinIDs; // include equivalence set
	vector<bool> vecBIfPeptideGroupUnique;
	// judge if RelatedProtein A is the subset of this RelatedProtein;
	bool mf_ifSubset(CRelatedProtein& A); 
	void Clear();
};

class CMergedProtein
{
public:
	void Clear();
	string m_strPrtoteinName;
	string m_strSubsetProteins;
	vector<string> m_vecExperiments;
	vector<double> m_vecMaxQuantiBAQOfExperiments;
	vector<double> m_vecReCaliBAQOfExperiments;
	vector<double> m_vecLFAQpepOfExperiments;
	vector<double> m_vecTop3OfExperiments;
	vector<double> m_vecLFAQproOfExperiments;
	vector<double> m_vecPredictedMolsByLFAQpepOfExperiments;
	vector<double> m_vecPredictedMolsByLFAQproOfExperiments;
	vector<double> m_vecPredictedMolsByiBAQOfExperiments;
	vector<double> m_vecPredictedMolsByTop3OfExperiments;
	map<string,vector <double>> m_mapExperimentsAndPeptidesIntensity;
};

class CMergedProteins
{
public:
	void mf_ExperimentsAnalysis();
	vector<CMergedProtein> m_vecMergedProteins;
	vector<string> m_strExpreriments;
	vector<double> m_vecMaxQuantiBAQCV;
	vector<double> m_vecReCaliBAQCV;
	vector<double> m_vecLFAQpepCV;
	vector<double> m_vecLFAQproCV;
	vector<double> m_vecTop3CV;
	
};
string AddDecorateStar(const string& strProjectName);

string Unicode2Multibyte(_TCHAR *chars);
map<string, string> GetParametersFromFile(const string &ParamFilePath);
string GetResultPath(string ParaFilePath);
bool CheckFilePath(string& filepath);
 
// about string
bool bPrefix(string str, string prefix);
string stringTrim(string str);
string removestringWhitespaces(string str);
string stringtoLower(string str);
string PrefixOfstring(string str, char delim);
string SuffixOfstring(string str, char delim);
bool stringSplit(string str, const string splitter, vector<string>& substrings);

double CalculatelogCV(const vector<double>&s);
double Average(double *Array, int Number);
double Average(vector<double> Array);
void Discretization(double *TestY, int m_iTestSampleNumber,int iBinNum=10);
void Discretization(vector<double> &TestY, int iBinNum=10);
double GetMaxValue(const vector<double> PeptidesIntensity);
double GetMaxValue(vector<double> PeptidesIntensity);
double GetMaxValue(double *PeptidesIntensity,int iNum);
double GetMinValue(vector<double> PeptidesIntensity);
double GetMinValue(double *PeptidesIntensity, int iNum);
double GetSum(vector<double> vec);

int GetMaxIndex(vector<double> vec);
double GetMedianIntensity(vector<double> PeptidesIntensity);
double IntegralFromVector(vector<double> Values, vector<double> times);

void GetQuantiles(vector<double> PeptidesIntensity,double &FirstQ,double &median,double &ThirdQ);

void Assertion(vector<double> &TestY);
void Assertion(double *TestY, int m_iTestSampleNumber);
void Scaling(vector<double> &TestY);
void Scaling(double *TestY, int m_iTestSampleNumber);
// Descending
void quickSort(double s[], int index[], int l, int r);
void quickSort(vector<size_t> &s, vector<size_t> &index, int l, int r);
void quickSort(vector<double> &s, vector<size_t> &index, int l, int r);

double spearsonCorrelation(double s1[], double s2[],int numberOfSample);
double spearsonCorrelation(vector<double> s1, double s2[]);
double spearsonCorrelation(vector<double> s1, vector<double> s2);
string strim(string &str);
vector<string> split(const string &str, const string& splitter);

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns, string sep="\t");
bool fStringToBool(string str,bool &bl);
string fInt2String(int i);
string fSize_t2String(size_t si);
string fDouble2String(double d);
bool is_file(const char* FILENAME);
bool is_dir(const char* path);

void SaveVector(string path, vector<double> vec);
void SaveVector(string path, vector<string> vec);
void SaveArray(string path, double * arr, int len);

// Judge if vector A is the subset of vector B;
bool ifSubset(vector<string>&A, vector<string>&B);
bool MoreThan(double a, double b);

double round(double number, unsigned int bits);