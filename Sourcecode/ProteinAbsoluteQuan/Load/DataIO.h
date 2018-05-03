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

#pragma once
#include"mzQuantMLHandler.h"
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include<iomanip>
using namespace std;

class CLoadIO
{
public:
	vector<CProtein> m_vecProteins;
	Cpeptides m_Peptides;

	int m_NumberOfPeptides;
	int m_iNumberOfExperiments;
	map<string, string> m_mapProteinIDAndSequence;
	map<string, int> m_mapPeptideSequenceAndAppearNumber;
	map<string, string> m_mapExperimentNameAndPeptideIntensityName;

	string mf_GetPeptidesSequenceFromID( int ID);
	double mf_GetPeptidesIntensityFromID( int ID, int iExperimentIndex);
	double mf_GetCombinedIntensityFromID(int ID);
	double mf_GetMedianIntensity( vector<double> PeptidesIntensity);
	void mf_saveProteins( CLoadParam param);
	bool mf_LoadProteinIDAndSequence( string strProteinFasterFilePath, FastaType fastatype);
	
	bool mf_LoadPeptides(string strPeptideFilePath);
	bool mf_LoadProteins( string strProteinFilePath, string strProteinFasterFilePath);
	//void mf_GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns,string sep="\t");

};

class CLoadMaxQuantIO :public CLoadIO
{
public:

	void mf_GetAttributesName(string ExperimentDesignPath);
	CLoadMaxQuantIO( const CLoadParam loadParam);
	
	//The columns of iBAQ_intrensity and Peptides_intensity are not fixed£¬so we get them based on the input file every time.
	int mf_GetPeptideAttributecolumn( string AttributeName);
	int mf_GetProteinAttributecolumn( map<string, int>mapAttrtibuteAndcolumns,
									int iDecoycolumnNum, int iPoteinalContaminantcolumnNum, 
									int iPeptideIDscolumnNum, int  iPeptideRazorcolumnNum, 
									vector<int> veciBAQcolumnNum);  //analysis the first row to get the attribute columns.
	
	map<string, string> m_mapExperimentNameAndiBAQIntensityName;

	bool mf_LoadPeptides( CLoadParam param);
	bool mf_LoadProteins( CLoadParam& param);
	bool mf_LoadPeptidesWithoutRedunPeptides( CLoadParam param);


};

class CloadmzQuantMLIO :public CLoadIO
{
public:
	///Default constructor
	CloadmzQuantMLIO();
	///Destructor
	virtual ~CloadmzQuantMLIO();

	/**
	@brief Loads a map from a mzQuantML file.
	@exception Exception::FileNotFound is thrown if the file could not be opened
	@exception Exception::ParseError is thrown if an error occurs during parsing
	*/
	void load( const CLoadParam &param);
	void mf_saveProteins( CLoadParam param);

protected:
	void parse_(const string & filename, XMLHandler * handler);

	// Remoce the Decoy protein;
	void mf_RemovDecoyProtein(CQuantInformation &m_qi, const CLoadParam &param);
	//Remove shared peptides;
	void mf_RemovedSharedPeptides(CQuantInformation &m_qi);

	void mf_CheckMatchedFasta( const CLoadParam &param);
	void mf_CheckEnoughInformation();
	CQuantInformation m_qi_;
};

class PeptideQuantInfo
{
public:
	string strPeptideSequence;
	string strProteinID;
	string strProteinGroups;
	vector<string> vecPrecursorCharges;
	vector<double> vecdPeptideIntensities; //MS1 intensity or MS2 intensity
	vector<int> vecRepNumOfPrecursorCharge;
	string strCondition;
	void Clear()
	{
		strProteinID = "";
		strProteinGroups = "";
		strCondition = "";
		vecPrecursorCharges.clear();
		vecdPeptideIntensities.clear();
		vecRepNumOfPrecursorCharge.clear();
	}
};

class CLoadPeakViewIO :public CLoadIO
{
public:
	void load(const CLoadParam &param);
	void mf_SaveProteinGroups(string path, const map<string, CRelatedProtein>& mapDistinctProIDAndRelatedInfo);
	void mf_SaveSharedPeptides(string path, const map<string, vector<string>>& mapSharedPeptidesAndProteins);
	void ScreenSharedPeptides(const CLoadParam &param, map<string, vector<string>>& mapSharedPeptidesAndProteins,
		map<string, string> &mapGUniquePeptideAndProteinGroups, map<string, vector<string>>& mapProteinGroupAndGUniquePeptides);
	void mf_saveProteins(CLoadParam param);

private:
	void mf_GetPeptidesAndProteins(const CLoadParam &param, map<string, vector<string>>& mapPeptidesAndProteinNames, map<string, vector<string>>& mapProteinAndPeptides);
	bool mf_GetPeakViewAttributeColumns(const map<string, int> &mapAttrtibuteAndcolumns, int &iPeptidecolumnNum, int &iProteinsNamecolumnNum,
		vector<int> &vecIntensitiesColumn,int &iPrecursorChargeColumnNum);
};





