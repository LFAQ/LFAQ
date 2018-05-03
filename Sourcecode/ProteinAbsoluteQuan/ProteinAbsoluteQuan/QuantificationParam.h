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

//File version: 2015-12-04
/*this class contain all parameters which used in the stage of training and predicting
*/
#pragma once
#include "BasicClass.h"

using namespace std;
extern int UniquePeptidesTrainThreshold;


class CQuantificationParam
{
public:
	string m_strProteinsPath;
	string m_strRegressionResult;
	string m_strIdentifySoftwareType;
	string m_strProteinFastaPath;
	FastaType m_eFastaType;
	string m_strIdentResultPath;
	string m_strQuantiResultPath;
	string m_strExperimentDesignPath;
	bool m_bIfCotainStandardProtein;
	string m_strIdentifierOfStandPro;
	string m_strStandardProteinsFilePath;

	string m_strProteinLFAQpepPath;
	string m_strCutPosition;        //the position of the protein where enzyme works
	bool m_blr;                     //cut to the left position or to the right, T-->left,  False-->right
	int m_imaxMissedClevage;        //the maxium missed clevage allowed
	int m_iMinPepLength;            //the required minimum length of the peptides, which is used to select peptides, known as L1 in [L1, L2]
	int m_iMaxPepLength;            //the required maxium length of the peptides, which is used to select peptides, known as L2 in [L1, L2]
	bool m_bIfCalculateiBAQ;        // return iBAQ from MaxQuant or iBAQ Calculated by the program;
	bool m_bIfCalculateTop3;

	double dDistributionCorrectFactor;
	double dMeanQfactor;
	double dStdQfactor;
	string strMaxPeptideQfactorType;
	string strDistributionLenghType;

	bool m_bIfExistDecoyProtein;
	bool m_bIfExistContanProtein;
	string m_strDecoyProteinIDPrefix;  //The prefix of decoy protein ID;
	string m_strContaminantfix; //The prefix of Contaminant protein ID;
	// for stepwise
	double m_dAlpha1;
	double m_dAlpha2;
	// for BART 
	double m_dAlpha;
	double m_dBeta;
	double m_iNumberOfTrees;
	int m_ik;
	string m_strRegressionMethod;
	bool m_bIfOptimizeParameters;

	void mf_setparameters(string parafilepath);
	int mf_CheckParameters();
};