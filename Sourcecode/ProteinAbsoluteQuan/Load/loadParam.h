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
#include"BasicClass.h"
class CLoadParam
{
public:
	CLoadParam();
	void SetParameters(string parafilepath);
	void CheckParameters();
	void ShowParam();

	DataType m_eDataType;
	FastaType m_eFastaType;

	string m_strIdentResultPath;
	string m_strPeptideFilePath;
	string m_strProteinFilePath;
	string m_strProteinFasterFilePath;
	string m_strExprimentDesignPath;
	string m_strQuantiResultPath;
	string m_strExtractedProteinsPath;

	bool m_bIfiBAQIntensityExist;
	bool m_bIfExistDecoyProtein;
	bool m_bIfExistContanProtein;
	
	string m_strDecoyProteinIDPrefix;  //The prefix of decoy protein ID;
	string m_strContaminantfix; //The prefix of Contaminant protein ID;
	string m_strGeneProteinsPath;
	string m_strGeneExpressionQuantificationResultPath;

	double m_dIsotopeDotProductThreshold;
	double m_dLibDotProductThreshold;
};
