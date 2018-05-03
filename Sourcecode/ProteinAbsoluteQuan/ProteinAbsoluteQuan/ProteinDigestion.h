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
#include"QuantificationParam.h"
using namespace std;
class CQuantificationParam;

class CProteinDigestion
{
public:
	CProteinDigestion(CQuantificationParam param);

	/*-----Get Method----------------*/
	vector <string>  getAllPossiblePeptides(string ProteinSequence); //so get all possible peptides without consider missing cut
	vector <string>  getAllPossiblePeptidesOpt(string ProteinSequence);//optimized method
	bool isInCutPosition(char aChar, string cutPosition);

	int mf_CalculateOptDigestionNumber(string ProteinSequence);
	int mf_CalculateDigestionNumberWithoutMissingCut(string ProteinSeqence);
private:
	string m_strCutPosition;//the position of the protein where enzyme works
	bool m_blr;//cut to the left position or to the right, T-->left,  False-->right
	int m_imaxMissedClevage;//the maxium missed clevage allowed
	int m_iPeptideCount;//total peptide that might be extracted from proteinSequence
	int m_iMinPepLength;//the required minimum length of the peptides, which is used to select peptides, known as L1 in [L1, L2]
	int m_iMaxPepLength;//the required maxium length of the peptides, which is used to select peptides, known as L2 in [L1, L2]

};

