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

#include"BARTRegression.h"
#include"stepwise.h"
// this class is in charge of regression to correct peptides intensity
class CRegression
{
public:
	void mf_RegressionRun(const vector<CProtein> &vTrainProteins, vector<CProtein> & vTestProteins, CQuantificationParam trainparam, string strExperimentName);
	void mf_StepwiseRegression(vector<CProtein> &proteins, CQuantificationParam params, string strExperimentName);
	void mf_BARTRegression(const vector<CProtein> &vTrainProteins, vector<CProtein> & vTestProteins, CQuantificationParam trainparam, string strExperimentName);
	

};
