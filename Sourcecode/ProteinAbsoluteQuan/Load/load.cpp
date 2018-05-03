/*
#  Copyright (C) 2015-2018 all rights reserved

#  This program is a free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/Licenses/
*/

// load.cpp : Define the entry point of the console application
#include "stdafx.h"
#include"DataIO.h"
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc <= 1)
	{
		cout << "ERROR: Usage: Load [parameters file]" << endl;
		cout << "Input parameters file\'s Format:" << endl;
		cout << "xxx.params" << endl;
		cout << "Version -0.1" << endl;
		system("pause");
		return -1;
	}
	clock_t begin, end;
	begin = clock();

	CLoadParam param;

	string ParaFilePath=Unicode2Multibyte(argv[1]);
	if (!CheckFilePath(ParaFilePath))
	{
		return 2;
	}

	string strResultPath = GetResultPath(ParaFilePath);
	if(!flog.mf_Init(strResultPath))
	{
		cout << "Cannot creat log file in " << strResultPath<< ". " << endl;
		return 3;
	}
	string strPrologue = AddDecorateStar("The beginning of Load module");
	cout << strPrologue;
	flog.mf_Input(strPrologue);
	param.SetParameters(ParaFilePath);
	param.CheckParameters();
    if (param.m_eDataType == MaxQuantTpye)
	{	//read maxquant results
		CLoadMaxQuantIO Proteinsio(param);
		Proteinsio.mf_LoadProteinIDAndSequence(param.m_strProteinFasterFilePath, param.m_eFastaType);
		Proteinsio.mf_LoadPeptidesWithoutRedunPeptides( param);
		Proteinsio.mf_LoadProteins( param);
		Proteinsio.mf_saveProteins( param);
	}
	else if (param.m_eDataType == mzQuantMLType)
	{
		CloadmzQuantMLIO Proteinsio;
		Proteinsio.mf_LoadProteinIDAndSequence(param.m_strProteinFasterFilePath, param.m_eFastaType);
		Proteinsio.load(param);
		Proteinsio.mf_saveProteins( param);
	}
	else if (param.m_eDataType == PeakViewType)
	{
		CLoadPeakViewIO PeakViewio;
		PeakViewio.load(param);
		PeakViewio.mf_saveProteins(param);
	}
	
	end = clock();
	int iRuntime = (end - begin) / CLOCKS_PER_SEC;
	cout << "The Load module took " << iRuntime << " seconds." << endl;
	string strRuntime = fInt2String(iRuntime);
	flog.mf_Input("The Load module took " +strRuntime+" seconds.\n");
	string strEnd = AddDecorateStar("The end of Load module");
	cout << strEnd;
	flog.mf_Input(strEnd);
	flog.mf_Destroy();
	return 0;
}



