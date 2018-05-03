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

#include"stdafx.h"
#include"loadParam.h"


CLoadParam::CLoadParam()
{
}
/*
load the parameters
*/
void  CLoadParam::SetParameters(string parafilepath)
{
	cout << "Setting Parameters from " << parafilepath<<"\n";
	flog.mf_Input("Setting Parameters from " + parafilepath+"\n");

	map <string, string>  mapArttributesInFile;
	map<string, string>::iterator mapAttributesInfileIter;
	mapArttributesInFile = GetParametersFromFile(parafilepath);
    string strFastaTypetemp;	
	string strLoadFileType;

	mapAttributesInfileIter = mapArttributesInFile.find("identificationfiletype");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		strLoadFileType = mapAttributesInfileIter->second;
		if (strLoadFileType == "maxquant")
		{
			m_eDataType = MaxQuantTpye;
		}
		else if (strLoadFileType == "mzQuantML")
		{
			m_eDataType = mzQuantMLType;
		}
		else if (strLoadFileType == "PeakView")
		{
			m_eDataType = PeakViewType;
		}
		else
		{
			cout << "Error:\tThe format of parameter identificationfiletype is not correct!" << endl;
			flog.mf_Input("Error:\tCannot support the identification result of " + strLoadFileType + "\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	else
	{
		cout << "Error:\tThe format of parameter identificationfiletype is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter identificationfiletype is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}

	if (m_eDataType == MaxQuantTpye)
	{
		mapAttributesInfileIter = mapArttributesInFile.find("input directory path");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			m_strIdentResultPath = mapAttributesInfileIter->second;
			if (m_eDataType == MaxQuantTpye)
			{
				m_strPeptideFilePath = m_strIdentResultPath + "\\peptides.txt";
				m_strProteinFilePath = m_strIdentResultPath + "\\proteinGroups.txt";
				m_strExprimentDesignPath = m_strIdentResultPath + "\\experimentalDesignTemplate.txt";
			}
		}
		else
		{
			cout << "Error:\tThe format of parameter identification  result file path is not correct." << endl;
			flog.mf_Input("Error:\tThe format of parameter identification  result file path is not correct.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	if (m_eDataType == PeakViewType || m_eDataType == mzQuantMLType)
	{
		mapAttributesInfileIter = mapArttributesInFile.find("input file path");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			m_strIdentResultPath = mapAttributesInfileIter->second;
		}
		else
		{
			cout << "Error:\tThe format of parameter \"Input file path\" is not correct." << endl;
			flog.mf_Input("Error:\tThe format of parameter \"Input file path\" is not correct.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	
	mapAttributesInfileIter = mapArttributesInFile.find("fastapath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strProteinFasterFilePath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Error:\tThe format of parameter fastapath is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter fastapath is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("identifierparsingrule");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			strFastaTypetemp = ">(.*?)\\s";
		}
		else
		{
			strFastaTypetemp = mapAttributesInfileIter->second;
		}
		m_eFastaType = strFastaTypetemp;
	}
	else
	{
		cout << "Error:\tThe format of parameter IdentifierParsingRule is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter IdentifierParsingRule is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttributesInfileIter = mapArttributesInFile.find("ifexistdecoyproteins");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_bIfExistDecoyProtein = false;
		}
		else
		{
			if (!fStringToBool(mapAttributesInfileIter->second, m_bIfExistDecoyProtein))
			{
				cout << "Error:\tCannot convert " << mapAttributesInfileIter->second << " to bool!" << endl;
				flog.mf_Input("Error:\tCannot convert " + mapAttributesInfileIter->second + " to bool!\n");
				flog.mf_Destroy();
				exit(1);
			}
		}		
	}
	else
	{
		cout << "Error:\tThe format of parameter ifexistdecoyproteins is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter ifexistdecoyproteins is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("prefixofdecoyprotein");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strDecoyProteinIDPrefix = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Error:\tThe format of parameter prefixofdecoyprotein is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter prefixofdecoyprotein is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}
	
	mapAttributesInfileIter = mapArttributesInFile.find("ifexistcontaminantproteins");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_bIfExistContanProtein = false;
		}
		else
		{
			if (!fStringToBool(mapAttributesInfileIter->second, m_bIfExistContanProtein))
			{
				cout << "Cannot convert " << mapAttributesInfileIter->second << " to bool!" << endl;
				flog.mf_Input("Error:\tCannot convert " + mapAttributesInfileIter->second + " to bool!\n");
				flog.mf_Destroy();
			}

		}
	}
	else
	{
		cout << "Error:\tThe format of parameter ifexistcontaminantproteins is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter ifexistcontaminantproteins is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttributesInfileIter = mapArttributesInFile.find("prefixofcontaminantprotein");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strContaminantfix = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Error:\tThe format of parameter prefixofcontaminantprotein is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter prefixofcontaminantprotein is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("resultpath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strQuantiResultPath = mapAttributesInfileIter->second;
		m_strExtractedProteinsPath = m_strQuantiResultPath + "\\ProteinsInfo.txt";
	}
	else
	{
		cout << "Error:\tThe format of parameter resultpath is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter resultpath is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}
	flog.mf_Input("Done\n");
}

void CLoadParam::CheckParameters()
{

	if (m_eDataType == MaxQuantTpye)
	{
		if (is_dir(m_strIdentResultPath.c_str())==-1)
		{
			cout << "Error:\tThe identification  result file path should be a directory!" << endl;
			flog.mf_Input("Error:\tThe identification  result file path should be a directory!\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	else if (m_eDataType == mzQuantMLType)
	{
		if (!is_file(m_strIdentResultPath.c_str()))
		{
			cout << "Error:\tThe identification  result file path should be a file!" << endl;
			flog.mf_Input("Error:\tThe identification  result file path should be a file!\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	if (!is_file(m_strProteinFasterFilePath.c_str()))
	{
		cout << "Error:\tThe path of protein fasta is not correct!" << endl;
		flog.mf_Input("Error:\tThe path of protein fasta is not correct!\n");
		flog.mf_Destroy();
		exit(1);
	}
	if (m_bIfExistDecoyProtein)
	{
		if (m_strDecoyProteinIDPrefix == "")
		{
			cout << "Error:\tPlease check the setting of the prefix of decoy protein!" << endl;
			flog.mf_Input("Error:\tPlease check the setting of the prefix of decoy protein!\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	if (m_bIfExistContanProtein)
	{
		if (m_strContaminantfix == "")
		{
			cout << "Error:\tPlease check the setting of the prefix of contaminant protein!" << endl;
			flog.mf_Input("Error:\tPlease check the setting of the prefix of contaminant protein!\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	if (is_dir(m_strQuantiResultPath.c_str())==-1)
	{
		cout << "Error:\tThe result file path must be a directory!";
		flog.mf_Input("Error:\tThe result file path must be a directory!\n");
		flog.mf_Destroy();
		exit(1);
	}

}

void CLoadParam::ShowParam()
{
	cout << "Peptides come from " << m_strPeptideFilePath << endl;
	cout << "Proteins come from " << m_strProteinFilePath << endl;
	cout << "And the corresponding fasta file is " << m_strProteinFasterFilePath << endl;
	cout << "We need Expriment Design File " << m_strExprimentDesignPath << endl;
}