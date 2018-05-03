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

#include"stdafx.h"
#include"QuantificationParam.h"


/*parameters
string m_strProteinsPath;
string m_strPeptidesInfoPath;
double m_dBestAlpha1;
double m_dBestAlpha2;
*/
int CQuantificationParam::mf_CheckParameters()
{
	if (!is_dir(m_strIdentResultPath.c_str()))
	{
		cout << "Error:\tThe \"Input directory path\" must be a directory!";
		flog.mf_Input("Error:\tThe \"Input directory path\" must be a directory!\n");
		flog.mf_Destroy();
		exit(1);
	}

	if (!is_file(m_strProteinFastaPath.c_str()))
	{
		cout << "Error:\tThe path of protein fasta is not correct!" << endl;
		flog.mf_Input("Error:\tThe path of protein fasta is not correct!\n");
		flog.mf_Destroy();
		exit(1);

	}

	if (!is_dir(m_strQuantiResultPath.c_str()))
	{
		cout << "Error:\tThe result file path must be a directory!";
		flog.mf_Input("Error:\tThe result file path must be a directory!\n");
		flog.mf_Destroy();
		exit(1);
	}


}

void CQuantificationParam::mf_setparameters(string parafilepath)
{
	cout << "Setting parameters\n";
	flog.mf_Input("Setting parameters\n");

	map <string, string>  mapArttributesInFile;
	map<string, string>::iterator mapAttributesInfileIter;
	mapArttributesInFile = GetParametersFromFile(parafilepath);

	mapAttributesInfileIter = mapArttributesInFile.find("identificationfiletype");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strIdentifySoftwareType = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Error:\tParameter \"IdentificationFileType\" is not correct" << endl;
		flog.mf_Input("Error:\tParameter \"IdentificationFileType\" is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("fastapath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strProteinFastaPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Error:\tParameter \"Fastapath\" is not correct" << endl;
		flog.mf_Input("Error:\tParameter \"Fastapath\" is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	string strFastaTypetemp;
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
		cout << "Error:\tParameter \"IdentifierParsingRule\" is not correct" << endl;
		flog.mf_Input("Error:\tParameter \"IdentifierParsingRule\" is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("resultpath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strQuantiResultPath = mapAttributesInfileIter->second;
		if (m_strIdentifySoftwareType == "mzQuantML")
		{
			m_strExperimentDesignPath = m_strQuantiResultPath+"\\mzQuantMLExperiment.config";
		}
		else if (m_strIdentifySoftwareType == "PeakView")
		{
			m_strExperimentDesignPath = m_strQuantiResultPath + "\\PeakViewExperiment.config";
		}
		else
		{
			mapAttributesInfileIter = mapArttributesInFile.find("input directory path");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				m_strIdentResultPath = mapAttributesInfileIter->second;
				if (m_strIdentifySoftwareType == "maxquant")
				{
					m_strExperimentDesignPath = m_strIdentResultPath + "\\experimentalDesignTemplate.txt";
				}
			}
			else
			{
				cout << "Error:\tParameter \"Input directory path\" is not correct" << endl;
				flog.mf_Input("Error:\tParameter \"Input directory path\" is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}
		}
	}
	else
	{
		cout << "Error:\tParameter \"ResultPath\" is not correct" << endl;
		flog.mf_Input("Error:\tParameter \"ResultPath\" is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("regressionmethod");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_strRegressionMethod = "BART";
		}
		else
		{
			m_strRegressionMethod = mapAttributesInfileIter->second;
		}		
		if (m_strRegressionMethod == "BART")
		{
			mapAttributesInfileIter = mapArttributesInFile.find("alpha");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				if (mapAttributesInfileIter->second == "")
				{
					m_dAlpha = 0.85;
				}
				else
				{
					m_dAlpha = atof(mapAttributesInfileIter->second.c_str());
				}				
			}
			else
			{
				cout << "Error:\tParameter \"alpha\" is not correct" << endl;
				flog.mf_Input("Error:\tParameter \"alpha\" is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}

			mapAttributesInfileIter = mapArttributesInFile.find("beta");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				if (mapAttributesInfileIter->second == "")
				{
					m_dBeta = 1.6;
				}
				else
				{
					m_dBeta = atof(mapAttributesInfileIter->second.c_str());
				}
			}
			else
			{
				cout << "Error:\tParameter \"beta\" is not correct" << endl;
				flog.mf_Input("Error:\tParameter beta is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}

			mapAttributesInfileIter = mapArttributesInFile.find("k");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				if (mapAttributesInfileIter->second == "")
				{
					m_ik = 2;
				}
				else
				{
					m_ik = atoi(mapAttributesInfileIter->second.c_str());
				}
				
			}
			else
			{
				cout << "Error:\tParameter k is not correct" << endl;
				flog.mf_Input("Error:\tParameter k is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}

			mapAttributesInfileIter = mapArttributesInFile.find("number of trees");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				if (mapAttributesInfileIter->second == "")
				{
					m_iNumberOfTrees = 200;
				}
				else
				{
					m_iNumberOfTrees = atoi(mapAttributesInfileIter->second.c_str());
				}
				
			}
			else
			{
				cout << "Error:\tParameter Number of trees is not correct" << endl;
				flog.mf_Input("Error:\tParameter Number of trees is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}
		}
		else if (m_strRegressionMethod == "stepwise")
		{
			mapAttributesInfileIter = mapArttributesInFile.find("alpha1");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				if (mapAttributesInfileIter->second == "")
				{
					m_dAlpha1 = 0.95;
				}
				else
				{
					m_dAlpha1 = atof(mapAttributesInfileIter->second.c_str());
				}				
			}
			else
			{
				cout << "Error:\tParameter alpha1 is not correct" << endl;
				flog.mf_Input("Error:\tParameter alpha1 is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}

			mapAttributesInfileIter = mapArttributesInFile.find("alpha2");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				if (mapAttributesInfileIter->second == "")
				{
					m_dAlpha2 = 0.95;
				}
				else
				{
					m_dAlpha2 = atof(mapAttributesInfileIter->second.c_str());
				}				
			}
			else
			{
				cout << "Error:\tParameter Alpha2 is not correct" << endl;
				flog.mf_Input("Error:\tParameter Alpha2 is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}
		}
		else
		{
			cout << "Error:\tParameter BART is not correct" << endl;
			flog.mf_Input("Error:\tParameter BART is not correct\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	else
	{
		cout << "Error:\tParameter RegressionMethod is not correct" << endl;
		flog.mf_Input("Error:\tParameter RegressionMethod is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	m_bIfOptimizeParameters = false;

	mapAttributesInfileIter = mapArttributesInFile.find("ifcotainstandardprotein");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_bIfCotainStandardProtein = false;
		}
		else
		{
			if (!fStringToBool(mapAttributesInfileIter->second, m_bIfCotainStandardProtein))
			{
				cout << "Error:\tCannot convert " << mapAttributesInfileIter->second << " to bool!" << endl;
				flog.mf_Input("Error:\tCannot convert " + mapAttributesInfileIter->second + " to bool!\n");
				flog.mf_Destroy();
				exit(1);
			}
		}

		if (m_bIfCotainStandardProtein)
		{
			mapAttributesInfileIter = mapArttributesInFile.find("identifierofstandardprotein");
			if (mapAttributesInfileIter != mapArttributesInFile.end())
			{
				m_strIdentifierOfStandPro = mapAttributesInfileIter->second;
			}
			else
			{
				cout << "Error:\tParameter IdentifierOfStandardProtein is not correct" << endl;
				flog.mf_Input("Error:\tParameter IdentifierOfStandardProtein is not correct\n");
				flog.mf_Destroy();
				exit(1);
			}
		}
	}
	else
	{
		cout << "Error:\tParameter IfCotainStandardProtein is not correct" << endl;
		flog.mf_Input("Error:\tParameter IfCotainStandardProtein is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("standardproteinsfilepath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_strStandardProteinsFilePath = "StandardProteins.txt";
		}
		else
		{
			m_strStandardProteinsFilePath = mapAttributesInfileIter->second;
		}		
	}
	else
	{
		cout << "Error:\tParameter StandardProteinsFilePath is not correct" << endl;
		flog.mf_Input("Error:\tParameter StandardProteinsFilePath is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("uniquepeptidestrainthreshold");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		UniquePeptidesTrainThreshold = atoi(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		UniquePeptidesTrainThreshold = 5;
	}

	if (m_strIdentifySoftwareType == "mzQuantML" || m_strIdentifySoftwareType == "PeakView")
	{
		mapAttributesInfileIter = mapArttributesInFile.find("input file path");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			m_strIdentResultPath = mapAttributesInfileIter->second;
		}
		else
		{
			cout << "Error:\tParameter \"Input file path\" is not correct" << endl;
			flog.mf_Input("Error:\tParameter \"Input file path\" is not correct\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	else
	{
		mapAttributesInfileIter = mapArttributesInFile.find("input directory path");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			m_strIdentResultPath = mapAttributesInfileIter->second;
		}
		else
		{
			cout << "Error:\tParameter \"Input directory path\" is not correct" << endl;
			flog.mf_Input("Error:\tParameter \"Input directory path\" is not correct\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	mapAttributesInfileIter = mapArttributesInFile.find("maxmissedcleavage");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_imaxMissedClevage = 1;
		}
		else
		{
			m_imaxMissedClevage = atoi(mapAttributesInfileIter->second.c_str());
		}		
	}
	else
	{
		cout << "Error:\tParameter MaxMissedCleavage is not correct" << endl;
		flog.mf_Input("Error:\tParameter MaxMissedCleavage is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("pepshortestlen");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_iMinPepLength = 6;
		}
		else
		{
			m_iMinPepLength = atoi(mapAttributesInfileIter->second.c_str());
		}		
	}
	else
	{
		cout << "Error:\tParameter PepShortestLen is not correct" << endl;
		flog.mf_Input("Error:\tParameter PepShortestLen is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("peplongestlen");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_iMaxPepLength = 38;
		}
		else
		{
			m_iMaxPepLength = atoi(mapAttributesInfileIter->second.c_str());
		}	
	}
	else
	{
		cout << "Error:\tParameter PepLongestLen is not correct" << endl;
		flog.mf_Input("Error:\tParameter PepLongestLen is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttributesInfileIter = mapArttributesInFile.find("ifcalculateibaq");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_bIfCalculateiBAQ = true;
		}
		else
		{
			if (!fStringToBool(mapAttributesInfileIter->second, m_bIfCalculateiBAQ))
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
		cout << "Error:\tParameter IfCalculateiBAQ is not correct" << endl;
		flog.mf_Input("Error:\tParameter IfCalculateiBAQ is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("ifcalculatetop3");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			m_bIfCalculateTop3 = true;
		}
		else
		{
			if (!fStringToBool(mapAttributesInfileIter->second, m_bIfCalculateTop3))
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
		cout << "Error:\tParameter IfCalculateTop3 is not correct" << endl;
		flog.mf_Input("Error:\tParameter IfCalculateTop3 is not correct\n");
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

	mapAttributesInfileIter = mapArttributesInFile.find("distributioncorrectfactor");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			dDistributionCorrectFactor = 10.0;
		}
		else
		{
			dDistributionCorrectFactor = atof(mapAttributesInfileIter->second.c_str());
		}		
	}
	else
	{
		dDistributionCorrectFactor = 10;
		cout << "Warning:\tParameter DistributionCorrectFactor is set 10 by default" << endl;
		flog.mf_Input("Waring:\tParameter DistributionCorrectFactor is set 10 by default\n");
	}

	mapAttributesInfileIter = mapArttributesInFile.find("maxpeptideqfactortype");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (mapAttributesInfileIter->second == "")
		{
			strMaxPeptideQfactorType = "ExponDistExpect";
		}
		else
		{
			strMaxPeptideQfactorType = mapAttributesInfileIter->second;
		}
	}
	else
	{
		strMaxPeptideQfactorType = "ExponDistExpect";
	}

	if (strMaxPeptideQfactorType == "ExponDistMonteCarlo" || strMaxPeptideQfactorType == "ExtremDistMonteCarlo")
	{
		mapAttributesInfileIter = mapArttributesInFile.find("meanqfactor");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			dMeanQfactor = atof(mapAttributesInfileIter->second.c_str());
		}
		else
		{
			cout << "Error:\tParameter MeanQfactor is not correct" << endl;
			flog.mf_Input("Error:\tParameter MeanQfactor is not correct\n");
			flog.mf_Destroy();
			exit(1);
		}
		mapAttributesInfileIter = mapArttributesInFile.find("stdqfactor");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			dStdQfactor = atof(mapAttributesInfileIter->second.c_str());
		}
		else
		{
			cout << "Parameter StdQfactor is not correct" << endl;
			flog.mf_Input("Error:\tParameter StdQfactor is not correct\n");
			flog.mf_Destroy();
			exit(1);
		}


		mapAttributesInfileIter = mapArttributesInFile.find("distributionlenghtype");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			strDistributionLenghType = mapAttributesInfileIter->second;
		}
		else
		{
			cout << "Error:\tParameter DistributionLenghType is not correct" << endl;
			flog.mf_Input("Error:\tParameter DistributionLenghType is not correct\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	m_strCutPosition = "KR";
	m_blr = false;
	m_strProteinLFAQpepPath = m_strQuantiResultPath + "\\ProteinResults";
	m_strProteinsPath = m_strQuantiResultPath + "\\ProteinsInfo.txt";
	m_strRegressionResult = m_strQuantiResultPath + "\\RegressionResults";
	cout << "Done\n";
	flog.mf_Input("Done\n");

}