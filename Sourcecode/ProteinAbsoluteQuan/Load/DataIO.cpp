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
#include"DataIO.h"


#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif

string CLoadIO::mf_GetPeptidesSequenceFromID( int ID)
{
	map<int, string>::iterator Iter;
	Iter = m_Peptides.m_mapPeptideIDAndSequence.find(ID);
	if (Iter != m_Peptides.m_mapPeptideIDAndSequence.end())
		return Iter->second;
	else
	{
		string str = fInt2String(ID);
		flog.mf_Input("\tWarning:\tCannot find peptide " + str +" in  filtered peptides! It may be shared peptide\n");
		return "null";
	}
}
double CLoadIO::mf_GetPeptidesIntensityFromID( int ID, int iExperimentIndex)
{
	map<int, vector<double>>::iterator Iter;
	Iter = m_Peptides.m_mapPeptideIDAndIntensitys.find(ID);
	if (Iter != m_Peptides.m_mapPeptideIDAndIntensitys.end())
		return Iter->second.at(iExperimentIndex);
	else
	{
		string str = fInt2String(ID);
		cout << "\tWarning:\tCannot find peptide's intensity " << ID << endl;
		flog.mf_Input("\tWarning:\tCannot find peptide "+str+"'s intensity. \n");
		return 0;
	}
}

double CLoadIO::mf_GetCombinedIntensityFromID(int ID)
{
	map<int, double>::iterator Iter;
	Iter = m_Peptides.m_mapPeptideIDAndCombinedIntensity.find(ID);
	if (Iter != m_Peptides.m_mapPeptideIDAndCombinedIntensity.end())
		return Iter->second;
	else
	{
		string str = fInt2String(ID);
		cout << "\tWarning:\tCannot find peptide " << ID << "'s combined intensity " << endl;
		flog.mf_Input("\tWarning:\tCannot find peptide " + str + "'s combined intensity. \n");
		return 0;
	}
}

bool CLoadIO::mf_LoadProteinIDAndSequence( string strProteinFasterFilePath, FastaType fastatype)
{
	char Buffer[BUFFERLENGTH];
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, strProteinFasterFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << strProteinFasterFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + strProteinFasterFilePath + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	char *pstr;
	string strIDTemp;
	string strSequenceTemp;
	string strTemp;
	fgets(Buffer, BUFFERLENGTH, pFile);
	while ((Buffer[0] == '\n') && (!feof(pFile)))  //allow empty row;
	{
		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	
	std::match_results<std::string::const_iterator> Matchresult;
	bool valid;
	string strFirstRow;
	pstr = Buffer;
	if (Buffer[0] == '>')  //get the protein id according to the regular expression fastatype;
	{
		strFirstRow = Buffer;
		valid = std::regex_search(strFirstRow, Matchresult, fastatype);
		if (valid == true)
		{
			strIDTemp = Matchresult[1];
		}
		else
		{
			cout << "Error:\tCannot parse the fasta file by the regular expression"<< endl;
			flog.mf_Input("Error:\tCannot parse the fasta file by the regular expression\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	fgets(Buffer, BUFFERLENGTH, pFile);
	while (Buffer[0] != '>')  // get the sequence of the protein id;
	{

		strTemp = Buffer;
		strSequenceTemp = strSequenceTemp + strTemp.substr(0, strTemp.size() - 1);
		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	while (!feof(pFile))
	{
		if (Buffer[0] == '\0')   //allow the blank lines
			continue;
		pstr = Buffer;

		if (Buffer[0] == '>')  //get the protein id according to the regular expression fastatype;
		{
			m_mapProteinIDAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));
			strSequenceTemp.clear();
			strFirstRow = Buffer;
			valid = std::regex_search(strFirstRow, Matchresult, fastatype);
			if (valid == true)
			{
				strIDTemp = Matchresult[1];
			}
			else
			{
				cout << "Error:\tCannot parse the fasta file by the regular expression" << endl;
				flog.mf_Input("Error:\tCannot parse the fasta file by the regular expression\n");
				flog.mf_Destroy();
				exit(1);
			}

			fgets(Buffer, BUFFERLENGTH, pFile);
		}
		else // get the sequence of the protein id;
		{
			strTemp.clear();
			strTemp = pstr;
			strSequenceTemp = strSequenceTemp + strTemp.substr(0, strTemp.size() - 1);
			fgets(Buffer, BUFFERLENGTH, pFile);
		}
	}
	m_mapProteinIDAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));

	fclose(pFile);
	return 1;
}

double CLoadIO::mf_GetMedianIntensity(vector<double> PeptidesIntensity)
{
	sort(PeptidesIntensity.begin(), PeptidesIntensity.end());
	int len = PeptidesIntensity.size();
	if (len == 0)             
		return 0.0;
	if (len % 2 == 0)
		return (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at((len - 2) / 2)) / 2;
	else
		return PeptidesIntensity.at(floor(len / 2));
}
void CLoadIO::mf_saveProteins( CLoadParam param)
{
	ofstream ofile(param.m_strExtractedProteinsPath.c_str());
	if (!ofile)
	{
		cout << "Error:\tCannot open " << param.m_strExtractedProteinsPath << endl;
		flog.mf_Input("Error:\tCannot open " +param.m_strExtractedProteinsPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "Proteins and peptides saved to " << param.m_strQuantiResultPath << endl;
		flog.mf_Input("Proteins and peptides saved to " +param.m_strQuantiResultPath + "\n");
	}
	ofile << "ProteinID\tMajority protein IDs\tPeptidesNumber\tvPeptidesSequence\tIntensity\t";
	map<string, string>::iterator mapExperimentIter;
	mapExperimentIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentIter++)
	{
		ofile << mapExperimentIter->second << "\t";
	}
	ofile << "bIfiBAQExist\t";
	if (param.m_bIfiBAQIntensityExist)
	{
		mapExperimentIter = m_mapExperimentNameAndPeptideIntensityName.begin();
		for (; mapExperimentIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentIter++)
		{
			ofile << "iBAQ " + mapExperimentIter->first + "\t";
		}

	}

	ofile << "\n";

	vector<CProtein>::iterator  ProteinIter;
	vector<string>::iterator vPeptideIter;
	vector<double>::iterator vPeptideIntensityIter;
	vector<double>::iterator vPeptideCombinedIntensityIter;
	int i = 0;
	map<string, vector<double>>::iterator mapExperimentAndPeptideIntensitysIter;
	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
	{
		ofile << ProteinIter->m_strProteinID << "\t" << ProteinIter->m_strProteinFullName << "\t";
		ofile << ProteinIter->m_iPeptidesNumber << "\t";
		for (vPeptideIter = ProteinIter->m_vPeptidesSequences.begin(); vPeptideIter != ProteinIter->m_vPeptidesSequences.end(); vPeptideIter++)
			ofile << *vPeptideIter << ";";
		ofile << "\t";
		vPeptideCombinedIntensityIter = ProteinIter->m_vPeptidesCombinedIntensity.begin();
		for (; vPeptideCombinedIntensityIter != ProteinIter->m_vPeptidesCombinedIntensity.end(); vPeptideCombinedIntensityIter++)
		{
			ofile << fixed << std::setprecision(6) << *vPeptideCombinedIntensityIter << ";";
		}
		ofile << "\t";
		mapExperimentAndPeptideIntensitysIter = ProteinIter->m_mapExperimentAndPeptidesIntensity.begin();
		for (; mapExperimentAndPeptideIntensitysIter != ProteinIter->m_mapExperimentAndPeptidesIntensity.end(); mapExperimentAndPeptideIntensitysIter++)
		{
			for (i = 0; i < mapExperimentAndPeptideIntensitysIter->second.size(); i++)
			{
				ofile << fixed << std::setprecision(6) << mapExperimentAndPeptideIntensitysIter->second.at(i) << ";";
			}
			ofile << "\t";
		}

		if (param.m_bIfiBAQIntensityExist)
		{
			ofile << 1 << "\t";
			for (i = 0; i < ProteinIter->m_vecdMaxQuantiBAQ.size(); i++)
			{
				ofile<<fixed <<std::setprecision(6)<< ProteinIter->m_vecdMaxQuantiBAQ.at(i) << "\t";
			}
		}
		else
		{
			ofile << 0 << "\t";
		}

		ofile << "\n";

	}
	ofile.close();
}

CLoadMaxQuantIO::CLoadMaxQuantIO(const CLoadParam loadParam)
{
	m_NumberOfPeptides = 0;
	mf_GetAttributesName(loadParam.m_strExprimentDesignPath);
}

/*mf_GetAttributesName
Because some attribute name change with the experiment design,
so we need determine the attribute name according to the experiment design file.
*/
void CLoadMaxQuantIO::mf_GetAttributesName(string ExperimentDesignPath)
{
	ifstream fin(ExperimentDesignPath.c_str());
	if (!fin)
	{
		cout << "Error:\tCannot open file " << ExperimentDesignPath << endl;
		flog.mf_Input("Error:\tCannot open file " +ExperimentDesignPath+"\n");
		flog.mf_Destroy();
		exit(2);
	}
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	string strLine;

	getline(fin, strLine);
	GetAttributesFromFirstRow(strLine, mapAttrtibuteAndcolumns, "\t");
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Experiment");
	int iExperimentColumn = mapAttrtibuteAndcolumnsIter->second;

	map<string, int> mapExperiments;
	map<string, int>::iterator mapExperimentsIter;
	vector<string> vecStrTemps;

	while (getline(fin, strLine))
	{
		vecStrTemps.clear();
		vecStrTemps = split(strLine, "\t");
		mapExperiments.insert(pair<string, int>(vecStrTemps[iExperimentColumn], 0));
	}

	m_iNumberOfExperiments = mapExperiments.size();

	string strPeptideIntensityName;
	string striBAQ_IntensityName;

	mapExperimentsIter = mapExperiments.begin();
	for (; mapExperimentsIter != mapExperiments.end(); mapExperimentsIter++)
	{
		strPeptideIntensityName = "Intensity " + mapExperimentsIter->first;
		m_mapExperimentNameAndPeptideIntensityName.insert(pair<string, string>(mapExperimentsIter->first, strPeptideIntensityName));
		striBAQ_IntensityName = "iBAQ " + mapExperimentsIter->first;
		m_mapExperimentNameAndiBAQIntensityName.insert(pair<string, string>(mapExperimentsIter->first, striBAQ_IntensityName));
	}						   

	fin.close();
}
/*mf_GetPeptideAttributecolumn
    determine the column number of peptide attributes which we need
	 "Sequence"��"id ", "Proteins",strPeptideIntensityName
*/
int CLoadMaxQuantIO::mf_GetPeptideAttributecolumn( std::string AttributeName)
{
	if (AttributeName == "Sequence")
		return 0;
	else if (AttributeName == "id")
		return 1;
	else if (AttributeName == "Proteins")
		return 2;
	int icolumn = 2;
	map<string, string>::iterator mapExperimentAndPeptideIntensityIter;
	mapExperimentAndPeptideIntensityIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentAndPeptideIntensityIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentAndPeptideIntensityIter++)
	{
		if (AttributeName == mapExperimentAndPeptideIntensityIter->second)
		{
			icolumn++;
			return icolumn;
		}
	}
	 
	cout<<"Error:\tCannot find the column: "<<AttributeName<<" in the input file"<<endl;
	flog.mf_Input("Error:\tCannot find the column: " + AttributeName + " in the input file.\n");
	return -1;

}

/*mf_LoadPeptides
extract peptides information from peptides.txt ,and save it in m_Peptides.
informations include:
"Sequence", "id", "Proteins", "Leading razor protein", peptide intensity;
deleting constraint:
1, column " Proteins"  has more than one protein;
*/

bool CLoadMaxQuantIO::mf_LoadPeptides(CLoadParam param)
{
	string strProteinsSequenceTemp;
	map<string, string>::iterator ProteinIDAndSequenceIter=m_mapProteinIDAndSequence.begin();
	int iPositionTemp = 0, nextPosition = 0;
	
	for (; ProteinIDAndSequenceIter != m_mapProteinIDAndSequence.end(); ProteinIDAndSequenceIter++)
	{
		strProteinsSequenceTemp = strProteinsSequenceTemp +"?"+ ProteinIDAndSequenceIter->second;
	}

	int ibeginTemp = 0,icutPosition = 0;
	string ProteinId1;

	cout << "Loading peptides Peptides.txt\n";
	flog.mf_Input("Loading peptides from Peptides.txt\n");
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, param.m_strPeptideFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << param.m_strPeptideFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strPeptideFilePath + ". The path of Peptides.txt is needed here.\n" );
		flog.mf_Destroy();
		exit(1);

	}
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	fgets(Buffer, BUFFERLENGTH, pFile);    

	int iSequencecolumnNum = 0, iIDcolumnNum = 0, iProteinsNamecolumnNum = 0, iIntensitycolumnNum = 0;
	bool bIfFindSequence = false, bIfFindID = false, bIfFindProteinNames = false, bIffindIntensity = false;
	pstr = Buffer;
	int icolumns = 0, count = 0;
	pstr1 = strstr(pstr, "\t");
	while (count<4&&pstr1!=NULL)
	{
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			flog.mf_Input("Error:\tThe format of \"peptides.txt\" is wrong.\n");
			flog.mf_Destroy();
			exit(1);
		}
		switch (mf_GetPeptideAttributecolumn(pstr))
		{
		case 0:
			iSequencecolumnNum = icolumns;
			bIfFindSequence = true;
			count++;
			break;
		case 1:
			iIDcolumnNum = icolumns;
			bIfFindID = true;
			count++;
			break;
		case 2:
			iProteinsNamecolumnNum = icolumns;
			bIfFindProteinNames = true;
			count++;
			break;
		case 3:
			iIntensitycolumnNum = icolumns;
			bIffindIntensity = true;
			count++;
			break;
		}
		icolumns++;
		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
	}
	if (count < 4)  
	{
		if (!bIfFindID)
		{
			cout << "Error:\tCannot find the column \"id\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"id\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIffindIntensity)
		{
			cout << "Error:\tCannot find the column PeptideIntensity in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column PeptideIntensity in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindProteinNames)
		{
			cout << "Error:\tCannot find the column \"Proteins\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Proteins\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindSequence)
		{
			cout << "Error:\tCannot find the column \"Sequence\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Sequence\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	icolumns = 0;
	count = 0;

	bool bIfDelete;
	string strSequenceTemp;
	double dIntensityTemp;
	int iIDTemp;
	map<string, string>::iterator mapProteinIDAndSequenceTempIter; 
	int iNormalProteinsNumber;  
	bool b_ifExistDecoyProtein = true;
	bool b_ifExistcontaminantProtein = true;
	string strProteinNameTemp;
	char * pstr2; 
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	while (!feof(pFile))
	{
		bIfDelete = false;
		count = 0;
		icolumns = 0;
		
		while (count < 4)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				flog.mf_Input("Error:\tThe format of peptides.txt  is wrong.\n");
				flog.mf_Destroy();
				exit(1);
			}
			if (icolumns == iSequencecolumnNum)
			{
				strSequenceTemp = pstr;
				count++;
			}
			if (icolumns == iIDcolumnNum)
			{
				iIDTemp = atoi(pstr);
				count++;
			}
			if (icolumns == iProteinsNamecolumnNum)
			{
				if (b_ifExistDecoyProtein)
					if (strstr(pstr, param.m_strDecoyProteinIDPrefix.c_str()) != NULL || (pstr == ""))
						bIfDelete = true;                      // delete the peptide, if its corresponding protein group contain decoy protein
				
				if (!b_ifExistcontaminantProtein)  //do not contain contamination proteins
				{  
					if (strstr(pstr, ";") != NULL) //remove the shared peptides 
						bIfDelete = true;
				}
				else   //may contain contamination proteins
				{
					iNormalProteinsNumber = 0;
					pstr2 = strstr(pstr, ";");
					while (pstr2 != NULL)
					{

						*pstr2 = '\0';
						strProteinNameTemp = pstr;
						strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length()); //change the constant 4 to the length of strContaminatefix by CC 20150727
						if (strProteinNameTemp != param.m_strContaminantfix)
							iNormalProteinsNumber++;
						pstr = pstr2 + 1;
						pstr2 = strstr(pstr, ";");
					}
					strProteinNameTemp = pstr;
					strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
					if ((strProteinNameTemp != param.m_strContaminantfix) && (strProteinNameTemp != ""))
						iNormalProteinsNumber++;

					if (iNormalProteinsNumber != 1)
						bIfDelete = true;
				}
				
				count++;
			}
			if (icolumns == iIntensitycolumnNum)
			{
				dIntensityTemp = atof(pstr);
				count++;
			}
			icolumns++;
			pstr = pstr1 + 1;
		}

		if (!bIfDelete)
		{
			ibeginTemp = 0;
			iPositionTemp = strProteinsSequenceTemp.find(strSequenceTemp, ibeginTemp);
			if (iPositionTemp == strProteinsSequenceTemp.npos)
			{
				cerr << "Error:\tThere is not peptide " << strSequenceTemp << " in Protein fasta" << endl;
				flog.mf_Input("Error:\tThere is not peptide " + strSequenceTemp +" in Protein fasta.\n");
				flog.mf_Destroy();
				exit(1);
			}
			else
			{
				icutPosition = strProteinsSequenceTemp.find("?", iPositionTemp + 1);
				if (icutPosition != strProteinsSequenceTemp.npos)
				{
					iPositionTemp = strProteinsSequenceTemp.find(strSequenceTemp, icutPosition + 1);
					if (iPositionTemp != strProteinsSequenceTemp.npos)
						bIfDelete = true;
				}
			}
		}
	 
		if (!bIfDelete)
		{
			m_Peptides.m_mapPeptideIDAndSequence.insert(pair<int, string>(iIDTemp, strSequenceTemp));	
		}
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	fclose(pFile);
	flog.mf_Input("Done\n");
	return 1;
}

/*mf_LoadPeptidesWithoutRedunPeptides
extract peptides information from peptides.txt ,and save it in m_Peptides. but without redundancy removal 
informations include:
"Sequence", "id", "Proteins", "Leading razor protein", peptide intensity;
deleting constraint:
1, column " Proteins"  has more than one protein;
*/
bool CLoadMaxQuantIO::mf_LoadPeptidesWithoutRedunPeptides( CLoadParam param)
{
	cout << "Loading peptides Peptides.txt\n";
	flog.mf_Input("Loading peptides from Peptides.txt\n");

	// open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, param.m_strPeptideFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << param.m_strPeptideFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strPeptideFilePath + ". The path of Peptides.txt is needed here.\n");
		flog.mf_Destroy();
		exit(1);
	}


	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	//get the columns of sequence、ID、Proteins、Intensity UPS2_yeast by analysising the first row;
	fgets(Buffer, BUFFERLENGTH, pFile);

	int iSequencecolumnNum = 0, iIDcolumnNum = 0;
	int iProteinsColumnNum = 0;
	int iCombinedIntensityColumnNum = 0;

	int iContaminantColumnNum = 0;
	int iReverseColumnNum = 0;
	vector<int> veciIntensitycolumnNum;
	bool bIfFindSequence = false, bIfFindID = false;
	bool bIfFindReverseNames = false;
	bool bIfFindProteins = false;
	bool bIfFindContaminant = false;
	vector<bool> vecbIffindIntensity;
	int icolumns = 0, count = 0;
	int iNumberOfcolumns = 6 + m_iNumberOfExperiments; // Todo
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;

	pstr = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Sequence");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iSequencecolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindSequence = true;
		count++;
	}
	else
	{
		bIfFindSequence = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("id");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iIDcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindID = true;
		count++;
	}
	else
	{
		bIfFindID = false;
	}


	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Reverse");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iReverseColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindReverseNames = true;
		count++;
	}
	else
	{
		bIfFindReverseNames = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Proteins");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinsColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindProteins = true;
		count++;
	}
	else
	{
		bIfFindProteins = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Potential contaminant");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iContaminantColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindContaminant = true;
		count++;
	}
	else
	{
		bIfFindContaminant = false;
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Intensity");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iCombinedIntensityColumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Error:\tCannot find the column \"Intensity\" in the peptides.txt." << endl;
		flog.mf_Input("Error:\tCannot find the column \"Intensity\" in the peptides.txt.\n");
		flog.mf_Destroy();
		exit(1);
	}

	map<string, string>::iterator mapExperimentAndPeptideIntensityIter;
	mapExperimentAndPeptideIntensityIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentAndPeptideIntensityIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentAndPeptideIntensityIter++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find(strim(mapExperimentAndPeptideIntensityIter->second));
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			veciIntensitycolumnNum.push_back(mapAttrtibuteAndcolumnsIter->second);
			vecbIffindIntensity.push_back(true);
			count++;
		}
		else
		{
			vecbIffindIntensity.push_back(false);
		}
	}

	if (count < iNumberOfcolumns)
	{
		if (!bIfFindID)
		{
			cout << "Error:\tCannot find the column \"id\" in the peptides.txt." << endl;
			flog.mf_Input("Error:\tCannot find the column \"id\" in the peptides.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindReverseNames)
		{
			cout << "Error:\tCannot find the column \"Reverse\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Reverse\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindProteins)
		{
			cout << "Error:\tCannot find the column \"Proteins\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Proteins\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindContaminant)
		{
			cout << "Error:\tCannot find the column \"Potential contaminant\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Potential contaminant\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindSequence)
		{
			cout << "Error:\tCannot find the column \"Sequence\" in the peptides.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Sequence\" in the peptides.txt\n");
			flog.mf_Destroy();
			exit(1);
		}
		for (int i = 0; i < vecbIffindIntensity.size(); i++)
		{
			if (!vecbIffindIntensity.at(i))
			{
				cout << "Error:\tCannot find the column PeptideIntensityName in the peptides.txt" << endl;
				flog.mf_Input("Error:\tCannot find the column PeptideIntensityName in the peptides.txt\n");
				flog.mf_Destroy();
				exit(1);

			}
		}

	}

	icolumns = 0;
	count = 0;
	bool bIfDelete;
	string strSequenceTemp;
	string strUniqueTemp;
	double dCombinedIntensitytemp;
	vector<double> vecdIntensityTemp;
	int iIDTemp;
	map<string, string>::iterator mapProteinIDAndSequenceTempIter;
	string strProteinNameTemp;
	string strzLeadRazorProteinsTemp;
	string strReverseTemp;
	string strContaminantTemp;
	char * pstr2;
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	while (!feof(pFile))
	{
		bIfDelete = false;
		count = 0;
		icolumns = 0;
		vecdIntensityTemp.clear();
		while (count < iNumberOfcolumns)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				flog.mf_Input("Error:\tThe format of peptides.txt is wrong.\n");
				flog.mf_Destroy();
				exit(1);
			}
			if (icolumns == iSequencecolumnNum)
			{
				strSequenceTemp = pstr;
				count++;
			}
			if (icolumns == iIDcolumnNum)
			{
				iIDTemp = atoi(pstr);
				count++;
			}
			if (icolumns == iReverseColumnNum)
			{
				strReverseTemp = pstr;
				if (strReverseTemp != "")
					bIfDelete = true;
				count++;
			}
			if (icolumns == iProteinsColumnNum)
			{
				/*if (param.m_bIfExistContanProtein)
				{
				iNormalProteinsNumber = 0;
				pstr2 = strstr(pstr, ";");
				while (pstr2 != NULL)
				{
				*pstr2 = '\0';
				strProteinNameTemp = pstr;
				strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
				if (strProteinNameTemp != param.m_strContaminantfix)
				iNormalProteinsNumber++;
				pstr = pstr2 + 1;
				pstr2 = strstr(pstr, ";");
				}
				strProteinNameTemp = pstr;
				strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
				if (strProteinNameTemp != param.m_strContaminantfix)
				iNormalProteinsNumber++;
				if (iNormalProteinsNumber == 0)
				{
				bIfDelete = true;
				}
				}
				else
				{
				strProteinNameTemp = pstr;
				if (strProteinNameTemp == "")
				bIfDelete = true;
				}
				*/

				/*
				if (param.m_bIfExistDecoyProtein)
				if (strstr(pstr, param.m_strDecoyProteinIDPrefix.c_str()) != NULL || (pstr == ""))
				bIfDelete = true;

				// delete the peptide, if its corresponding protein group contain decoy protein
				if (!param.m_bIfExistContanProtein)
				{
				if (strstr(pstr, ";") != NULL) //remove the shared peptides
				bIfDelete = true;
				}
				else
				{
				iNormalProteinsNumber = 0;
				pstr2 = strstr(pstr, ";");
				while (pstr2 != NULL)
				{
				*pstr2 = '\0';
				strProteinNameTemp = pstr;
				strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
				if (strProteinNameTemp != param.m_strContaminantfix)
				iNormalProteinsNumber++;
				pstr = pstr2 + 1;
				pstr2 = strstr(pstr, ";");
				}
				strProteinNameTemp = pstr;
				strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
				if (strProteinNameTemp != param.m_strContaminantfix)
				iNormalProteinsNumber++;

				if (iNormalProteinsNumber != 1)
				bIfDelete = true;
				}
				*/

				count++;
			}
			if (icolumns == iContaminantColumnNum)
			{
				strContaminantTemp = pstr;
				if (strContaminantTemp != "")
					bIfDelete = true;
				count++;
			}
			if (icolumns == iCombinedIntensityColumnNum)
			{
				dCombinedIntensitytemp = atof(pstr);
				count++;
			}
			for (int i = 0; i < veciIntensitycolumnNum.size(); i++)
			{
				if (icolumns == veciIntensitycolumnNum.at(i))
				{
					vecdIntensityTemp.push_back(atof(pstr));
					count++;
				}
			}
			icolumns++;
			pstr = pstr1 + 1;
		}

		if (!bIfDelete)
		{
			m_Peptides.m_mapPeptideIDAndIntensitys.insert(pair<int, vector<double>>(iIDTemp, vecdIntensityTemp));
			m_Peptides.m_mapPeptideIDAndSequence.insert(pair<int, string>(iIDTemp, strSequenceTemp));
			m_Peptides.m_mapPeptideIDAndCombinedIntensity.insert(pair<int, double>(iIDTemp, dCombinedIntensitytemp));
		}
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	fclose(pFile);
	cout << "\t" << m_Peptides.m_mapPeptideIDAndSequence.size() << " peptides loaded\n";
	int iSizeTemp = m_Peptides.m_mapPeptideIDAndSequence.size();
	string str = fInt2String(iSizeTemp);
	flog.mf_Input("\t" + str + " peptides loaded\n");
	flog.mf_Input("Done\n");
	return 1;
}

/*mf_LoadProteins
extract proteins information from proteinGroups.txt ,and save it in m_vecProteins.
informations include:
"Majority protein IDs", "Unique peptides", "iBAQ UPS2_...",
"Reverse", "Potential contaminant", "Peptide IDs", "Peptide is razor",
deleting constraint:
1, column " Reverse" contain "+";
3, proteinID is P07339 and P08311;
4, "Peptide is razor" do not contain true;
5, peptide intensity equal or less than 0, delete it; 
   if the number of protein's peptides don't  greater than 0, delete it;
*/
bool CLoadMaxQuantIO::mf_LoadProteins(CLoadParam& param)
{
	cout << "Loading protein from  proteinGroups.txt\n";
	flog.mf_Input("Loading protein from  proteinGroups.txt\n");

	//open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, param.m_strProteinFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << param.m_strProteinFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strProteinFilePath + "\n");
		flog.mf_Destroy();
		exit(1);
	}

	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	int icolumns = 0;
	int count = 0;

	int iMajorProteinIDcolumnNum = 0;
	vector<int> veciBAQcolumnNum;
	int iDecoycolumnNum = 0;
	int iPeptideIDscolumnNum = 0;
	int iPeptideRazorcolumnNum = 0;
	int iOnlyIdentBySiteColumnNum = 0;

	bool bIfFindProteinID = false;
	vector<bool> vecbIfFindiBAQ;
	bool bIFFindDecoy = false;
	bool bIfFindPoteinalContaminant = false;
	bool bIfFindPeptideIDs = false;
	bool bIfFindPeptideRazor = false;
	bool bIfFindOnlySite = false;

	int iAttributeNumbers = 5 +  m_iNumberOfExperiments;// the number of Attributes extracted from proteingroups.txt 


	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Majority protein IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iMajorProteinIDcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindProteinID = true;
		count++;
	}
	else
	{
		bIfFindProteinID = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Reverse");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iDecoycolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIFFindDecoy = true;
		count++;
	}
	else
	{
		bIFFindDecoy = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Peptide IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptideIDscolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindPeptideIDs = true;
		count++;
	}
	else
	{
		bIfFindPeptideIDs = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Peptide is razor");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptideRazorcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindPeptideRazor = true;
		count++;
	}
	else
	{
		bIfFindPeptideRazor = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Only identified by site");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iOnlyIdentBySiteColumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindOnlySite = true;
		count++;
	}
	else
	{
		bIfFindOnlySite = false;
	}

	map<string, string>::iterator mapExperimentNameAndiBAQIntensityNameIter;
	mapExperimentNameAndiBAQIntensityNameIter = m_mapExperimentNameAndiBAQIntensityName.begin();
	for (; mapExperimentNameAndiBAQIntensityNameIter != m_mapExperimentNameAndiBAQIntensityName.end(); mapExperimentNameAndiBAQIntensityNameIter++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find(strim(mapExperimentNameAndiBAQIntensityNameIter->second));
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			param.m_bIfiBAQIntensityExist = true;
			veciBAQcolumnNum.push_back(mapAttrtibuteAndcolumnsIter->second);
			vecbIfFindiBAQ.push_back(true);
			count++;
		}
		else
		{
			vecbIfFindiBAQ.push_back(false);
		}
	}


	if (!param.m_bIfiBAQIntensityExist)  iAttributeNumbers -= m_iNumberOfExperiments;

	if (count < iAttributeNumbers)
	{
		if (!bIfFindPeptideIDs)
		{
			cout << "Error:\tCannot find the column \"Peptide IDs\" in the proteinGroups.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Peptide IDs\" in the proteinGroups.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindPeptideRazor)
		{
			cout << "Error:\tCannot find the column \"Peptide is razor\" in the proteinGroups.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Peptide is razor\" in the proteinGroups.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindPoteinalContaminant)
		{
			cout << "Error:\tCannot find the column \"Potential contaminant\" in the proteinGroups.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Potential contaminant\" in the proteinGroups.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindProteinID)
		{
			cout << "Error:\tCannot find the column \"Majority protein IDs\" in the proteinGroups.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Majority protein IDs\" in the proteinGroups.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIFFindDecoy)
		{
			cout << "Error:\tCannot find the column \"Reverse\" in the proteinGroups.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Reverse\" in the proteinGroups.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (!bIfFindOnlySite)
		{
			cout << "Error:\tCannot find the column \"Only identified by site\" in the proteinGroups.txt" << endl;
			flog.mf_Input("Error:\tCannot find the column \"Only identified by site\" in the proteinGroups.txt.\n");
			flog.mf_Destroy();
			exit(1);
		}
		if (param.m_bIfiBAQIntensityExist)
		{
			for (int i = 0; i < vecbIfFindiBAQ.size(); i++)
			{
				if (!vecbIfFindiBAQ.at(i))
				{
					cout << "Error:\tCannot find the column iBAQPeptideIntensityName in the proteinGroups.txt" << endl;
					flog.mf_Input("Error:\tCannot find the column iBAQPeptideIntensityName in the proteinGroups.txt.\n");
					flog.mf_Destroy();
					exit(1);

				}
			}
		}

	}

	icolumns = 0;
	count = 0;
	bool bIfDelete;
	CProtein cprotein;
	string strProteinIDTemp;
	string strSequenceTemp;
	vector<double> vecdiBAQ_IntensityTemp;
	string strPeptidesIDTemp;
	string strPeptidesIDRazor;
	string strOnlyIdentBySite;
	string strtemp;
	int itemp;
	double dIntensityTemp = 0.0;
	vector<double> vecdIntensityTemp;
	int begin1 = 0, end1 = 0;
	int begin2 = 0, end2 = 0;
	int NumberOfSharePeptides = 0;
	map<string, string>::iterator ProteinIDSequenceIter;

	int iNormalProteinsNumber;
	bool b_ifExistDecoyProtein = true;
	bool b_ifExistcontaminantProtein = true;
	string strProteinNameTemp;
	char * pstr2;
	bool bIfStandProFirst;;
	int i = 0;
	bool bIfiBAQAllZeros;
	double diBAQTemp;

	map<string, string>::iterator mapExperimentNamesIter;
	map<string, vector<double>>::iterator mapExperimentAndIntensitysIter;
	int iExperimentIndex = 0;


	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	map<string, bool>::iterator mapSequenceAndIfUseIter;
	while (!feof(pFile))
	{
		bIfDelete = false;
		bIfiBAQAllZeros = true;
		count = 0;
		icolumns = 0;
		vecdiBAQ_IntensityTemp.clear();
		while (count < iAttributeNumbers)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				flog.mf_Input("Error:\tThe format of \"proteinGroups.txt\" is wrong.\n");
				flog.mf_Destroy();
				exit(1);
			}

			if (icolumns == iMajorProteinIDcolumnNum)
			{				//delete the protein group if it contains decoy protein or it only contains contaminate proteins. If left, take the first protein as its proxy.

				cprotein.m_strProteinFullName = pstr;
				if (b_ifExistDecoyProtein)
				if (strstr(pstr, param.m_strDecoyProteinIDPrefix.c_str()) != NULL)
					bIfDelete = true;                      // delete the peptide, if its corresponding protein group contain reverse protein

				bIfStandProFirst = false;
				if (!b_ifExistcontaminantProtein)
				{
					if (bIfStandProFirst)
					{
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							while (pstr2 != NULL)
							{
								*pstr2 = '\0';
								strProteinIDTemp = pstr;
								if ((strProteinIDTemp.find("ups") < strProteinIDTemp.size()) || (strProteinIDTemp.find("UPS") < strProteinIDTemp.size()))
								{
									break;
								}
								pstr = pstr2 + 1;
								pstr2 = strstr(pstr, ";");
							}
						}
						else
						{
							strProteinIDTemp = pstr;
						}
					}
					else
					{
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							*pstr2 = '\0';
							strProteinIDTemp = pstr;
							pstr = pstr2 + 1;
						}
						else
						{
							strProteinIDTemp = pstr;
						}
					}

				}
				else
				{
					iNormalProteinsNumber = 0;
					if (bIfStandProFirst)
					{
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							while (pstr2 != NULL)
							{

								*pstr2 = '\0';
								strProteinNameTemp = pstr;
								if ((strProteinNameTemp.find("ups") < strProteinNameTemp.size()) || (strProteinNameTemp.find("UPS") < strProteinNameTemp.size()))
								{
									strProteinIDTemp = pstr;
								}
								strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
								if (strProteinNameTemp != param.m_strContaminantfix)
									iNormalProteinsNumber++;

								pstr = pstr2 + 1;
								pstr2 = strstr(pstr, ";");
							}
							strProteinNameTemp = pstr;
							if ((strProteinNameTemp.find("ups") < strProteinNameTemp.size()) || (strProteinNameTemp.find("UPS") < strProteinNameTemp.size()))
							{
								strProteinIDTemp = pstr;
							}
							strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
							if (strProteinNameTemp != param.m_strContaminantfix)
								iNormalProteinsNumber++;
						}
						else
						{   //only contain one protein
							strProteinNameTemp = pstr;
							strProteinIDTemp = pstr;
							strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
							if (strProteinNameTemp != param.m_strContaminantfix)
								iNormalProteinsNumber++;
						}

						if (iNormalProteinsNumber == 0)
							bIfDelete = true;
					}
					else
					{ // no UPS2, get the first ID 
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							while (pstr2 != NULL)
							{
								*pstr2 = '\0';
								strProteinNameTemp = pstr;
								strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
								if (strProteinNameTemp != param.m_strContaminantfix)
								{
									strProteinIDTemp = pstr;
									iNormalProteinsNumber++;
									break;
								}
								pstr = pstr2 + 1;
								pstr2 = strstr(pstr, ";");
							}
							strProteinNameTemp = pstr;
							strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
							if ((strProteinNameTemp != param.m_strContaminantfix) && (iNormalProteinsNumber == 0))
							{
								strProteinIDTemp = pstr;
								iNormalProteinsNumber++;
							}
						}
						else
						{   //only contain one protein
							strProteinNameTemp = pstr;
							strProteinIDTemp = pstr;
							strProteinNameTemp = strProteinNameTemp.substr(0, param.m_strContaminantfix.length());
							if (strProteinNameTemp != param.m_strContaminantfix)
								iNormalProteinsNumber++;
						}
						if (iNormalProteinsNumber == 0)
							bIfDelete = true;
					}
				}
				count++;
			}

			if (icolumns == iOnlyIdentBySiteColumnNum)
			{
				strOnlyIdentBySite = pstr;
				if (strOnlyIdentBySite != ""){
					bIfDelete = true;
				}
				count++;
			}
			for (i = 0; i < veciBAQcolumnNum.size(); i++)
			{
				if (icolumns == veciBAQcolumnNum.at(i))
				{
					diBAQTemp = atof(pstr);
					if (diBAQTemp != 0.0)
					{
						bIfiBAQAllZeros = false;
					}
					vecdiBAQ_IntensityTemp.push_back(atof(pstr));
					count++;
				}
			}

			if (icolumns == iDecoycolumnNum)
			{
				count++;
			}
			if (icolumns == iPeptideIDscolumnNum)
			{
				strPeptidesIDTemp = pstr;
				count++;
			}
			if (icolumns == iPeptideRazorcolumnNum)
			{
				strPeptidesIDRazor = pstr;
				count++;
			}
			icolumns++;
			pstr = pstr1 + 1;
		}
		if (bIfiBAQAllZeros == true)
			bIfDelete = true;

		if (!bIfDelete)
		{
			cprotein.m_strProteinID = strProteinIDTemp;
			int itemp1 = strProteinIDTemp.find("|");
			int itemp2 = strProteinIDTemp.find("|", itemp1 + 1);
			strProteinIDTemp = strProteinIDTemp.substr(itemp1 + 1, itemp2 - itemp1 - 1);
			if (strProteinIDTemp == "P07339ups" || strProteinIDTemp == "P08311ups")
			{  //delete P07339 and P08311 in UPS2;add by gzhq 20150709
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				continue;
			}
			ProteinIDSequenceIter = m_mapProteinIDAndSequence.find(cprotein.m_strProteinID);
			if (ProteinIDSequenceIter != m_mapProteinIDAndSequence.end())
			{
				cprotein.m_strProteinSequence = ProteinIDSequenceIter->second;
			}
			else
			{
				cout << "\tWarning:\tCannot find protein: " << cprotein.m_strProteinID << " in fasta file." << endl;
				flog.mf_Input("\tWarning:\tCannot find protein: " + cprotein.m_strProteinID + " in fasta file.\n");
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				continue;
			}
			for (i = 0; i < vecdiBAQ_IntensityTemp.size(); i++)
			{
				cprotein.m_vecdMaxQuantiBAQ.push_back(vecdiBAQ_IntensityTemp.at(i));
			}

			begin1 = 0;
			begin2 = 0;

			end1 = strPeptidesIDRazor.find(";", begin1 + 1);
			strtemp.clear();
			while (end1 != strPeptidesIDRazor.npos)
			{
				if (begin1 == 0)
					strtemp = strPeptidesIDRazor.substr(begin1, end1 - begin1);
				else
					strtemp = strPeptidesIDRazor.substr(begin1 + 1, end1 - begin1 - 1);
				end2 = strPeptidesIDTemp.find(";", begin2 + 1);
				if (begin2 == 0)
					itemp = atoi(strPeptidesIDTemp.substr(begin2, end2 - begin2).c_str());
				else
					itemp = atoi(strPeptidesIDTemp.substr(begin2 + 1, end2 - begin2 - end2 - 1).c_str());
				begin1 = end1;
				begin2 = end2;
				transform(strtemp.begin(), strtemp.end(), strtemp.begin(), ::tolower);//change to lowcase
				if (strtemp == "true")
				{
					cprotein.m_vPeptidesIDs.push_back(itemp);
					strtemp = mf_GetPeptidesSequenceFromID(itemp);
					if (strtemp == "null")
					{
						cprotein.m_vPeptidesIDs.pop_back();
						end1 = strPeptidesIDRazor.find(";", begin1 + 1);
						continue;
					}
					mapExperimentNamesIter = m_mapExperimentNameAndPeptideIntensityName.begin();
					iExperimentIndex = 0;
					for (; mapExperimentNamesIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentNamesIter++)
					{
						dIntensityTemp = mf_GetPeptidesIntensityFromID(itemp, iExperimentIndex);
						iExperimentIndex++;
						mapExperimentAndIntensitysIter = cprotein.m_mapExperimentAndPeptidesIntensity.find(mapExperimentNamesIter->first);
						if (mapExperimentAndIntensitysIter == cprotein.m_mapExperimentAndPeptidesIntensity.end())
						{
							vecdIntensityTemp.clear();
							vecdIntensityTemp.push_back(dIntensityTemp);
							cprotein.m_mapExperimentAndPeptidesIntensity.insert(pair<string, vector<double>>(mapExperimentNamesIter->first, vecdIntensityTemp));
						}
						else
						{
							mapExperimentAndIntensitysIter->second.push_back(dIntensityTemp);
						}
					}//end for experiments

					dIntensityTemp = mf_GetCombinedIntensityFromID(itemp);
					cprotein.m_vPeptidesCombinedIntensity.push_back(dIntensityTemp);
					cprotein.m_vPeptidesSequences.push_back(strtemp);

				}

				end1 = strPeptidesIDRazor.find(";", begin1 + 1);

			}// end while
			if (begin2 == 0)
			{  //only one peptide
				strtemp = strPeptidesIDRazor.substr(begin1, strPeptidesIDRazor.size() - begin1);
				itemp = atoi(strPeptidesIDTemp.substr(begin2, strPeptidesIDTemp.size() - begin2).c_str());
			}
			else
			{
				strtemp = strPeptidesIDRazor.substr(begin1 + 1, strPeptidesIDRazor.size() - begin1 - 1);
				itemp = atoi(strPeptidesIDTemp.substr(begin2 + 1, strPeptidesIDTemp.size() - begin2 - 1).c_str());
			}

			transform(strtemp.begin(), strtemp.end(), strtemp.begin(), ::tolower);
			if (strtemp == "true")
			{
				strtemp = mf_GetPeptidesSequenceFromID(itemp);
				if (strtemp != "null")
				{
					cprotein.m_vPeptidesIDs.push_back(itemp);
					mapExperimentNamesIter = m_mapExperimentNameAndPeptideIntensityName.begin();
					iExperimentIndex = 0;
					for (; mapExperimentNamesIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentNamesIter++)
					{
						dIntensityTemp = mf_GetPeptidesIntensityFromID(itemp, iExperimentIndex);
						iExperimentIndex++;
						mapExperimentAndIntensitysIter = cprotein.m_mapExperimentAndPeptidesIntensity.find(mapExperimentNamesIter->first);
						if (mapExperimentAndIntensitysIter == cprotein.m_mapExperimentAndPeptidesIntensity.end())
						{
							vecdIntensityTemp.clear();
							vecdIntensityTemp.push_back(dIntensityTemp);
							cprotein.m_mapExperimentAndPeptidesIntensity.insert(pair<string, vector<double>>(mapExperimentNamesIter->first, vecdIntensityTemp));
						}
						else
						{
							mapExperimentAndIntensitysIter->second.push_back(dIntensityTemp);
						}
					} //end for experiments
					
					dIntensityTemp = mf_GetCombinedIntensityFromID(itemp);
					cprotein.m_vPeptidesCombinedIntensity.push_back(dIntensityTemp);
					cprotein.m_vPeptidesSequences.push_back(strtemp);

				}
			}


			cprotein.m_iPeptidesNumber = cprotein.m_vPeptidesSequences.size();
			if (cprotein.m_vPeptidesSequences.size() > 0)
			{
				m_NumberOfPeptides += cprotein.m_iPeptidesNumber;
				m_vecProteins.push_back(cprotein);
			}
			cprotein.Clear();

		} // end for if (strtemp == "true")

		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	fclose(pFile);
	if (m_vecProteins.size() == 0)
	{
		cout << "Error:\tPlease check the fasta file or the input file!" << endl;
		flog.mf_Input("Error:\tPlease check the fasta file or the input file!\n");
		flog.mf_Destroy();
		exit(1);
	}
	cout << "\t" << m_vecProteins.size() << " Proteins  and " << m_NumberOfPeptides << " peptides Loaded" << endl;
	int iSizeTemp = m_vecProteins.size();
	string strProteinSize = fInt2String(iSizeTemp);
	string strNumberOfPeptides = fInt2String(m_NumberOfPeptides);
	flog.mf_Input("\t" + strProteinSize + " Proteins  and " + strNumberOfPeptides + " peptides Loaded\n");
	flog.mf_Input("Done\n");
	return 1;
}

class XMLCleaner_
{
public:
	explicit XMLCleaner_(XMLHandler * handler) :
		p_(handler)
	{

	}

	~XMLCleaner_()
	{
		p_->reset();
	}

private:
	XMLHandler * p_;
};

CloadmzQuantMLIO::CloadmzQuantMLIO()
{
}

CloadmzQuantMLIO::~CloadmzQuantMLIO()
{
}

void CloadmzQuantMLIO::load( const CLoadParam &param)
{
	int begin = 0;
	begin = param.m_strIdentResultPath.find_last_of(".");
	string SuffixTemp;
	if (begin>0)
		SuffixTemp = param.m_strIdentResultPath.substr(begin + 1, param.m_strIdentResultPath.size() - begin);
	else
	{
		cout << "Error:\tPlease set a correct mzq file path!\n";
		flog.mf_Input("Error:\tPlease set a correct mzq file path!\n");
		flog.mf_Destroy();
		exit(1);
	}
	if (SuffixTemp != "mzq")
	{
		cout << "Error:\tPlease set a specific mzq file path!\n";
		flog.mf_Input("Error:\tPlease set a specific mzq file path!\n");
		flog.mf_Destroy();
		exit(1);
	}
	fstream fTest(param.m_strIdentResultPath);
	if (!fTest)
	{
		cout << "Error:\tPlease set a correct mzq file path!\n";
		flog.mf_Input("Error:\tPlease set a correct mzq file path!\n");
		flog.mf_Destroy();
		exit(1);
	}else
	{
		fTest.close();
	}
	mzQuantMLHandler handler(m_qi_,param.m_strIdentResultPath);
	parse_(param.m_strIdentResultPath, &handler);
	if (param.m_bIfExistDecoyProtein)
	{
		mf_RemovDecoyProtein(m_qi_,param);
	}
	mf_RemovedSharedPeptides(m_qi_);
	mf_CheckMatchedFasta(param);
	mf_CheckEnoughInformation();

}
void CloadmzQuantMLIO::mf_RemovDecoyProtein(CQuantInformation &m_qi, const CLoadParam& param)
{
	map<string, vector<string>> mapProAndPepTemp;

	map<string, vector<string>>::iterator iter;
	for (iter = m_qi.mapProteinAndPeptideIDs.begin(); iter != m_qi.mapProteinAndPeptideIDs.end(); iter++)
	{
		if (!bPrefix(iter->first,param.m_strDecoyProteinIDPrefix))
			mapProAndPepTemp.insert(pair<string, vector<string>>(iter->first, iter->second));
	}
	m_qi.mapProteinAndPeptideIDs.clear();
	for (iter = mapProAndPepTemp.begin(); iter != mapProAndPepTemp.end(); iter++)
	{
		m_qi.mapProteinAndPeptideIDs.insert(pair<string, vector<string>>(iter->first, iter->second));
	}
}

void CloadmzQuantMLIO::mf_RemovedSharedPeptides(CQuantInformation &m_qi)
{
	map<string, vector<string>> mapPeptideIdsAndProteinIds;
	map<string, vector<string>>::iterator mapProteinAndPeptideIDsiter;
	map<string, vector<string>>::iterator mapPeptideIdsAndProteinIdsiter;
	vector<string>::iterator vecPeptideidIter;
	vector<string> vecProteinIdsTemp;
	mapProteinAndPeptideIDsiter = m_qi.mapProteinAndPeptideIDs.begin();
	for (; mapProteinAndPeptideIDsiter != m_qi.mapProteinAndPeptideIDs.end(); mapProteinAndPeptideIDsiter++)
	{
		vecPeptideidIter = mapProteinAndPeptideIDsiter->second.begin();
		for (; vecPeptideidIter != mapProteinAndPeptideIDsiter->second.end(); vecPeptideidIter++)
		{
			mapPeptideIdsAndProteinIdsiter = mapPeptideIdsAndProteinIds.find(*vecPeptideidIter);
			if (mapPeptideIdsAndProteinIdsiter == mapPeptideIdsAndProteinIds.end())
			{// new peptide
				vecProteinIdsTemp.clear();
				vecProteinIdsTemp.push_back(mapProteinAndPeptideIDsiter->first);
				mapPeptideIdsAndProteinIds.insert(pair<string, vector<string>>(*vecPeptideidIter, vecProteinIdsTemp));

			}
			else
			{
				mapPeptideIdsAndProteinIdsiter->second.push_back(mapProteinAndPeptideIDsiter->first);
			}

		}
	}



	map<string, vector<string>> mapProteinAndPeptideIDsTemp;
	mapProteinAndPeptideIDsiter = m_qi.mapProteinAndPeptideIDs.begin();
	for (; mapProteinAndPeptideIDsiter != m_qi.mapProteinAndPeptideIDs.end(); mapProteinAndPeptideIDsiter++)
	{
		mapProteinAndPeptideIDsTemp.insert(pair<string, vector<string>>(mapProteinAndPeptideIDsiter->first, mapProteinAndPeptideIDsiter->second));
	}

	m_qi.mapProteinAndPeptideIDs.clear();
	mapProteinAndPeptideIDsiter = mapProteinAndPeptideIDsTemp.begin();
	vector<string> vecPeptideTemp;
	for (; mapProteinAndPeptideIDsiter != mapProteinAndPeptideIDsTemp.end(); mapProteinAndPeptideIDsiter++)
	{
		vecPeptideidIter = mapProteinAndPeptideIDsiter->second.begin();
		for (; vecPeptideidIter != mapProteinAndPeptideIDsiter->second.end(); vecPeptideidIter++)
		{
			if (mapPeptideIdsAndProteinIds[*vecPeptideidIter].size() == 1)
			{
				if (m_qi.mapProteinAndPeptideIDs.find(mapProteinAndPeptideIDsiter->first) == m_qi.mapProteinAndPeptideIDs.end())
				{// first peptide
					vecPeptideTemp.clear();
					vecPeptideTemp.push_back(*vecPeptideidIter);
					m_qi.mapProteinAndPeptideIDs.insert(pair<string, vector<string>>(mapProteinAndPeptideIDsiter->first, vecPeptideTemp));

				}
				else
				{
					m_qi.mapProteinAndPeptideIDs[mapProteinAndPeptideIDsiter->first].push_back(*vecPeptideidIter);
				}
			}
		}
	}
}
void CloadmzQuantMLIO::mf_CheckMatchedFasta(const CLoadParam &param)
{
	map<string, vector<string>>::iterator mapProteinAndPeptideIdIter;
	map<string, string>::iterator mapProteinIdAndSequenceIter;
	
	
	mapProteinAndPeptideIdIter = m_qi_.mapProteinAndPeptideIDs.begin();
	for (; mapProteinAndPeptideIdIter != m_qi_.mapProteinAndPeptideIDs.end(); mapProteinAndPeptideIdIter++)
	{
		mapProteinIdAndSequenceIter = m_mapProteinIDAndSequence.find(mapProteinAndPeptideIdIter->first);
		if (mapProteinIdAndSequenceIter == m_mapProteinIDAndSequence.end())
		{
			cout << "Error:\tThe identifiers extracted from the fasta file do not match the identifiers extracted from the mzQuantML file! Cannot find protein " << mapProteinAndPeptideIdIter->first<< " in fasta file. Please check it again." << endl;
			flog.mf_Input("Error:\tThe identifiers extracted from the fasta file do not match the identifiers extracted from the mzQuantML file! Cannot find protein " + mapProteinAndPeptideIdIter->first + " in fasta file. Please check it again.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
}
void CloadmzQuantMLIO::mf_CheckEnoughInformation()
{
	if (m_qi_.mapProteinAndPeptideIDs.size()==0)
	{
		cout << "Error:\tCannot find information about proteins and their corresponding peptides in mzQuantML file!" << endl;
		flog.mf_Input("Error:\tCannot find information about proteins and their corresponding peptides in mzQuantML file!\n" );
		flog.mf_Destroy();
		exit(1);
	}
	if (m_qi_.mapPeptideAndExperimenIntensity.size() == 0)
	{
		cout << "Error:\tCannot find information about peptides and their corresponding intensities in mzQuantML file!" << endl;
		flog.mf_Input("Error:\tCannot find information about peptides and their corresponding intensities in mzQuantML file!\n");
		flog.mf_Destroy();
		exit(1);
	}
	if (m_qi_.mapPeptideIdAndSequence.size() == 0)
	{
		cout << "Error:\tCannot find information about peptide sequences in mzQuantML file!" << endl;
		flog.mf_Input("Error:\tCannot find information about peptide sequences in mzQuantML file!\n");
		flog.mf_Destroy();
		exit(1);
	}
}
void CloadmzQuantMLIO::parse_(const string & filename, XMLHandler * handler)
{
	// ensure handler->reset() is called to save memory (in case the XMLFile reader, e.g. FatureXMLFile, is used again)
	XMLCleaner_ clean(handler);

	// initialize parser
	try
	{
		xercesc::XMLPlatformUtils::Initialize();
	}
	catch (const xercesc::XMLException & toCatch)
	{
		//throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + StringManager().convert(toCatch.getMessage()));
	}

	xercesc::SAX2XMLReader * parser = xercesc::XMLReaderFactory::createXMLReader();
	parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces, false);
	parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes, false);

	parser->setContentHandler(handler);
	parser->setErrorHandler(handler);

	xercesc::InputSource * source;

	const char *pChar = filename.c_str();
	XMLCh * pXMLch = xercesc::XMLString::transcode(pChar);
	source = new xercesc::LocalFileInputSource(pXMLch);

	// detailed information given http://xerces.apache.org/xerces-c/apiDocs-3/classParser.html#a624fc687a49b917c11ef632367568b60
	parser->parse(*source);

	delete(parser);
	delete source;
}

void CloadmzQuantMLIO::mf_saveProteins(CLoadParam param)
{
	ofstream ofile(param.m_strExtractedProteinsPath.c_str());
	if (!ofile)
	{
		cout << "Error:\tCannot open " << param.m_strExtractedProteinsPath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strExtractedProteinsPath+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "Proteins and peptides saved to " << param.m_strQuantiResultPath << endl;
		flog.mf_Input("Proteins and peptides saved to " + param.m_strQuantiResultPath+"\n");
	}
	ofile << "ProteinID\tMajority protein IDs\tPeptidesNumber\tvPeptidesSequence\t";
	vector<string>::iterator vecExperimentIter;
	vecExperimentIter = m_qi_.vecExperimentNames.begin();
	for (; vecExperimentIter != m_qi_.vecExperimentNames.end(); vecExperimentIter++)
	{
		ofile << "Intensity " << *vecExperimentIter << "\t";
	}
	ofile << "Intensity" << "\t";
	ofile << "bIfiBAQExist\t";
	ofile << "\n";

	map<string, vector<string>>::iterator mapProteinAndPeptideIDsIter;
	vector<string>::iterator vecPeptidesIter;
	vector<double> dvCombinedPepIntensities;
	mapProteinAndPeptideIDsIter = m_qi_.mapProteinAndPeptideIDs.begin();
	for (; mapProteinAndPeptideIDsIter != m_qi_.mapProteinAndPeptideIDs.end(); mapProteinAndPeptideIDsIter++)
	{
		ofile << mapProteinAndPeptideIDsIter->first << "\t";
		ofile << mapProteinAndPeptideIDsIter->first << "\t";
		ofile << mapProteinAndPeptideIDsIter->second.size() << "\t";
		dvCombinedPepIntensities.clear();
		for (vecPeptidesIter = mapProteinAndPeptideIDsIter->second.begin(); vecPeptidesIter != mapProteinAndPeptideIDsIter->second.end(); vecPeptidesIter++)
		{
			ofile << m_qi_.mapPeptideIdAndSequence[*vecPeptidesIter] << ";";
			dvCombinedPepIntensities.push_back(0.0);
		}
		ofile << "\t";
		int i = 0;
		int j = 0;
		for (i = 0; i < m_qi_.vecExperimentNames.size(); i++)
		{
			for (j=0,vecPeptidesIter = mapProteinAndPeptideIDsIter->second.begin(); vecPeptidesIter != mapProteinAndPeptideIDsIter->second.end(); vecPeptidesIter++,j++)
			{
				dvCombinedPepIntensities[j] += m_qi_.mapPeptideAndExperimenIntensity[*vecPeptidesIter][m_qi_.vecExperimentNames[i]];
				ofile << m_qi_.mapPeptideAndExperimenIntensity[*vecPeptidesIter][m_qi_.vecExperimentNames[i]] << ";";
			}
			ofile << "\t";
		}
		for (j = 0, vecPeptidesIter = mapProteinAndPeptideIDsIter->second.begin(); vecPeptidesIter != mapProteinAndPeptideIDsIter->second.end(); vecPeptidesIter++, j++)
		{
			ofile << dvCombinedPepIntensities[j] << ";";
		}
		ofile << "\t";
		ofile << 0 << "\t";
		ofile << endl;
	}
	ofile.close();
	// Save necessary experiment information
	ofstream fout(param.m_strQuantiResultPath + "\\mzQuantMLExperiment.config");
	fout << "Experiment Number\t" << m_qi_.vecExperimentNames.size() << endl;
	for (int i = 0; i < m_qi_.vecExperimentNames.size(); i++)
	{
		fout << "Experiment id\t" << m_qi_.vecExperimentNames.at(i)<<endl;
	}
	fout.close();
}

void CLoadPeakViewIO::mf_SaveProteinGroups(string path, const map<string, CRelatedProtein>& mapDistinctProIDAndRelatedInfo)
{
	ofstream oFile(path);
	if (!oFile)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " + path + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "\t\tProtein groups saved to " << path << "\n";
		flog.mf_Input("\t\tProteingroups saved to " + path + "\n");
	}
	oFile << "Represent Protein ID\tProtein IDs\tNumber of Proteins\tPeptides\tPeptide is Group Unique\tPeptide Number(All)\t\
Peptide Number(group unique)\n";

	map<string, CRelatedProtein>::const_iterator mapDistinctProIDAndRelatedInfoIter;
	mapDistinctProIDAndRelatedInfoIter = mapDistinctProIDAndRelatedInfo.begin();
	int i = 0;
	int j = 0;
	int iGuniquePeptideNum;
	for (; mapDistinctProIDAndRelatedInfoIter != mapDistinctProIDAndRelatedInfo.end(); mapDistinctProIDAndRelatedInfoIter++)
	{
		oFile << mapDistinctProIDAndRelatedInfoIter->first << "\t";
		
		if (mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.size() < 1)
		{
			cout << "Error:\tProtein " << mapDistinctProIDAndRelatedInfoIter->first << " do not have subset proteins.\n";
		}
		else
		{
			i = 0;
			for (; i < mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.size() - 1; i++)
			{
				oFile << mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.at(i) << ";";
			}
			oFile << mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.at(i) << "\t";
		}
		oFile << mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.size() << "\t";
		j = 0;
		for (j = 0; j < mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.size()-1; j++)
		{
			oFile << mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(j) << ";";
		}
		oFile << mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(j) << "\t";

		iGuniquePeptideNum = 0;
		j = 0;
		for (; j < mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.size()-1; j++)
		{
			if (mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.at(j))
			{
				iGuniquePeptideNum++;
				oFile << "True;";
			}
			else
			{
				oFile << "False;";
			}
		}
		if (mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.at(j))
		{
			iGuniquePeptideNum++;
			oFile << "True\t";
		}
		else
		{
			oFile << "False\t";
		}
		oFile << mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.size() << "\t";
		oFile << iGuniquePeptideNum << "\t";
		oFile << endl; 			
	}

	oFile.close();
}

void CLoadPeakViewIO::mf_SaveSharedPeptides(string path, const map<string, vector<string>>& mapSharedPeptidesAndProteins)
{
	ofstream oFile(path);
	if (!oFile)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " + path + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "\t\tShared peptides saved to " << path << "\n";
		flog.mf_Input("\t\tShared peptides saved to " + path + "\n");
	}

	oFile << "Shared Peptide\tProtein Groups" << endl;
	map<string, vector<string>>::const_iterator mapSharedPeptidesAndProteinsIter;
	mapSharedPeptidesAndProteinsIter = mapSharedPeptidesAndProteins.begin();
	for (; mapSharedPeptidesAndProteinsIter != mapSharedPeptidesAndProteins.end(); mapSharedPeptidesAndProteinsIter++)
	{
		oFile << mapSharedPeptidesAndProteinsIter->first << "\t";
		for (int i = 0; i < mapSharedPeptidesAndProteinsIter->second.size(); i++)
		{
			oFile << mapSharedPeptidesAndProteinsIter->second.at(i) << ";";
		}
		oFile << endl;
	}

	oFile.close();

}

bool CLoadPeakViewIO::mf_GetPeakViewAttributeColumns(const map<string, int> &mapAttrtibuteAndcolumns, int &iPeptidecolumnNum, int &iProteinsNamecolumnNum,
	vector<int> &vecIntensitiesColumn, int &iPrecursorChargeColumnNum)
{
	map<string, int>::const_iterator mapAttrtibuteAndcolumnsIter;
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Peptide");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptidecolumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"Peptide\" in the SWATH 2.0 result file." << endl;
		flog.mf_Input("Error:\tCannot find column \"Peptide\" in the SWATH 2.0 result file.\n");
		return false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Protein");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinsNamecolumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"Protein\" in the SWATH 2.0 result file." << endl;
		flog.mf_Input("Error:\tCannot find column \"Protein\" in the SWATH 2.0 result file.\n");
		return false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Precursor Charge");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPrecursorChargeColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Error:\tCannot find column \"Precursor Charge\" in the SWATH 2.0 result file." << endl;
		flog.mf_Input("Error:\tCannot find column \"Precursor Charge\" in the SWATH 2.0 result file.\n");
		return false;
	}

	int i = 5;
	for (; i < mapAttrtibuteAndcolumns.size(); i++)
	{
		vecIntensitiesColumn.push_back(i);
	}

	return true;

}
void CLoadPeakViewIO::mf_GetPeptidesAndProteins(const CLoadParam &param,
	map<string, vector<string>>& mapPeptidesAndProteinNames, map<string, vector<string>>& mapProteinAndPeptides)
{
	// open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, param.m_strIdentResultPath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << param.m_strIdentResultPath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strIdentResultPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}

	int iSequencecolumnNum = 0, iProteinsNamecolumnNum = 0;
	int iPrecursorChargeColumnNum = 0;
	vector<int> vecIntensitiesColumn;

	char Buffer[BUFFERLENGTH];
	string strLine;
	char *pstr;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	strLine = Buffer;
    GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns, ",");
	if (!mf_GetPeakViewAttributeColumns(mapAttrtibuteAndcolumns, iSequencecolumnNum, iProteinsNamecolumnNum, 
		vecIntensitiesColumn, iPrecursorChargeColumnNum))
	{
		flog.mf_Destroy();
		exit(1);
	}

	vector<string> vecStrTemps;
	string strPeptideSequenceTemp;
	string strProteinNameTemp;
	map<string, vector<string>>::iterator mapPeptidesAndProteinNamesIter;
	map<string, vector<string>>::iterator mapProteinAndPeptidesIter;
	vector<string> vecStrProteinNamesTemp;
	bool bIfAppearedBefore;

	fgets(Buffer, BUFFERLENGTH, pFile);
	while (!feof(pFile))
	{
		strLine = Buffer;
		if (strLine == "\n")
		{
			fgets(Buffer, BUFFERLENGTH, pFile);
			continue;
		}
		vecStrTemps.clear();
		vecStrTemps = split(strLine, ",");
		strPeptideSequenceTemp = vecStrTemps.at(iSequencecolumnNum);
		strProteinNameTemp = vecStrTemps.at(iProteinsNamecolumnNum);

		if (strProteinNameTemp == "[ RT-Cal protein ]")
		{
			fgets(Buffer, BUFFERLENGTH, pFile);
			continue;
		}
		// delete modification 
		int iBegin = 0, iEnd = 0;
		iBegin = strPeptideSequenceTemp.find("[");
		while (iBegin != strPeptideSequenceTemp.npos)
		{
			iEnd = strPeptideSequenceTemp.find("]", iBegin + 1);
			strPeptideSequenceTemp.erase(iBegin, iEnd - iBegin + 1);
			if (iEnd < strPeptideSequenceTemp.size() - 1)
			{
				iBegin = strPeptideSequenceTemp.find("[", iBegin);
			}
			else
			{
				iBegin = strPeptideSequenceTemp.npos;
			}
		}

		//mapPeptidesAndProteinNames
		mapPeptidesAndProteinNamesIter = mapPeptidesAndProteinNames.find(strPeptideSequenceTemp);
		if (mapPeptidesAndProteinNamesIter == mapPeptidesAndProteinNames.end())
		{ // new peptide
			vecStrProteinNamesTemp.clear();
			vecStrProteinNamesTemp.push_back(strProteinNameTemp);
			mapPeptidesAndProteinNames.insert(pair<string, vector<string>>(strPeptideSequenceTemp, vecStrProteinNamesTemp));
		}
		else
		{
			bIfAppearedBefore = false;
			for (int i = 0; i < mapPeptidesAndProteinNamesIter->second.size(); i++)
			{
				if (mapPeptidesAndProteinNamesIter->second.at(i) == strProteinNameTemp)
				{
					bIfAppearedBefore = true;
					break;
				}
			}
			if (!bIfAppearedBefore)
			{
				mapPeptidesAndProteinNamesIter->second.push_back(strProteinNameTemp);
			}
		}

		//mapProteinAndPeptides
		vector<string> vecPeptideSequencesTemp;
		mapProteinAndPeptidesIter = mapProteinAndPeptides.find(strProteinNameTemp);
		if (mapProteinAndPeptidesIter == mapProteinAndPeptides.end())
		{// new peptide
			vecPeptideSequencesTemp.clear();
			vecPeptideSequencesTemp.push_back(strPeptideSequenceTemp);
			mapProteinAndPeptides.insert(pair<string, vector<string>>(strProteinNameTemp, vecPeptideSequencesTemp));
		}
		else
		{
			bIfAppearedBefore = false;
			for (int i = 0; i < mapProteinAndPeptidesIter->second.size(); i++)
			{
				if (mapProteinAndPeptidesIter->second.at(i) == strPeptideSequenceTemp)
				{
					bIfAppearedBefore = true;
					break;
				}
			}
			if (!bIfAppearedBefore)
			{
				mapProteinAndPeptidesIter->second.push_back(strPeptideSequenceTemp);
			}

		}

		fgets(Buffer, BUFFERLENGTH, pFile);
	} // end of reading file

	fclose(pFile);

}

void CLoadPeakViewIO::ScreenSharedPeptides(const CLoadParam &param, map<string, vector<string>>& mapSharedPeptidesAndProteins,
	map<string, string> &mapGUniquePeptideAndProteinGroups, map<string, vector<string>>& mapProteinGroupAndGUniquePeptides)
{
	cout << "\tScreening shared peptides\n";
	flog.mf_Input("\tScreening shared peptides\n");

	map<string, vector<string>> mapPeptidesAndProteinNames;
	map<string, vector<string>>::iterator mapPeptidesAndProteinNamesIter;
	map<string, vector<string>> mapProteinAndPeptides;
	map<string, vector<string>>::iterator mapProteinAndPeptidesIter;
	mf_GetPeptidesAndProteins(param, mapPeptidesAndProteinNames, mapProteinAndPeptides);
	cout <<"\t\t" << mapPeptidesAndProteinNames.size() << " peptides loaded" << endl;
	cout << "\t\t" << mapProteinAndPeptides.size() << " proteins loaded" << endl;
	flog.mf_Input("\t\t" + fInt2String(mapPeptidesAndProteinNames.size()) + " peptides loaded\n");
	flog.mf_Input("\t\t" + fInt2String(mapProteinAndPeptides.size()) + " proteins loaded\n");

	vector<CRelatedProtein> vecSharedProteins;
	map<string, CRelatedProtein> mapDistinctProIDAndRelatedInfo;
	map<string, CRelatedProtein>::iterator mapDistinctProIDAndRelatedInfoIter;
	CRelatedProtein RelatedProeinTemp;

	//DistinctProteins
	mapPeptidesAndProteinNamesIter = mapPeptidesAndProteinNames.begin();
	for (; mapPeptidesAndProteinNamesIter != mapPeptidesAndProteinNames.end(); mapPeptidesAndProteinNamesIter++)
	{
		if (mapPeptidesAndProteinNamesIter->second.size() == 1)
		{
			mapDistinctProIDAndRelatedInfoIter = mapDistinctProIDAndRelatedInfo.find(mapPeptidesAndProteinNamesIter->second.at(0));
			if (mapDistinctProIDAndRelatedInfoIter == mapDistinctProIDAndRelatedInfo.end())
			{// new protein
				RelatedProeinTemp.Clear();
				RelatedProeinTemp.strProteinID = mapPeptidesAndProteinNamesIter->second.at(0);
				RelatedProeinTemp.vecPeptides.push_back(mapPeptidesAndProteinNamesIter->first);
				RelatedProeinTemp.vecBIfPeptideGroupUnique.push_back(true);
				mapDistinctProIDAndRelatedInfo.insert(pair<string, CRelatedProtein>(mapPeptidesAndProteinNamesIter->second.at(0), RelatedProeinTemp));
			}
			else // exist protein
			{
				mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.push_back(mapPeptidesAndProteinNamesIter->first);
				mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.push_back(true);
			}
		}
	}

	//vecSharedProteins
	bool bIfUniquePeptide;
	mapProteinAndPeptidesIter = mapProteinAndPeptides.begin();
	for (; mapProteinAndPeptidesIter != mapProteinAndPeptides.end(); mapProteinAndPeptidesIter++)
	{
		mapDistinctProIDAndRelatedInfoIter = mapDistinctProIDAndRelatedInfo.find(mapProteinAndPeptidesIter->first);
		if (mapDistinctProIDAndRelatedInfoIter == mapDistinctProIDAndRelatedInfo.end())
		{ // not distinct protein
			RelatedProeinTemp.Clear();
			RelatedProeinTemp.strProteinID = mapProteinAndPeptidesIter->first;
			for (int i = 0; i < mapProteinAndPeptidesIter->second.size(); i++)
			{
				RelatedProeinTemp.vecPeptides.push_back(mapProteinAndPeptidesIter->second.at(i));
				RelatedProeinTemp.vecBIfPeptideGroupUnique.push_back(false);
			}
			vecSharedProteins.push_back(RelatedProeinTemp);
		}
		else // complete non protein unique peptides for DistinctProteins
		{
			for (int i = 0; i < mapProteinAndPeptidesIter->second.size(); i++)
			{
				bIfUniquePeptide = false;
				for (int j = 0; j < mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.size(); j++)
				{
					if (mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(j) == mapProteinAndPeptidesIter->second.at(i))
					{
						bIfUniquePeptide = true;
						break;
					}
				}
				if (!bIfUniquePeptide)
				{
					mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.push_back(mapProteinAndPeptidesIter->second.at(i));
					mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.push_back(false);
				}
			}
		}
	}

	mapDistinctProIDAndRelatedInfoIter = mapDistinctProIDAndRelatedInfo.begin();
	for (; mapDistinctProIDAndRelatedInfoIter != mapDistinctProIDAndRelatedInfo.end(); mapDistinctProIDAndRelatedInfoIter++)
	{
		mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.push_back(mapDistinctProIDAndRelatedInfoIter->first);
		for (int i = 0; i < vecSharedProteins.size(); i++)
		{
			if (mapDistinctProIDAndRelatedInfoIter->second.mf_ifSubset(vecSharedProteins.at(i)))
			{
				mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.push_back(vecSharedProteins.at(i).strProteinID);
				vecSharedProteins.at(i).bIfSubSet = true;
			}

		}
	}

	mapDistinctProIDAndRelatedInfoIter = mapDistinctProIDAndRelatedInfo.begin();
	for (; mapDistinctProIDAndRelatedInfoIter != mapDistinctProIDAndRelatedInfo.end(); mapDistinctProIDAndRelatedInfoIter++)
	{
		if (mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.size() > 1)
		{
			for (int i = 0; i < mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.size(); i++)
			{
				if (mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.at(i) == false)
				{
					mapPeptidesAndProteinNamesIter = mapPeptidesAndProteinNames.find(mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i));
					if (mapPeptidesAndProteinNamesIter == mapPeptidesAndProteinNames.end())
					{
						cout << "Error:\tCannot find peptide " << mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i) << endl;
						flog.mf_Input("Error:\tCannot find peptide " + mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i) + "\n");
					}
					if (ifSubset(mapPeptidesAndProteinNamesIter->second, mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs))
					{
						mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.at(i) = true;
					}
				}
			}

		}
	}

	for (int i = 0; i < vecSharedProteins.size(); i++)
	{
		vecSharedProteins.at(i).vecSubsetProteinIDs.push_back(vecSharedProteins.at(i).strProteinID);
		for (int j = i + 1; j < vecSharedProteins.size(); j++)
		{
			if (!vecSharedProteins.at(i).bIfSubSet)
			{
				if (vecSharedProteins.at(i).mf_ifSubset(vecSharedProteins.at(j)))
				{
					vecSharedProteins.at(i).vecSubsetProteinIDs.push_back(vecSharedProteins.at(j).strProteinID);
					vecSharedProteins.at(j).bIfSubSet = true;
				}
			}
		}
	}

	bool bIfCantainGroupUniquePeptide;
	for (int i = 0; i < vecSharedProteins.size(); i++)
	{
		if (vecSharedProteins.at(i).bIfSubSet == false && vecSharedProteins.at(i).vecSubsetProteinIDs.size()>1)
		{
			bIfCantainGroupUniquePeptide = false;
			RelatedProeinTemp.Clear();
			RelatedProeinTemp.strProteinID = vecSharedProteins.at(i).strProteinID;
			for (int j = 0; j < vecSharedProteins.at(i).vecSubsetProteinIDs.size(); j++)
			{
				RelatedProeinTemp.vecSubsetProteinIDs.push_back(vecSharedProteins.at(i).vecSubsetProteinIDs.at(j));
			}
			for (int j = 0; j < vecSharedProteins.at(i).vecPeptides.size(); j++)
			{
				RelatedProeinTemp.vecPeptides.push_back(vecSharedProteins.at(i).vecPeptides.at(j));
				mapPeptidesAndProteinNamesIter = mapPeptidesAndProteinNames.find(vecSharedProteins.at(i).vecPeptides.at(j));
				if (mapPeptidesAndProteinNamesIter == mapPeptidesAndProteinNames.end())
				{
					cout << "Error:\tCannot find peptide " << vecSharedProteins.at(i).vecPeptides.at(j) << endl;
					flog.mf_Input("Error:\tCannot find peptide " + vecSharedProteins.at(i).vecPeptides.at(j) + "\n");
				}
				if (ifSubset(mapPeptidesAndProteinNamesIter->second, vecSharedProteins.at(i).vecSubsetProteinIDs))
				{
					RelatedProeinTemp.vecBIfPeptideGroupUnique.push_back(true);
					vecSharedProteins.at(i).vecBIfPeptideGroupUnique.at(j) = true;
					bIfCantainGroupUniquePeptide = true;
				}
				else
				{
					RelatedProeinTemp.vecBIfPeptideGroupUnique.push_back(false);
				}
			}
			if (bIfCantainGroupUniquePeptide)
			{
				mapDistinctProIDAndRelatedInfo.insert(pair<string, CRelatedProtein>(vecSharedProteins.at(i).strProteinID, RelatedProeinTemp));
			}
		}
	}

	vector<string>vecGUniquePeptidesTemp;
	vector<string>vecProteinsTemp;
	map<string, vector<string>>::iterator mapSharedPeptidesAndProteinsIter;

	mapSharedPeptidesAndProteins.clear();
	mapProteinGroupAndGUniquePeptides.clear();
	mapGUniquePeptideAndProteinGroups.clear();
	string strProteinGroupTemp;
	mapDistinctProIDAndRelatedInfoIter = mapDistinctProIDAndRelatedInfo.begin();
	for (; mapDistinctProIDAndRelatedInfoIter != mapDistinctProIDAndRelatedInfo.end(); mapDistinctProIDAndRelatedInfoIter++)
	{
		vecGUniquePeptidesTemp.clear();
		for (int i = 0; i < mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.size(); i++)
		{
			if (mapDistinctProIDAndRelatedInfoIter->second.vecBIfPeptideGroupUnique.at(i))
			{
				strProteinGroupTemp = "";
				size_t j = 0;
				for (j = 0; j < mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.size()-1; j++)
				{
					strProteinGroupTemp += mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.at(j) + ";";
				}
				strProteinGroupTemp += mapDistinctProIDAndRelatedInfoIter->second.vecSubsetProteinIDs.at(j);
				vecGUniquePeptidesTemp.push_back(mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i));
				mapGUniquePeptideAndProteinGroups.insert(pair<string, string>(mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i), strProteinGroupTemp));
			}
			else
			{
				mapSharedPeptidesAndProteinsIter = mapSharedPeptidesAndProteins.find(mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i));
				if (mapSharedPeptidesAndProteinsIter == mapSharedPeptidesAndProteins.end())
				{// new shared peptide
					vecProteinsTemp.clear();
					vecProteinsTemp.push_back(mapDistinctProIDAndRelatedInfoIter->first);
					mapSharedPeptidesAndProteins.insert(pair<string, vector<string>>(mapDistinctProIDAndRelatedInfoIter->second.vecPeptides.at(i), vecProteinsTemp));
				}
				else
				{ // exist shared peptide
					mapSharedPeptidesAndProteinsIter->second.push_back(mapDistinctProIDAndRelatedInfoIter->first);
				}
			}
		}
		mapProteinGroupAndGUniquePeptides.insert(pair<string, vector<string>>(mapDistinctProIDAndRelatedInfoIter->first, vecGUniquePeptidesTemp));

	}

	for (int i = 0; i < vecSharedProteins.size(); i++)
	{
		if (vecSharedProteins.at(i).bIfSubSet == false)
		{
			for (int j = 0; j < vecSharedProteins.at(i).vecBIfPeptideGroupUnique.size(); j++)
			{
				if (vecSharedProteins.at(i).vecBIfPeptideGroupUnique.at(j) == false)
				{
					mapSharedPeptidesAndProteinsIter = mapSharedPeptidesAndProteins.find(vecSharedProteins.at(i).vecPeptides.at(j));
					if (mapSharedPeptidesAndProteinsIter == mapSharedPeptidesAndProteins.end())
					{// new shared peptide
						vecProteinsTemp.clear();
						vecProteinsTemp.push_back(vecSharedProteins.at(i).strProteinID);
						mapSharedPeptidesAndProteins.insert(pair<string, vector<string>>(vecSharedProteins.at(i).vecPeptides.at(j), vecProteinsTemp));
					}
					else
					{ // exist shared peptide
						mapSharedPeptidesAndProteinsIter->second.push_back(vecSharedProteins.at(i).strProteinID);
					}
				}
			}
		}
	}

	string strProteinGroupsPath = param.m_strQuantiResultPath + "\\ProteinGroups.txt";
	mf_SaveProteinGroups(strProteinGroupsPath, mapDistinctProIDAndRelatedInfo);
	string strSharedPeptidesPath = param.m_strQuantiResultPath + "\\SharedPeptides.txt";
	mf_SaveSharedPeptides(strSharedPeptidesPath, mapSharedPeptidesAndProteins);

	cout << "\t\t" << mapProteinGroupAndGUniquePeptides.size() << " protein groups obtained\n";
	flog.mf_Input("\t\t" + fInt2String(mapProteinGroupAndGUniquePeptides.size()) + " protein groups obtained\n");
	flog.mf_Input("\tDone\n");
}


 //For simplify, assume there is only one experiment in the PeakView result file.
void CLoadPeakViewIO::load(const CLoadParam &param)
{
	cout << "Loading peptides and proteins from the SWATH 2.0 result\n";
	flog.mf_Input("Loading peptides and proteins from the SWATH 2.0 result\n");

	map<string, vector<string>> mapSharedPeptidesAndProteins;
	map<string, vector<string>>::iterator mapSharedPeptidesAndProteinsIter;

	map<string, vector<string>> mapProteinGroupAndGUniquePeptides;
	map<string, vector<string>>::iterator mapProteinGroupAndGUniquePeptidesIter;
	map<string, string> mapGUniquePeptideAndProteinGroups;
	map<string, string>::iterator mapGUniquePeptideAndProteinGroupsIter;
	ScreenSharedPeptides(param, mapSharedPeptidesAndProteins, mapGUniquePeptideAndProteinGroups, mapProteinGroupAndGUniquePeptides);


	// open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, param.m_strIdentResultPath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << param.m_strIdentResultPath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strIdentResultPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}

	char Buffer[BUFFERLENGTH];
	string strLine;
	char *pstr;
	//get the columns of sequence、ID、Proteins、Intensity UPS2_yeast by analysising the first row;
	int iSequencecolumnNum = 0, iProteinsNamecolumnNum = 0;
	int iPrecursorChargeColumnNum = 0;

	vector<int> veciIntensitycolumnNum;
	int icolumns = 0, count = 0;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	stringstream ss;
	// get the columes of peptides attributes according to the first row
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	strLine = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns, ",");
	if (!mf_GetPeakViewAttributeColumns(mapAttrtibuteAndcolumns, iSequencecolumnNum, iProteinsNamecolumnNum,
		veciIntensitycolumnNum,iPrecursorChargeColumnNum))
	{
		flog.mf_Destroy();
		exit(1);
	}

	map<string, PeptideQuantInfo> mapPeptidesAndQuantInfo;
	map<string, PeptideQuantInfo>::iterator mapPeptidesAndQuantInfoIter;
	PeptideQuantInfo PeptideQuantInfoTemp;
	string strPeptideSequenceTemp;
	string strPrecursorChargeTemp;
	string strProteinNameTemp;
	string strFragmentType;
	double dPeptideIntensityTemp;

	CProtein cprotein;
	vector<string> vecStrTemps;
	string strSharedPeptideSequence;
	vector<double> vecPeptideIntensitiesInRep;
	
	vector<double> vecMS2IntensitiesTemp;
	vector<string>vecProteinsTemp;
	while (!feof(pFile))
	{
		fgets(Buffer, BUFFERLENGTH, pFile);
		strLine = Buffer;
		if (strLine == "\n")
		{
			continue;
		}
		vecStrTemps.clear();
		vecStrTemps = split(strLine, ",");
		if (vecStrTemps.size() != mapAttrtibuteAndcolumns.size())
		{
			flog.mf_Input("Warning:\tThe columns number of line "+strLine+" do not equal to the columns number of the header line\n");
		}
		strPeptideSequenceTemp = vecStrTemps.at(iSequencecolumnNum);
		strPrecursorChargeTemp = vecStrTemps.at(iPrecursorChargeColumnNum);
		strProteinNameTemp = vecStrTemps.at(iProteinsNamecolumnNum);
		if (param.m_bIfExistDecoyProtein)
		{
			if (strProteinNameTemp.find(param.m_strDecoyProteinIDPrefix) != strProteinNameTemp.npos)
			{
				continue;
			}
		}
		vecPeptideIntensitiesInRep.clear();

		// delete modification 
		int iBegin = 0, iEnd = 0;
		iBegin = strPeptideSequenceTemp.find("[");
		while (iBegin != strPeptideSequenceTemp.npos)
		{
			iEnd = strPeptideSequenceTemp.find("]", iBegin + 1);
			strPeptideSequenceTemp.erase(iBegin, iEnd-iBegin+1);
			if (iEnd < strPeptideSequenceTemp.size() - 1)
			{
				iBegin = strPeptideSequenceTemp.find("[", iBegin);
			}
			else
			{
				iBegin = strPeptideSequenceTemp.npos;
			}		
		}

		for (int i = 0; i < veciIntensitycolumnNum.size(); i++)
		{
			dPeptideIntensityTemp = atof(vecStrTemps.at(veciIntensitycolumnNum.at(i)).c_str());
			vecPeptideIntensitiesInRep.push_back(dPeptideIntensityTemp);
		}
		dPeptideIntensityTemp = Average(vecPeptideIntensitiesInRep);

		vecMS2IntensitiesTemp.clear();
		mapGUniquePeptideAndProteinGroupsIter = mapGUniquePeptideAndProteinGroups.find(strPeptideSequenceTemp);
		if (mapGUniquePeptideAndProteinGroupsIter != mapGUniquePeptideAndProteinGroups.end())
		{
			vecProteinsTemp = split(mapGUniquePeptideAndProteinGroupsIter->second, ";");
		}
		if (mapGUniquePeptideAndProteinGroupsIter != mapGUniquePeptideAndProteinGroups.end() && vecProteinsTemp[0] == strProteinNameTemp)
		{

			mapPeptidesAndQuantInfoIter = mapPeptidesAndQuantInfo.find(strPeptideSequenceTemp);
			if (mapPeptidesAndQuantInfoIter == mapPeptidesAndQuantInfo.end())
			{// new peptide
				PeptideQuantInfoTemp.Clear();
				PeptideQuantInfoTemp.strProteinID = strProteinNameTemp;
				PeptideQuantInfoTemp.strProteinGroups = mapGUniquePeptideAndProteinGroupsIter->second;
				PeptideQuantInfoTemp.vecPrecursorCharges.push_back(strPeptideSequenceTemp);
				PeptideQuantInfoTemp.vecdPeptideIntensities.push_back(dPeptideIntensityTemp);
				mapPeptidesAndQuantInfo.insert(pair<string, PeptideQuantInfo>(strPeptideSequenceTemp, PeptideQuantInfoTemp));
			}
			else
			{ // old peptide, new precursor charge
				mapPeptidesAndQuantInfoIter->second.vecPrecursorCharges.push_back(strPeptideSequenceTemp);
				mapPeptidesAndQuantInfoIter->second.vecdPeptideIntensities.push_back(dPeptideIntensityTemp);
			}
		}
	}
	fclose(pFile);

	vector<double> vecdIntensityTemp;
	map<string, int> mapProteinsIndexs;
	map<string, int>::iterator mapProteinsIndexsIter;
	int iProteinIndex = 0;
	mapPeptidesAndQuantInfoIter = mapPeptidesAndQuantInfo.begin();
	for (; mapPeptidesAndQuantInfoIter != mapPeptidesAndQuantInfo.end(); mapPeptidesAndQuantInfoIter++)
	{
		dPeptideIntensityTemp = 0.0;
		for (int i = 0; i<mapPeptidesAndQuantInfoIter->second.vecdPeptideIntensities.size(); i++)
		{
			dPeptideIntensityTemp += mapPeptidesAndQuantInfoIter->second.vecdPeptideIntensities.at(i);
		}
		mapProteinsIndexsIter = mapProteinsIndexs.find(mapPeptidesAndQuantInfoIter->second.strProteinID);
		if (mapProteinsIndexsIter == mapProteinsIndexs.end())
		{// new protein
			mapGUniquePeptideAndProteinGroupsIter = mapGUniquePeptideAndProteinGroups.find(mapPeptidesAndQuantInfoIter->first);
			if (mapGUniquePeptideAndProteinGroupsIter == mapGUniquePeptideAndProteinGroups.end())
			{
				flog.mf_Input("Error:\tCannot find peptide " + mapPeptidesAndQuantInfoIter->first + " in unique peptides.\n");
				flog.mf_Destroy();
				exit(1);
			}
			cprotein.Clear();
			vecdIntensityTemp.clear();
			vecdIntensityTemp.push_back(dPeptideIntensityTemp);
			cprotein.m_mapExperimentAndPeptidesIntensity.insert(\
                pair<string, vector<double>>("AssumeOnlyOneExperiment", vecdIntensityTemp));
			cprotein.m_strProteinID = mapPeptidesAndQuantInfoIter->second.strProteinID;
			cprotein.m_strProteinFullName = mapGUniquePeptideAndProteinGroupsIter->second;
			cprotein.m_vPeptidesSequences.push_back(mapPeptidesAndQuantInfoIter->first);
            m_vecProteins.push_back(cprotein);
			mapProteinsIndexs.insert(pair<string, int>(\
                mapPeptidesAndQuantInfoIter->second.strProteinID, iProteinIndex));
			iProteinIndex++;

		}
		else
		{ // old protein
			m_vecProteins.at(mapProteinsIndexsIter->second).m_vPeptidesSequences.push_back(mapPeptidesAndQuantInfoIter->first);
			m_vecProteins.at(mapProteinsIndexsIter->second).m_mapExperimentAndPeptidesIntensity["AssumeOnlyOneExperiment"].push_back(dPeptideIntensityTemp);
		}
	}
    
	flog.mf_Input("Done\n");
}

void CLoadPeakViewIO::mf_saveProteins(CLoadParam param)
{
	ofstream ofile(param.m_strExtractedProteinsPath.c_str());
	if (!ofile)
	{
		cout << "Error:\tCannot open " << param.m_strExtractedProteinsPath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strExtractedProteinsPath + "\n");
		flog.mf_Destroy();
		exit(1);
	}
	else
	{
		cout << "Proteins and peptides saved to " << param.m_strExtractedProteinsPath << endl;
		flog.mf_Input("Proteins and peptides saved to " + param.m_strExtractedProteinsPath + "\n");

	}
	ofile << "ProteinID\tMajority protein IDs\tPeptidesNumber\tvPeptidesSequence\t";
	ofile << "Intensity ExperimentOnlyOne\tIntensity\t";
	ofile << "bIfiBAQExist\t";
	ofile << "\n";

	vector<CProtein>::iterator  ProteinIter;
	vector<string>::iterator vPeptideIter;
	vector<double>::iterator vPeptideIntensityIter;
	vector<double> dvCombinedPepIntensities;

	int i = 0;
	map<string, vector<double>>::iterator mapExperimentAndPeptideIntensitysIter;
	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
	{
		ofile << ProteinIter->m_strProteinID << "\t" << ProteinIter->m_strProteinFullName << "\t";
		ofile << ProteinIter->m_vPeptidesSequences.size() << "\t";
		dvCombinedPepIntensities.clear();

		for (vPeptideIter = ProteinIter->m_vPeptidesSequences.begin(); vPeptideIter != ProteinIter->m_vPeptidesSequences.end(); vPeptideIter++)
		{
			ofile << *vPeptideIter << ";";
			dvCombinedPepIntensities.push_back(0.0);
		}
			
		ofile << "\t";

		mapExperimentAndPeptideIntensitysIter = ProteinIter->m_mapExperimentAndPeptidesIntensity.begin();
		for (; mapExperimentAndPeptideIntensitysIter != ProteinIter->m_mapExperimentAndPeptidesIntensity.end(); mapExperimentAndPeptideIntensitysIter++)
		{
			for (i = 0; i < mapExperimentAndPeptideIntensitysIter->second.size(); i++)
			{
                ofile << fixed << std::setprecision(6) << mapExperimentAndPeptideIntensitysIter->second.at(i) << ";";
				dvCombinedPepIntensities[i] += mapExperimentAndPeptideIntensitysIter->second.at(i);
			}
			ofile << "\t";
		}
		mapExperimentAndPeptideIntensitysIter = ProteinIter->m_mapExperimentAndPeptidesIntensity.begin();
		for (i = 0; i < mapExperimentAndPeptideIntensitysIter->second.size(); i++)
		{
			ofile << fixed << std::setprecision(6) << dvCombinedPepIntensities[i] << ";";
		}

		ofile << "\t";
		ofile << 0 << "\t";
		ofile << "\n";

	}
	ofile.close();

	// Save necessary experiment information
	ofstream fout(param.m_strQuantiResultPath+ "\\PeakViewExperiment.config");
	fout << "Experiment Number\t" << 1 << endl;
	fout << "Experiment id\tExperimentOnlyOne" << endl;
	fout.close();
}

