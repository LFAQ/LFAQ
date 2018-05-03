#include"stdafx.h"
#include"LoadpQuantProcess.h"
#include"loadParam.h"
#include"DataIO.h"
#include"ProteinDigestion.h"

//���㵰�׵�IBAQֵ
bool CProteinWorker::mf_CalculateProteinIntensity(Cpeptides peptides, map<string, string> &mapMsMsIDAndPeptidesSequenceWithMod, CLoadParam param)
{
	cout << "\tLoading Protein!\n\n";

	//���ļ�
	FILE * pinFile;
	pinFile = fopen(param.m_strProteinFilePath.c_str(), "r");
	if (pinFile == NULL)
	{
		cout << "Error when opening \"" << param.m_strProteinFilePath << endl;
		return false;
	}
	ofstream ofile(param.m_strProteinsPath.c_str());
	if (!ofile)
	{
		cout << "Error when opening \"" << param.m_strProteinsPath << endl;
		return false;
	}
	ofile << "m_strProteinID\tm_dPeptidesIntensitySum \tm_iNumberOfTheoreticEnzyme\tcprotein.m_dIBAQ\n";
	CLoadPfindIO loadpfindio;
	map<string, string> mapProteinIDAndSequenceTemp;

	//int iPeptidesNumberTemp;
	// load protein fasta
	loadpfindio.mf_LoadProteinFasta(mapProteinIDAndSequenceTemp, param.m_strProteinFasterFilePath);

	//׼���������õ�����ʱ����
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	fgets(Buffer, BUFFERLENGTH, pinFile);
	string strProteinIDTemp;
	string strProteinIDPrefixTemp;
	string strPeptideSequenceTemp;

	CProtein cprotein;
	vector<string>::iterator peptideSequenceIter;
	map<string, string>::iterator ProteinIDSequenceIter;
	string  strWithFlankingRegionTemp;
	CPeptideAttribute peptideAttributeTemp;
	PepAttributeWorker PeptideAttributeWorkerTemp;
	map<string, double>::iterator mapPeptideSequenceAndIntensityIter;
	CProteinDigestion proteinDigestion(param);
	map<string, string>::iterator MsMsIdAndPeptideSequenceWithModIter;
	string strMsMsIdTemp;

	while (!feof(pinFile))
	{
		pstr = Buffer;
		if (pstr[0] == '@')
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				//��ȡ��������
				pstr = pstr1 + 1;
				pstr1 = strstr(pstr, "\t");
				if (pstr1 == NULL)
				{
					cout << "Wrong Protein ID type!" << endl;
					exit(0);
				}
				*pstr1 = '\0';
				strProteinIDTemp = pstr;
				strProteinIDPrefixTemp = strProteinIDTemp.substr(0, param.m_strReverseProteinIDPrefix.size());
				transform(strProteinIDPrefixTemp.begin(), strProteinIDPrefixTemp.end(), strProteinIDPrefixTemp.begin(), ::tolower);
				// delete Reverse protein and contaminant protein
				if (strProteinIDPrefixTemp == param.m_strReverseProteinIDPrefix || strProteinIDPrefixTemp == param.m_strContaminantProteinIDPrefix)
				{
					fgets(Buffer, BUFFERLENGTH, pinFile);
					continue;
				}

				ProteinIDSequenceIter = mapProteinIDAndSequenceTemp.find(strProteinIDTemp);
				if (ProteinIDSequenceIter != mapProteinIDAndSequenceTemp.end())
				{
					cprotein.m_strProteinSequence = ProteinIDSequenceIter->second;
				}
				else
				{
					cout << "��fasta�ļ����Ҳ���������" << strProteinIDTemp << endl;
					exit(0);
				}


				//Ĭ������4��֮��������Ķ����е���
				for (int i = 0; i < 4; i++)
				{
					fgets(Buffer, BUFFERLENGTH, pinFile);
				}

				fgets(Buffer, BUFFERLENGTH, pinFile);
				pstr = Buffer;
				while ((pstr[0] != '@') && (!feof(pinFile)))
				{
					pstr1 = strstr(pstr, "\t");
					if (pstr1 != NULL)
					{
						*pstr1 = '\0';
						strMsMsIdTemp = pstr;
						MsMsIdAndPeptideSequenceWithModIter = mapMsMsIDAndPeptidesSequenceWithMod.find(strMsMsIdTemp);
						if (MsMsIdAndPeptideSequenceWithModIter == mapMsMsIDAndPeptidesSequenceWithMod.end())
						{
							cerr << "Can not find msms " << strMsMsIdTemp << endl;
							exit(0);
						}
						else
						{
							strPeptideSequenceTemp = MsMsIdAndPeptideSequenceWithModIter->second;
						}
						// enzyme number
						cprotein.m_iNumberOfTheoreticEnzyme = proteinDigestion.mf_CalculateOptDigestionNumber(cprotein.m_strProteinSequence);
						////�����Ķ�����
						//strWithFlankingRegionTemp = cprotein.mf_GetPeptidesAdjacentSequence(strPeptideSequenceTemp);
						//cprotein.m_vPeptidesSequencesWithFlankingRegion.push_back(strWithFlankingRegionTemp);                                         //�����Ķ�IDȡ����Ӧ�Ķ����к�����15λ������					
						//peptideAttributeTemp.Clear();
						//peptideAttributeTemp.m_vecAttributes = PeptideAttributeWorkerTemp.mf_GetAttributeFromSequence(strPeptideSequenceTemp, strWithFlankingRegionTemp);
						//cprotein.m_vPeptidesAttributes.push_back(peptideAttributeTemp);

						//Ѱ����Ӧ�Ķε�intgensity��
						mapPeptideSequenceAndIntensityIter = peptides.m_mapPeptideSequenceAndIntensity.find(strPeptideSequenceTemp);
						if (mapPeptideSequenceAndIntensityIter == peptides.m_mapPeptideSequenceAndIntensity.end())
						{
							cerr << "Can not find peptide " << strPeptideSequenceTemp << endl;
							exit(0);
						}
						else
						{
							cprotein.m_vPeptidesNativeIntensity.push_back(mapPeptideSequenceAndIntensityIter->second);
							cprotein.m_vPeptidesIntensity.push_back(mapPeptideSequenceAndIntensityIter->second);
						}
						//iPeptidesNumberTemp++;

					}
				fgets(Buffer, BUFFERLENGTH, pinFile);
				pstr = Buffer;
			}//end while

		}// end if 

		//pQuant�Ľ����û�и���LFQ��intensityֵ��
		cprotein.m_dLFQ = 0.0;
		cprotein.m_strProteinID = strProteinIDTemp;
		cprotein.m_strProteinFullName = strProteinIDTemp;

		mf_CalculateProteinIBAQ_Intensity(cprotein);

		// save the result 

			ofile << cprotein.m_strProteinID << "\t" << cprotein.m_dPeptidesIntensitySum << "\t" << cprotein.m_iNumberOfTheoreticEnzyme << "\t" << cprotein.m_dIBAQ << endl;

		cprotein.Clear();

	} // end if(pstr[0] == '@')  һ������ 
	if ((Buffer[0] != '@') && (!feof(pinFile)))
		fgets(Buffer, BUFFERLENGTH, pinFile);
}
	fclose(pinFile);
	ofile.close();
	return 1;
}

void  CProteinWorker::mf_CalculateProteinIBAQ_Intensity(CProtein &cprotein)
{

	for (size_t i = 0; i < cprotein.m_vPeptidesNativeIntensity.size(); i++)
		{
			cprotein.m_dPeptidesIntensitySum += cprotein.m_vPeptidesNativeIntensity.at(i);
		}
	cprotein.m_iPeptidesNumber = cprotein.m_vPeptidesSequences.size();
	if (cprotein.m_iNumberOfTheoreticEnzyme>0)
		cprotein.m_dIBAQ = cprotein.m_dPeptidesIntensitySum / cprotein.m_iNumberOfTheoreticEnzyme;
	else
	{
		cerr << "protein " << cprotein.m_strProteinID << " do not have Digension peptides." << endl;
		exit(0);
	}
}

