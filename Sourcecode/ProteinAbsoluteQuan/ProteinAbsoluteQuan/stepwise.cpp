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
#include"stepwise.h"


void Cstepwise::mf_StepwiseInit( vector<CProtein> proteins)
{
	cout << "\tStepwise Initing! \n";
	flog.mf_Input("\tStepwise Initing! \n");
	int i;
	m_iNumberOfFeature = NumberOfFeature;
	vector<CProtein>::iterator ProteinIter;
	m_NumberOfTrainSamples = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		m_NumberOfTrainSamples += ProteinIter->m_vPeptidesSequences.size();
	}
	m_px = new double*[m_NumberOfTrainSamples];         // dynamic allocation memory
	for (i = 0; i < m_NumberOfTrainSamples; i++) { m_px[i] = new double[m_iNumberOfFeature + 1]; }
	if (m_px == NULL)
	{
		cout << "Error:\tmemory allocation failure for m_px" << endl;
		flog.mf_Input("Error:\tmemory allocation failure for m_px.\n");
	}
	m_pTestXY = new double*[m_NumberOfTrainSamples];         // dynamic allocation memory
	for (i = 0; i < m_NumberOfTrainSamples; i++) { m_pTestXY[i] = new double[m_iNumberOfFeature + 1]; }
	m_pTrainXForCV = new double*[m_NumberOfTrainSamples];         // dynamic allocation memory
	for (i = 0; i < m_NumberOfTrainSamples; i++) { m_pTrainXForCV[i] = new double[m_iNumberOfFeature + 1]; }
	m_pTrainPredictY = new double[m_NumberOfTrainSamples];
	m_pTestPredictY = new double[m_NumberOfTrainSamples];
	m_pNativeTestPredictY = new double[m_NumberOfTrainSamples]; //added by CC 20150527
	m_psortIndex = new int[m_NumberOfTrainSamples];
	m_pLeaveXIndex = new int[m_NumberOfTrainSamples];
	m_pr = new double*[m_iNumberOfFeature + 1];
	for (i = 0; i < m_iNumberOfFeature + 1; i++) { m_pr[i] = new double[m_iNumberOfFeature + 1]; }
	m_pxx = new double[m_iNumberOfFeature + 1];
	m_pCoefficient = new double[m_iNumberOfFeature + 1];
	m_pv = new double[m_iNumberOfFeature + 1];
	m_ps = new double[m_iNumberOfFeature + 1];
	m_pye = new double[m_NumberOfTrainSamples];
	m_pyr = new double[m_NumberOfTrainSamples];
}

bool Cstepwise::mf_Normalize(int i)
{
	if (i == 0)
		return mf_Normalize(m_px, m_NumberOfLimitedSamples, m_iNumberOfFeature);
	else if (i == 1)
		return mf_Normalize(m_pTestXY, m_NumberOfTrainSamples, m_iNumberOfFeature);
	else
		return false;
}

//using ZScore； 
bool Cstepwise::mf_Normalize(double **p, int m, int n)          
{
	if (n < 1 || m < 3) 
		return false;

	double *meanXY = new double[n + 1];
	double *sdeXY = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
		meanXY[i] = 0.0;
	for (int i = 0; i < m; i++)   
		for (int j = 0; j < n + 1; j++)
			meanXY[j] += p[i][j];
	for (int j = 0; j < n + 1; j++)
		meanXY[j] = meanXY[j] / m;

	for (int i = 0; i < n + 1; i++)
		sdeXY[i] = 0.0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n + 1; j++)
			sdeXY[j] += (p[i][j] - meanXY[j])*(p[i][j] - meanXY[j]);
	for (int j = 0; j < n + 1; j++)
		sdeXY[j] = sqrt(sdeXY[j] / (m - 1));

	m_dMeanY = meanXY[n];       
	m_dSdeY = sdeXY[n];
	//standarization
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n + 1; j++)
			p[i][j] = (p[i][j] - meanXY[j]) / sdeXY[j];
	delete[]meanXY;
	delete[] sdeXY;
	return true;
}


void Cstepwise::mf_StepwiseRegression( double ff1, double ff2, double es)
{
	cout << "\tstepwise begin:\n";
	flog.mf_Input("\tstepwise begin:\n");
	int i, j, ii, m, imi, imx, l, it;
	double z, phi, sd, vmi, vmx, q, fmi, fmx;
	m_dthresh1 = ff1;  m_dthresh2 = ff2;  m_dEps = es;
	m = m_iNumberOfFeature + 1; q = 0.0;
	double *Tempsd = new double[m_iNumberOfFeature + 1];
	for (j = 0; j <= m_iNumberOfFeature; j++)
	{
		z = 0.0;
		for (i = 0; i <= m_NumberOfLimitedSamples - 1; i++)
			z = z + m_px[i][j];                   
		m_pxx[j] = z / m_NumberOfLimitedSamples;
	}
	for (i = 0; i <= m_iNumberOfFeature; i++)
		for (j = 0; j <= i; j++)
		{
			z = 0.0;
			for (ii = 0; ii <= m_NumberOfLimitedSamples - 1; ii++)
				z = z + (m_px[ii][i] - m_pxx[i])*(m_px[ii][j] - m_pxx[j]);
			m_pr[i][j] = z;
		}
	for (i = 0; i <= m_iNumberOfFeature; i++)  Tempsd[i] = sqrt(m_pr[i][i]);// the standard deviation of ith feature
	for (i = 0; i <= m_iNumberOfFeature; i++)
		for (j = 0; j <= i; j++)
		{
			m_pr[i][j] = m_pr[i][j] / (Tempsd[i] * Tempsd[j]);
			m_pr[j][i] = m_pr[i][j];
		}
	phi = m_NumberOfLimitedSamples - 1.0;
	sd = Tempsd[m_iNumberOfFeature] / sqrt(m_NumberOfLimitedSamples - 1.0);  //the standard devivation of last feature;
	it = 1;
	while (it == 1)
	{
		it = 0;
		vmi = 1.0e+35; vmx = 0.0;
		imi = -1; imx = -1;
		for (i = 0; i <= m_iNumberOfFeature; i++)
		{
			m_pv[i] = 0.0;
			m_pCoefficient[i] = 0.0; 
			m_ps[i] = 0.0;
		}
		for (i = 0; i <= m_iNumberOfFeature - 1; i++)
			if (m_pr[i][i] >= m_dEps)
			{
				m_pv[i] = m_pr[i][m_iNumberOfFeature] * m_pr[m_iNumberOfFeature][i] / m_pr[i][i];
				if (m_pv[i] >= 0.0)
				{
					if (m_pv[i]>vmx) { vmx = m_pv[i]; imx = i; }
				}
				else
				{
					m_pCoefficient[i] = m_pr[i][m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature] / Tempsd[i];
					m_ps[i] = sqrt(m_pr[i][i])*sd / Tempsd[i];
					if (fabs(m_pv[i])<vmi)
					{
						vmi = fabs(m_pv[i]); imi = i;
					}
				}
			}
		if (phi != m_iNumberOfFeature - 1.0)
		{
			z = 0.0;
			for (i = 0; i <= m_iNumberOfFeature - 1; i++)  
				z = z + m_pCoefficient[i] * m_pxx[i];
			m_pCoefficient[m_iNumberOfFeature] = m_pxx[m_iNumberOfFeature] - z;
			m_ps[m_iNumberOfFeature] = sd; 
			m_pv[m_iNumberOfFeature] = q;
		}
		else
		{
			m_pCoefficient[m_iNumberOfFeature] = m_pxx[m_iNumberOfFeature]; m_ps[m_iNumberOfFeature] = sd;
		}
		fmi = vmi*phi / m_pr[m_iNumberOfFeature][m_iNumberOfFeature];
		fmx = (phi - 1.0)*vmx / (m_pr[m_iNumberOfFeature][m_iNumberOfFeature] - vmx);
		if ((fmi<m_dthresh2) || (fmx >= m_dthresh1))
		{
			if (fmi<m_dthresh2)  
			{ 
				phi = phi + 1.0;
				l = imi; 
			    cout << "\t\tDelete variable " << l << "\n"; 
				flog.mf_Input("\t\tDelete variable " + fInt2String(l)+"\n");
			}
			else 
			{ 
				phi = phi - 1.0; 
				l = imx; 
				cout << "\t\tAdd variable " << l << "\n";
				flog.mf_Input("\t\tAdd variable " +fInt2String(l)+"\n");
			}
			for (i = 0; i <= m_iNumberOfFeature; i++)
				if (i != l)
					for (j = 0; j <= m_iNumberOfFeature; j++)
						if (j != l)
							m_pr[i][j] = m_pr[i][j] - (m_pr[l][j] / m_pr[l][l])*m_pr[i][l];
			for (j = 0; j <= m_iNumberOfFeature; j++)
				if (j != l) m_pr[l][j] = m_pr[l][j] / m_pr[l][l];
			for (i = 0; i <= m_iNumberOfFeature; i++)
				if (i != l)  m_pr[i][l] = -m_pr[i][l] / m_pr[l][l];
			m_pr[l][l] = 1.0 / m_pr[l][l];
			q = m_pr[m_iNumberOfFeature][m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature];
			sd = sqrt(m_pr[m_iNumberOfFeature][m_iNumberOfFeature] / phi)*Tempsd[m_iNumberOfFeature];
			m_dC = sqrt(1.0 - m_pr[m_iNumberOfFeature][m_iNumberOfFeature]);
			m_dF = (phi*(1.0 - m_pr[m_iNumberOfFeature][m_iNumberOfFeature])) / ((m_NumberOfLimitedSamples - phi - 1.0)*m_pr[m_iNumberOfFeature][m_iNumberOfFeature]);
			it = 1;
		}
	}
	for (i = 0; i <= m_NumberOfLimitedSamples - 1; i++)
	{
		z = 0.0;
		for (j = 0; j <= m_iNumberOfFeature - 1; j++) z = z + m_pCoefficient[j] * m_px[i][j];
		m_pye[i] = m_pCoefficient[m_iNumberOfFeature] + z; m_pyr[i] = m_px[i][m_iNumberOfFeature] - m_pye[i];
	}
	delete[]Tempsd;
}

void Cstepwise::mf_saveStepwiseXY(string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "The training set saved to" << path << endl;
	flog.mf_Input("The training set saved to " + path +"\n");
	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;
	for (int i = 0; i<m_NumberOfLimitedSamples; i++)
	{
		for (int j = 0; j <= m_iNumberOfFeature; j++)
			ofile << m_px[i][j] << "\t";
		ofile << "\n";
	}
	ofile.close();

}

void Cstepwise::mf_saveRegressionResult(string path) 
{
	int i, j;

	ofstream fout(path.c_str(), ios::out);
	if (!fout)
	{
		cout << "Error:\tCannot open " << path << endl;
		flog.mf_Input("Error:\tCannot open " + path +"\n");
		flog.mf_Destroy();
		exit(1);
	}
	fout << endl;
	fout << "thresh1 = " << m_dthresh1 << "    " << "thresh2 = " << m_dthresh2 << endl;
	fout << "Eps = " << m_dEps << endl;
	fout << "square of deviance and residual sum of squares:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << "Ex(" << i << ")=" << m_pxx[i] << "    ";
	}
	fout << "\nEy=" << m_pxx[m_iNumberOfFeature] << endl;
	fout << "regression coefficients:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << m_pCoefficient[i] << "\t";
	}
	fout << m_pCoefficient[i];
	fout << endl << endl;
	fout << "sum of squares of partial regression:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << m_pv[i] << "\t";
	}
	fout << endl << "residual sum of squares = " << m_pv[m_iNumberOfFeature] << endl << endl;
	fout << "the standard deviation of regression coefficients:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << m_ps[i] << "\t";
	}
	fout << endl << "the estimated SD = " << m_ps[m_iNumberOfFeature]  << endl;
	fout << "multiple correlation coefficientsC = " << m_dC << endl << endl;
	fout << "the value of F-test = " << m_dF << endl << endl;
	fout << "estimated values of the conditional expectation of independent variables and ";
	fout<<" residuals of observed variables:" << endl;
	for (i = 0; i<m_NumberOfLimitedSamples; i++)
	{
		fout << m_pye[i] << "\t" << m_pyr[i] << endl;;
	}
	fout << endl;
	fout << "m_NumberOfLimitedSamples residuals of observed variables:" << endl;
	for (i = 0; i <= m_iNumberOfFeature; i++)
	{
		for (j = 0; j <= m_iNumberOfFeature; j++)
		{
			fout << m_pr[i][j] << " ";
		}
		fout << endl;
	}

	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
	{
		m_pTrainPredictY[i] = 0.0;
		for (int j = 0; j < m_iNumberOfFeature; j++)
			m_pTrainPredictY[i] += m_pCoefficient[j] * m_px[i][j];

		m_pTrainPredictY[i] += m_pCoefficient[m_iNumberOfFeature];
	}

	m_dRMSE = 0.0;
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
		m_dRMSE += (m_px[i][m_iNumberOfFeature] - m_pTrainPredictY[i])*(m_px[i][m_iNumberOfFeature] - m_pTrainPredictY[i]);
	m_dRMSE = sqrt(m_dRMSE / m_NumberOfLimitedSamples);
	fout << "The RMSE of stepwise regression is " << m_dRMSE << endl;

	m_dCorrelation = 0.0;
	double meanReal = 0.0;
	double meanTest = 0.0;
	double sd1 = 0.0, sd2 = 0.0;

	meanTest = Average(m_pTrainPredictY, m_NumberOfLimitedSamples);
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
		meanReal += m_px[i][m_iNumberOfFeature];
	meanReal = meanReal / m_NumberOfLimitedSamples;
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
		m_dCorrelation += (m_px[i][m_iNumberOfFeature] - meanReal)*(m_pTrainPredictY[i] - meanTest);
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
	{
		sd1 += (m_px[i][m_iNumberOfFeature] - meanReal)*(m_px[i][m_iNumberOfFeature] - meanReal);
		sd2 += (m_pTrainPredictY[i] - meanTest)*(m_pTrainPredictY[i] - meanTest);
	}
	m_dCorrelation = m_dCorrelation / sqrt(sd1*sd2);
	fout << "The Correlation between predicted Y and real Y  is " << m_dCorrelation << endl;
	cout << "\tThe Correlation between predicted Y and real Y  is " << m_dCorrelation << endl;
	flog.mf_Input("\tThe Correlation between predicted Y and real Y is " +fDouble2String(m_dCorrelation)+"\n");
	fout.close();

	cout << "\tThe regression result saved to " << path << endl;
	flog.mf_Input("\tThe regression result saved to " + path + "\n");

}
void  Cstepwise::mf_SaveFeatures( string strFeaturePath)
{
	ofstream oFile(strFeaturePath);
	if (!oFile)
	{
		flog.mf_Input("Error\tCannot open " + strFeaturePath +"\n");
		cout << "Error:\tCannot open " << strFeaturePath << endl;
	}

	oFile << "Feature Index\t Coefficience\n";
	for (int i = 1; i<m_iNumberOfFeature; i++)
	{
		oFile <<i<<"\t"<< m_pCoefficient[i] <<endl;
	}
	oFile.close();

}
bool Cstepwise::mf_ConstructXY(vector<CProtein> proteins)
{
	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	m_NumberOfLimitedSamples = 0;  
	m_iNumberOfTrainProteins = 0;

	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_bIfInTrainSet == true)
		{
			m_iNumberOfTrainProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.size(); j++)
				{
					m_px[m_NumberOfLimitedSamples][j] = ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.at(j);
				}

				m_px[m_NumberOfLimitedSamples][j] = ProteinIter->m_vPeptideExperimentalQFactors.at(n);
				m_NumberOfLimitedSamples++;
			}

		}

	}
	cout << "\tConstructed the training set using " << m_NumberOfLimitedSamples << " peptides in "<<m_iNumberOfTrainProteins <<" proteins"<< endl; 
	flog.mf_Input("\tConstructed the training set using " +fInt2String(m_NumberOfLimitedSamples)+" peptides in "+fInt2String(m_iNumberOfTrainProteins) + " proteins\n"); 
	
	return true;
}

bool Cstepwise::mf_loadModel( string m_strRegressionResult)
{
	cout << "\tLoading model\n";
	flog.mf_Input("\tLoading model");

	ifstream fin(m_strRegressionResult.c_str());
	if (!fin)
	{
		cout << "Error:\tCannot open " << m_strRegressionResult << endl;
		flog.mf_Input("Error:\tCannot open " + m_strRegressionResult+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0, i = 0;

	while (getline(fin, strTemp1))
	{
		strTemp2 = strTemp1.substr(0, 8);
		if (strTemp2 == "regression coefficient")
		{
			getline(fin, strTemp1);

			iTemp1 = 0;
			iTemp2 = strTemp1.find("\t", iTemp1);
			while (iTemp2 != strTemp1.npos)
			{

				strTemp2 = strTemp1.substr(iTemp1, iTemp2 - iTemp1);
				m_pCoefficient[i] = atof(strTemp2.c_str());
				i++;
				iTemp1 = iTemp2 + 1;
				iTemp2 = strTemp1.find("\t", iTemp1);
			}
			strTemp2 = strTemp1.substr(iTemp1, strTemp1.size() - iTemp1);
			m_pCoefficient[i] = atof(strTemp2.c_str());
			cout << "The Model has " << i << " coefficients.\n";
			flog.mf_Input("The Model has " +fInt2String(i)+ " coefficients.\n");

			return true;
		}
	}
	flog.mf_Input("...Done\n");
	return false;
}
void Cstepwise::mf_saveTestXY( string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "The testing set saved to " << path << endl;
	flog.mf_Input("The testing set saved to " +path +"\n");

	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;

	for (int i = 0; i<m_iTestSampleNumber; i++)
	{
		for (int j = 0; j <= m_iNumberOfFeature; j++)
			ofile << m_pTestXY[i][j] << "\t";
		ofile << "\n";
	}

	ofile.close();

}

void Cstepwise::mf_saveStepwiseX(vector<CProtein> proteins, string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "The training set saved to " << path << endl;
	flog.mf_Input("The training set saved to " + path +"\n");

	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;

	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{

		for (j = 0; j < ProteinIter->m_iPeptidesNumber; j++)
		{
			ofile << ProteinIter->m_strProteinID << "\t";
			ofile << ProteinIter->m_vPeptidesSequences.at(j) << "\t";
			for (t = 0; t <= m_iNumberOfFeature; t++)
				ofile << m_pTestXY[i][t] << "\t";
			ofile << endl;
			i++;
		}
	}
	ofile.close();

}

void Cstepwise::mf_LoadTestXLFAQpep( string LFAQpepPathByOtherMethod)
{
	cout << "\tLoading LFAQpeps!\n\n";
	flog.mf_Input("\tLoading LFAQpeps!\n\n");

	ifstream fin(LFAQpepPathByOtherMethod.c_str());
	if (!fin)
	{
		cout << "Error:\tCannot open " << LFAQpepPathByOtherMethod << endl;
		flog.mf_Input("Error:\tCannot open " + LFAQpepPathByOtherMethod +"\n");
		flog.mf_Destroy();
		exit(1);
	}
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0, i = 0;
	int iSamplesTemp = 0;

	while (getline(fin, strTemp1))
	{
		m_pTestPredictY[iSamplesTemp] = atof(strTemp1.c_str());
		m_pNativeTestPredictY[iSamplesTemp] = m_pTestPredictY[iSamplesTemp];
		iSamplesTemp++;
	}
	if (iSamplesTemp != m_iTestSampleNumber)
	{
		cout << "Error:\tThe Number of LFAQpeps does not correspond to TestX's" << endl;
		flog.mf_Input("Error:\tThe Number of LFAQpeps does not correspond to TestX's\n");
		flog.mf_Destroy();
		exit(1);
	}

}


/*
Calculate peptides quantification coefficient 
*/
void Cstepwise::mf_Prediction(vector<CProtein> &proteins)
{		
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		m_pTestPredictY[i] = 0.0;
		for (int j = 0; j < m_iNumberOfFeature; j++)
		{
			m_pTestPredictY[i] += m_pCoefficient[j] * m_pTestXY[i][j];  

		}
		m_pTestPredictY[i] += m_pCoefficient[m_iNumberOfFeature];
	
	}
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		m_pNativeTestPredictY[i] = m_pTestPredictY[i]; 
	}

	//need discretization
	Discretization(m_pTestPredictY, m_iTestSampleNumber,10);
	//Assertion(m_pTestPredictY, m_iTestSampleNumber);

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	int index = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			proteinIter->m_vPeptideQfactors.clear();
			for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
			{
				proteinIter->m_vPeptideQfactors.push_back(m_pTestPredictY[index]);
				index++;
			}
		}
			
	}

}

bool Cstepwise::mf_ConstructTestX(vector<CProtein> proteins)
{ //add this function by gzhq 20150603
	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	m_iTestSampleNumber = 0;
	m_iNumberOfTestProteins = 0;
	map<string, CPeptideAttribute> mapPeptideSequenceAndAttribute;

	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			m_iNumberOfTestProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.size(); j++)
				{
					m_pTestXY[m_iTestSampleNumber][j] = ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.at(j);
				}
				m_pTestXY[m_iTestSampleNumber][j] = ProteinIter->m_vPeptideExperimentalQFactors.at(n);  
				m_iTestSampleNumber++;
			}
		}
	}
	cout << "Constructed the testing set using " << m_iTestSampleNumber <<" peptides in "<<m_iNumberOfTestProteins<<" proteins" <<endl; 
	flog.mf_Input("Constructed the testing set using " +fInt2String(m_iTestSampleNumber)+ " peptides in " +fInt2String(m_iNumberOfTestProteins)+ " proteins\n");
	return true;
}

void Cstepwise::mf_Predict( vector<CProtein> &proteins, string RegressionResultPath)
{
	mf_loadModel(RegressionResultPath);
	mf_ConstructTestX(proteins);

	int iPosTemp = RegressionResultPath.find_last_of("\\");
	iPosTemp = RegressionResultPath.find_last_of("\\", iPosTemp - 1);
	mf_Prediction(proteins);

}

Cstepwise:: ~Cstepwise()
{
	int i;
	for (i = 0; i < m_NumberOfTrainSamples; i++) { delete[m_iNumberOfFeature + 1] m_px[i]; m_px[i] = NULL; }
	delete[m_NumberOfTrainSamples] m_px;
	m_px = NULL;
	for (i = 0; i < m_NumberOfTrainSamples; i++) { delete[m_iNumberOfFeature + 1] m_pTestXY[i]; m_pTestXY[i] = NULL; }
	for (i = 0; i < m_NumberOfTrainSamples; i++) { delete[m_iNumberOfFeature + 1] m_pTrainXForCV[i]; m_pTrainXForCV[i] = NULL; }
	delete[] m_pTrainXForCV;
	m_pTrainXForCV = NULL;

	delete[] m_pLeaveXIndex;
	delete[]m_pTrainPredictY;
	delete[]m_psortIndex;
	for (i = 0; i < m_iNumberOfFeature; i++) { delete[m_iNumberOfFeature + 1] m_pr[i]; m_pr[i] = NULL; }
	delete[m_NumberOfTrainSamples] m_pr;
	m_pr = NULL;
	delete[] m_pxx;
	delete[m_iNumberOfFeature+1]m_pCoefficient;
	delete[]m_pv;
	delete[]m_ps;
	delete[]m_pye;
	delete[]m_pyr;

	for (i = 0; i < m_iTestSampleNumber; i++) { delete[m_iNumberOfFeature + 1] m_pTestXY[i]; m_pTestXY[i] = NULL; }
	delete[m_iTestSampleNumber] m_pTestXY;
	m_pTestXY = NULL;
	delete[]m_pTestPredictY;
	delete[]m_pNativeTestPredictY; 

	cout << "*********************Stepwise end！***************************\n";

}


