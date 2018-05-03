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
#include"BARTRegression.h"
using std::cout;
using std::endl;

CBART::CBART()
{
	m_iFeatureNumber=0;
	m_iTrainSampleNumber=0;
	m_dAutoCorrelation=0.0;
	m_dMiny=0.0;
	m_dMaxy=0.0;
	m_iTestSampleNumber=0;
	m_iNumberOfTrainProteins=0;
	m_iNumberOfTestProteins=0;
	m_dRss=0.0;  
	m_dRestemp=0.0; //a residual
}
void CBART::mf_xitofile(std::ofstream& fnm, xinfo xi)
{
	fnm << xi.size() << endl;
	for (uint i = 0; i< xi.size(); i++) {
		fnm << xi[i].size() << endl;
		for (uint j = 0; j<xi[i].size(); j++) fnm << xi[i][j] << "  ";
		fnm << endl;
	}
	fnm << endl;
}
int CBART::mf_GetNumberOfTrainProteins(vector<CProtein> &proteins)
{
	vector<CProtein>::iterator ProteinIter;
	int iNumberOfTrainProteins = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			iNumberOfTrainProteins++;

		} // end if PeptidesThreshold
	} // end for proteins
	return iNumberOfTrainProteins;
}
// construct the train set 
bool CBART::mf_ConstructXY(vector<CProtein> proteins, string strExperimentName)
{
	cout << "Constructing the training set of BART\n";
	flog.mf_Input("Constructing the training set of BART\n");
	
	m_dMiny = 1e15; //using range of y to calibrate prior for bottom node mu's
	m_dMaxy = -1.0;
	vector<CProtein>::iterator  ProteinIter;
	size_t j = 0; //feature number
	int n = 0; //peptide number of per protein
	double dytemp;

	m_iTrainSampleNumber = 0;  
	m_iNumberOfTrainProteins = 0;

	m_vTrainx.clear();
	m_vTrainy.clear();
	m_cAllys.clear();
	
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_bIfInTrainSet==true)
		{
			m_iNumberOfTrainProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.size(); j++)
				{
					m_vTrainx.push_back(ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.at(j));
				}
				dytemp = ProteinIter->m_vPeptideExperimentalQFactors.at(n);
				if (dytemp < m_dMiny) m_dMiny = dytemp;
				if (dytemp > m_dMaxy) m_dMaxy = dytemp;
				m_cAllys.sy += dytemp; // sum of y
				m_cAllys.sy2 += dytemp*dytemp; // sum of y^2
				m_vTrainy.push_back(dytemp);
				m_iTrainSampleNumber++;
			}

		} 

	} // end for proteins

	if (m_iTrainSampleNumber < 1) 
	{
		cout << "Error:\tFeatureNumber<1\n";
		flog.mf_Input("Error:\tFeatureNumber<1 when constructing the train set\n");
		flog.mf_Destroy();
		exit(1);
		return false;
	}

	m_iFeatureNumber = m_vTrainx.size() / m_iTrainSampleNumber;
	if (m_vTrainx.size() != m_iTrainSampleNumber*m_iFeatureNumber)
	{
		cout << "Error:\tThere is something wrong with the dimension of the training set.\n";
		flog.mf_Input("Error:\tThere is something wrong with the dimension of the training set.\n");
		flog.mf_Destroy();
		exit(1);
	}

	m_cAllys.n = m_iTrainSampleNumber;
	cout << "\tThe number of samples in the training set is " << m_iTrainSampleNumber << endl;
	cout << "\tThe first and last of y is: ";
	cout << m_vTrainy[0] << ", " << m_vTrainy[m_iTrainSampleNumber - 1] << endl;
	string strTrainNumber = fInt2String(m_iTrainSampleNumber);
	flog.mf_Input("\tThe number of samples in the training set is " + strTrainNumber + "\n");
	flog.mf_Input("\tThe first and last of y is ");
	string strDoubleTemp = fDouble2String(m_vTrainy[0]);
	string strDoubleTemp2 = fDouble2String(m_vTrainy[m_iTrainSampleNumber - 1]);
	flog.mf_Input(strDoubleTemp+","+strDoubleTemp2+"\n");
	m_dMeany = m_cAllys.sy / m_iTrainSampleNumber;
	m_dStdy = sqrt((m_cAllys.sy2 - m_iTrainSampleNumber*m_dMeany*m_dMeany) / (m_iTrainSampleNumber - 1));
	cout << "\tThe mean of y is " << m_dMeany<<"\n";
	cout<<"\tThe standard deviation of y is " << m_dStdy << endl;
	cout << "\tConstructed the training set using " << m_iTrainSampleNumber << " samples in " << m_iNumberOfTrainProteins<<" proteins." << endl; 
	flog.mf_Input("\tThe number of samples in the training set is " + fDouble2String(m_dMeany) + "\n");
	flog.mf_Input("\tThe first and last of y is " + fDouble2String(m_dStdy) + "\n");
	flog.mf_Input("\tConstructed the training set using " + fInt2String(m_iTrainSampleNumber)+ " samples in " +fInt2String(m_iNumberOfTrainProteins )+ " proteins.\n" );

#ifdef Debug_Mode
	string BARTTrainYPath = m_strQuantiResultPath + "\\BARTTrainY_"+strExperimentName+".txt";
	SaveVector(BARTTrainYPath, m_vTrainy);
#endif
	cout << "\tThe number of features in the training set is "<< m_iFeatureNumber << endl;
	cout << "\tThe first row of training set: " << m_vTrainx[0] << " ...  " << m_vTrainx[m_iFeatureNumber - 1] << endl;
	cout << "\tThe last row of training set: " << m_vTrainx[(m_iTrainSampleNumber - 1)*m_iFeatureNumber] << " ...  " << m_vTrainx[m_iTrainSampleNumber*m_iFeatureNumber - 1] << endl;
	flog.mf_Input("\tThe number of features in the training set is " + fInt2String(m_iFeatureNumber));
	flog.mf_Input("\tThe first row of training set: " +fDouble2String(m_vTrainx[0]) + " ...  " + fDouble2String(m_vTrainx[m_iFeatureNumber - 1]) + "\n");
	flog.mf_Input("\tThe last row of training set: " + fDouble2String(m_vTrainx[(m_iTrainSampleNumber - 1)*m_iFeatureNumber])+ " ...  " +fDouble2String( m_vTrainx[m_iTrainSampleNumber*m_iFeatureNumber - 1])+"\n");
	flog.mf_Input("Done\n");
	return true;
}

bool CBART::mf_ConstructTestX(vector<CProtein> proteins, string strExperimentName)
{
	//data for predictions
	cout << "Constructing the testing set of BART\n";
	flog.mf_Input("Constructing the testing set of BART\n");

	vector<CProtein>::iterator  ProteinIter;
	size_t j = 0; //feature number
	int n = 0; //peptide number of per protein
	double dytemp;
	m_vTesty.clear();
	
	// Get the number of test peptides
	m_iTestSampleNumber = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			m_iTestSampleNumber += ProteinIter->m_iPeptidesNumber;
		} // end if PeptidesThreshold
	} // end for proteins
	try
	{
		m_pTestx = NULL;
		m_pTestx = new  double[m_iFeatureNumber*m_iTestSampleNumber];
	}
	catch (bad_alloc &memExp)
	{
		const char * pchar = memExp.what();
		string strTemp = pchar;
		cout << "Error:\tCannot allocate memory for the testing set of BART. The reason maybe out of memory!" << strTemp << "\n";
		flog.mf_Input("Error:\tCannot allocate memory for the testing set of BART. The reason maybe out of memory!" + strTemp + "\n");
		flog.mf_Destroy();
		exit(-1);
	}

	m_iNumberOfTestProteins = 0;
	int iTestxindex = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			m_iNumberOfTestProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.size(); j++)
				{
					m_pTestx[iTestxindex] = ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.at(j);
					iTestxindex++;
				}
				dytemp = ProteinIter->m_vPeptideExperimentalQFactors.at(n);
				m_vTesty.push_back(dytemp);
			}
		} // end if PeptidesThreshold
	} // end for proteins

#ifdef Debug_Mode
	string str = m_strQuantiResultPath + "\\BARTTestY.txt";
	SaveVector(str, m_vTesty);
#endif

	dip.n = 0;
	if (m_iTestSampleNumber)
	{
		dip.n = m_iTestSampleNumber;
		dip.p = m_iFeatureNumber; 
		dip.x = m_pTestx;
		dip.y = &m_vTesty[0]; 
		cout << "\tConstructed the testing set using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins" << endl;
		flog.mf_Input("\tConstructed the testing set using " +fInt2String(m_iTestSampleNumber)+ " peptides in " +fInt2String( m_iNumberOfTestProteins) + " proteins\n");
	}
	cout << "\tThe number of samples in the testing set is: "<< dip.n << endl;
	flog.mf_Input("\tThe number of samples in the testing set is: " + fSize_t2String(dip.n) + "\n");

	if (dip.n) {
		cout << "\tThe first row of the testing set: " << dip.x[0] << " ...  " << dip.x[m_iFeatureNumber - 1] << endl;
		cout << "\tThe last row of the testing set: " << dip.x[(dip.n - 1)*m_iFeatureNumber] << " ...  " << dip.x[dip.n*m_iFeatureNumber - 1] << endl;
		flog.mf_Input("\tThe first row of the testing set: " +fDouble2String(dip.x[0]) + " ...  " + fDouble2String(dip.x[m_iFeatureNumber - 1])+"\n");
		flog.mf_Input("\tThe last row of the testing set: " + fDouble2String(dip.x[(dip.n - 1)*m_iFeatureNumber]) + " ...  " + fDouble2String(dip.x[dip.n*m_iFeatureNumber - 1]) + "\n");
	}
	flog.mf_Input("Done\n");
	return true;
}

void CBART::mf_SelectFeature( vector<tree> trees,size_t burn, map<size_t, size_t> &FeatureCount)
{
	tree::cnpv npv;
	size_t inodeIndex = 0;
	map<size_t, size_t>::iterator FeatureCountIter;
	map<size_t, size_t>::iterator mapFeatureCountIter;
	vector<size_t> SortedFeatureCount;
	vector<size_t> SortedFeatureIndex;
	
	ofstream SelectedFeatureFile(m_strSelectedFeaturePath);
	if (!SelectedFeatureFile.is_open())
	{
		cout << "Error:\tCannot open SelectedFeatures.txt" << endl;
		flog.mf_Input("Error:\tCannot open SelectedFeatures.txt");
		flog.mf_Destroy();
		exit(-1);
	}
	for (size_t i = 0; i < trees.size(); i++)
	{
		npv.clear();
		trees[i].getnodes(npv);
		for (inodeIndex = 0; inodeIndex < npv.size(); inodeIndex++)
		{
			if (npv[inodeIndex]->ntype() != 'b')
			{
				FeatureCountIter = FeatureCount.find(npv[inodeIndex]->getv());
				if (FeatureCountIter != FeatureCount.end())
				{
					FeatureCountIter->second++;
				}
				else
				{
					FeatureCount.insert(pair<size_t, size_t>(npv[inodeIndex]->getv(), 1));
				}
			}
		}
	}
	
	mapFeatureCountIter = FeatureCount.begin();
	for (; mapFeatureCountIter != FeatureCount.end(); mapFeatureCountIter++)
	{
		SortedFeatureIndex.push_back(mapFeatureCountIter->first);
		SortedFeatureCount.push_back(mapFeatureCountIter->second);
	}

	quickSort(SortedFeatureCount, SortedFeatureIndex, 0, SortedFeatureIndex.size() - 1);
	for (size_t i = 0; i < SortedFeatureIndex.size(); i++)
	{
		SelectedFeatureFile << SortedFeatureIndex.at(i) << "\t" << SortedFeatureCount.at(i) << endl;
	}
	SelectedFeatureFile.close();
}
bool CBART::mf_BartRegressionRun(vector<CProtein>&proteins, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{

	cout << "Bayes addictive regression trees\n";
	cout << "\tThe parameters of BART:\n";
	cout << "\t\tburn = " << burn<<"\n";
	cout << "\t\tnd = " << nd << "\n";
	cout << "\t\tnumber of trees = " << m << "\n";
	cout << "\t\tlambda = " << lambda << "\n";
	cout << "\t\tnu = " << nu << "\n";
	cout << "\t\tkfac = " << kfac << "\n";
	flog.mf_Input("Bayes addictive regression trees\n");
	flog.mf_Input("\tThe parameters of BART:\n");
	flog.mf_Input("\t\tburn = " + fSize_t2String(burn) + "\n");
	flog.mf_Input("\t\tnd = " +fSize_t2String(nd) + "\n");
	flog.mf_Input("\t\tnumber of trees = " +fSize_t2String(m)+ "\n");
	flog.mf_Input("\t\tlambda = " + fDouble2String(lambda) + "\n");
	flog.mf_Input("\t\tnu = " + fDouble2String(nu) + "\n");
	flog.mf_Input("\t\tkfac = " +fDouble2String(kfac) + "\n");

#ifdef Debug_Mode
	string strTrainXYPath = m_strQuantiResultPath + "//TrainXY.txt";
	mf_saveTrainXY(strTrainXYPath);
	string strCheckTestX = m_strQuantiResultPath+"//TestXY.txt";
	mf_saveTestXY(strCheckTestX);
#endif

	double* pAllfit;//sum of fit of all trees
	double* pResidual; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	double* pFittemp; //fit of current tree
	//in sample fit
	double* pPredictTrainMean; //posterior mean of in-sample fit, sum draws,then divide
	//out of sample fit
	double* pPredictTesttemp; //temporary fit vector to compute prediction
	double* pPredictTestMean; //posterior mean for prediction

	//x cutpoints
	size_t nc = 100; //100 equally spaced cutpoints from min to max.
	makexinfo(m_iFeatureNumber, m_iTrainSampleNumber, &m_vTrainx[0], xi, nc);

	//trees
	tree cTreeTemp;
	cTreeTemp.setm(m_dMeany / m);
	m_vTrees.clear();
	for (size_t i = 0; i < m; i++)
	{
		m_vTrees.push_back(cTreeTemp);
	}

	//dinfo
	pAllfit = new double[m_iTrainSampleNumber]; //sum of fit of all trees
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pAllfit[i] = m_dMeany;
	pResidual = new double[m_iTrainSampleNumber]; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	pFittemp = new double[m_iTrainSampleNumber]; //fit of current tree
	di.n = m_iTrainSampleNumber; di.p = m_iFeatureNumber; di.x = &m_vTrainx[0]; di.y = pResidual; //the y for each draw will be the residual 

	pi.pb = .5; //prob of birth given  birth/death
	pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
	pi.beta = beta; //2 for bart means it is harder to build big trees.
	pi.tau = (m_dMaxy - m_dMiny) / (2 * kfac*sqrt((double)m));
	pi.sigma = m_dStdy;

	cout << "\t\talpha = " << pi.alpha << "\n";
	cout << "\t\tbeta = " << pi.beta << "\n";
	cout << "\t\tsigma = " << pi.sigma << "\n";
	cout << "\t\ttau = " << pi.tau << "\n";
	flog.mf_Input("\t\talpha = " +fDouble2String(pi.alpha)+ "\n");
	flog.mf_Input("\t\tbeta = " + fDouble2String(pi.beta) + "\n");
	flog.mf_Input("\t\tsigma = " +fDouble2String(pi.sigma) + "\n");
	flog.mf_Input("\t\ttau = " + fDouble2String(pi.tau) + "\n");

	//--------------------------------------------------
	//storage for ouput
	//in sample fit
	pPredictTrainMean = new double[m_iTrainSampleNumber]; //posterior mean of in-sample fit, sum draws,then divide
	for (size_t i = 0; i<m_iTrainSampleNumber; i++)
		pPredictTrainMean[i] = 0.0;

	//out of sample fit
	pPredictTestMean = new double[dip.n];
	pPredictTesttemp = new double[dip.n];
	for (size_t i = 0; i<dip.n; i++) 
		pPredictTestMean[i] = 0.0;

	//for sigma draw
	double rss;  //residual sum of squares
	double restemp; //a residual

	//--------------------------------------------------
	//mcmc
	//random number generation
	uint seed = 99;
	RNG gen(seed); //this one random number generator is used in all draws

	cout << "\tMCMC:\n";
	flog.mf_Input("\tMCMC:\n");
	clock_t tp;
	tp = clock();
	for (size_t i = 0; i<(nd + burn); i++) 
	{
		if (i % 100 == 0)
		{
			cout << "\t\ti: " << i << endl;
			flog.mf_Input("\t\ti: " +fInt2String(i)+"\n");
		}
			
		//draw trees
		for (size_t j = 0; j<m; j++) 
		{
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
			{
				pAllfit[k] = pAllfit[k] - pFittemp[k];
				pResidual[k] = m_vTrainy[k] - pAllfit[k];
			}
			bd(m_vTrees[j], xi, di, pi, gen);
			drmu(m_vTrees[j], xi, di, pi, gen);
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
				pAllfit[k] += pFittemp[k];
		}
		//draw sigma
		rss = 0.0;
		for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
		{ 
			restemp = m_vTrainy[k] - pAllfit[k];
			rss += restemp*restemp;
		}
		pi.sigma = sqrt((nu*lambda + rss) / gen.chi_square(nu + m_iTrainSampleNumber));
		if (i >= burn) 
		{
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
				pPredictTrainMean[k] += pAllfit[k];
			if (dip.n) 
			{
				for (size_t j = 0; j<m; j++) 
				{
					fit(m_vTrees[j], xi, dip, pPredictTesttemp);
					for (size_t k = 0; k<dip.n; k++)
						pPredictTestMean[k] += pPredictTesttemp[k];
				}
			}
		}
	}
	tp = clock() - tp;
	double thetime = (double)(tp) / (double)(CLOCKS_PER_SEC);
	cout << "Solved the BART model using MCMC took " << thetime <<" seconds"<< endl;
	flog.mf_Input("Solved the BART model using MCMC took " + fDouble2String(thetime) + " seconds\n");

	m_vTrainPredicty.clear();
	for (size_t k = 0; k < m_iTrainSampleNumber; k++)
	{
		pPredictTrainMean[k] /= (double)nd;
		m_vTrainPredicty.push_back(pPredictTrainMean[k]);
	}
	m_dAutoCorrelation = spearsonCorrelation(m_vTrainPredicty, m_vTrainy);

#ifdef Debug_Mode
	string strBARTTrainPredictedy = m_strQuantiResultPath + "\\BARTTrainPredictedy.txt";
	SaveVector(strBARTTrainPredictedy, m_vTrainPredicty);
#endif

	if (dip.n) {
		for (size_t k = 0; k < dip.n; k++)
		{
			pPredictTestMean[k] /= (double)nd;
			m_vTestPredicty.push_back(pPredictTestMean[k]);
		}
	}

	map<size_t, size_t> mapFeatureCount;
	map<size_t, size_t>::iterator mapFeatureCountIter;
	vector<size_t> SortedFeatureCount;
	vector<size_t> SortedFeatureIndex;

	mf_SelectFeature(m_vTrees,burn, mapFeatureCount);

	mapFeatureCountIter = mapFeatureCount.begin();
	for (; mapFeatureCountIter != mapFeatureCount.end(); mapFeatureCountIter++)
	{
		SortedFeatureIndex.push_back(mapFeatureCountIter->first);
		SortedFeatureCount.push_back(mapFeatureCountIter->second);
	}

	quickSort(SortedFeatureCount, SortedFeatureIndex, 0, SortedFeatureIndex.size() - 1);

	Discretization(pPredictTestMean, m_iTestSampleNumber,10);
	//Assertion(pPredictTestMean, m_iTestSampleNumber);
	//Scaling(pPredictTestMean, m_iTestSampleNumber);

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
				proteinIter->m_vPeptideQfactors.push_back(pPredictTestMean[index]);
				index++;
			}
		}
		
	}

	delete []pAllfit;
	pAllfit = NULL;
	delete []pResidual;
	pResidual = NULL;
	delete []pFittemp;
	pFittemp = NULL;
	delete []pPredictTrainMean;
	pPredictTrainMean = NULL;
	delete []pPredictTestMean;
	pPredictTestMean = NULL;
	delete []pPredictTesttemp;
	pPredictTesttemp = NULL;
	delete[]m_pTestx;
	m_pTestx = NULL;
	return true;
}


void CBART::mf_MakeIndexPartition(size_t SampleNumber, vector<size_t> &IndexPartation)
{
	for (size_t i = 0; i < SampleNumber; i++)
	{
		IndexPartation.push_back(floor(rand()%10));
	}
}

void CBART::mf_ConstructCVData( vector<CProtein> proteins, vector<size_t>IndexPartation, size_t fold, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double>&TestyForCV)
{
	m_dMiny = INFINITY; //using range of y to calibrate prior for bottom node mu's
	m_dMaxy = -INFINITY;

	vector<CProtein>::iterator  ProteinIter;
	size_t j = 0; //feature num
	int n = 0; //peptide num per protein
	double dytemp;

	TrainxForCV.clear();
	TrainYForCV.clear();
	TestxForCV.clear();
	m_iTrainSampleNumber = 0;  
	m_iTestSampleNumber = 0;
	m_cAllys.clear();
	m_iNumberOfTestProteins = 0;
	m_iNumberOfTrainProteins = 0;
	int i = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{		
		if (ProteinIter->m_bIfInTrainSet==true)
		{
			if (IndexPartation.at(i) == fold)
			{
				for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
				{
					for (j = 0; j < ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.size(); j++)
					{
						TestxForCV.push_back(ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.at(j));
					}
					dytemp = ProteinIter->m_vPeptideExperimentalQFactors.at(n);
					TestyForCV.push_back(dytemp);
					m_iTestSampleNumber++;
				}
				m_iNumberOfTestProteins++;
			}
			else
			{
				for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
				{
					for (j = 0; j < ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.size(); j++)
					{
						TrainxForCV.push_back(ProteinIter->m_vPeptidesFeatures.at(n).m_vecFeatures.at(j));
					}
					dytemp = ProteinIter->m_vPeptideExperimentalQFactors.at(n);
					if (dytemp < m_dMiny) m_dMiny = dytemp;
					if (dytemp > m_dMaxy) m_dMaxy = dytemp;
					m_cAllys.sy += dytemp; // sum of y
					m_cAllys.sy2 += dytemp*dytemp; // sum of y^2
					TrainYForCV.push_back(dytemp);
					m_iTrainSampleNumber++;
				}
				m_iNumberOfTrainProteins++;
			}
			i++;

		} // end if PeptidesThreshold

	} // end for proteins

	m_iTrainSampleNumber = TrainYForCV.size();

	m_cAllys.n = m_iTrainSampleNumber;
	m_dMeany = m_cAllys.sy / m_iTrainSampleNumber;
	m_dStdy = sqrt((m_cAllys.sy2 - m_iTrainSampleNumber*m_dMeany*m_dMeany) / (m_iTrainSampleNumber - 1));

	m_iFeatureNumber = TrainxForCV.size() / m_iTrainSampleNumber;

	if (TrainxForCV.size() != m_iTrainSampleNumber*m_iFeatureNumber) 
	{
		cout << "Error: input x file has wrong number of values\n";
		flog.mf_Input("Error:\t input x file has wrong number of values\n");
	}

	if (m_iTestSampleNumber)
	{
		dip.n = m_iTestSampleNumber;
		dip.p = m_iFeatureNumber; 
		dip.x = &TestxForCV[0];
		dip.y = &TestyForCV[0]; 
		cout << "Constructed the testing set using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins" << endl;
		flog.mf_Input("Constructed the testing set using "+ fInt2String(m_iTestSampleNumber) + " peptides in " + fInt2String(m_iNumberOfTestProteins) + " proteins\n");
	}

}
void CBART::mf_TrainAndPredictForCV( vector<CProtein> &proteins, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double> &PredictedY, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{
	double* pAllfit;//sum of fit of all trees
	double* pResidual; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	double* pFittemp; //fit of current tree
	//in sample fit
	double* pPredictTrainMean; //posterior mean of in-sample fit, sum draws,then divide
	//out of sample fit
	double* pPredictTesttemp; //temporary fit vector to compute prediction
	double* pPredictTestMean; //posterior mean for prediction

	//x cutpoints
	size_t nc = 100; //100 equally spaced cutpoints from min to max.
	makexinfo(m_iFeatureNumber, m_iTrainSampleNumber, &TrainxForCV[0], xi, nc);

	//trees
	tree cTreeTemp;
	cTreeTemp.setm(m_dMeany / m);
	m_vTrees.clear();
	for (size_t i = 0; i < m; i++)
	{
		m_vTrees.push_back(cTreeTemp);
	}

	//dinfo
	pAllfit = new double[m_iTrainSampleNumber]; //sum of fit of all trees
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pAllfit[i] = m_dMeany;
	pResidual = new double[m_iTrainSampleNumber]; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	pFittemp = new double[m_iTrainSampleNumber]; //fit of current tree
	di.n = m_iTrainSampleNumber; di.p = m_iFeatureNumber; di.x = &TrainxForCV[0]; di.y = pResidual; //the y for each draw will be the residual 

	pi.pb = .5; //prob of birth given  birth/death
	pi.alpha =alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
	pi.beta = beta; //2 for bart means it is harder to build big trees.

	pi.tau = (m_dMaxy - m_dMiny) / (2 * kfac*sqrt((double)m));
	pi.sigma = m_dStdy;	
	cout << "\talpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	cout << "\tsigma, tau: " << pi.sigma << ", " << pi.tau << endl;
	flog.mf_Input("\talpha, beta: " + fDouble2String(pi.alpha) + ", " + fDouble2String(pi.beta) + "\n");
	flog.mf_Input("\tsigma, tau: " + fDouble2String(pi.sigma) + ", " + fDouble2String(pi.tau) + "\n");

	//--------------------------------------------------
	//storage for ouput in sample fit
	pPredictTrainMean = new double[m_iTrainSampleNumber]; //posterior mean of in-sample fit, sum draws,then divide
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pPredictTrainMean[i] = 0.0;

	//out of sample fit
	pPredictTestMean = new double[dip.n];
	pPredictTesttemp = new double[dip.n];
	for (size_t i = 0; i<dip.n; i++) pPredictTestMean[i] = 0.0;

	//for sigma draw
	double rss;  //residual sum of squares
	double restemp; //a residual

	//--------------------------------------------------
	//mcmc
	//random number generation
	uint seed = 99;
	RNG gen(seed); //this one random number generator is used in all draws

	cout << "\tMCMC:\n";
	flog.mf_Input("\tMCMC:\n");
	clock_t tp;
	tp = clock();
	for (size_t i = 0; i<(nd + burn); i++)
	{
		if (i % 100 == 0)
		{
			cout << "\t\ti: " << i << endl;
			flog.mf_Input("\t\ti: " + fInt2String(i) + "\n");
		}
		//draw trees
		for (size_t j = 0; j<m; j++)
		{
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
			{
				pAllfit[k] = pAllfit[k] - pFittemp[k];
				pResidual[k] = TrainYForCV[k] - pAllfit[k];
			}
			bd(m_vTrees[j], xi, di, pi, gen);
			drmu(m_vTrees[j], xi, di, pi, gen);
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) pAllfit[k] += pFittemp[k];
		}
		//draw sigma
		rss = 0.0;
		for (size_t k = 0; k<m_iTrainSampleNumber; k++)
		{
			restemp = TrainYForCV[k] - pAllfit[k];
			rss += restemp*restemp;
		}
		pi.sigma = sqrt((nu*lambda + rss) / gen.chi_square(nu + m_iTrainSampleNumber));
		if (i >= burn)
		{
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
				pPredictTrainMean[k] += pAllfit[k];
			if (dip.n)
			{
				for (size_t j = 0; j<m; j++)
				{
					fit(m_vTrees[j], xi, dip, pPredictTesttemp);
					for (size_t k = 0; k<dip.n; k++) pPredictTestMean[k] += pPredictTesttemp[k];
				}
			}
		}
	}
	tp = clock() - tp;
	double thetime = (double)(tp) / (double)(CLOCKS_PER_SEC);

 	for (size_t k = 0; k<m_iTrainSampleNumber; k++) pPredictTrainMean[k] /= (double)nd;
	if (dip.n) 
	{
		for (size_t k = 0; k < dip.n; k++)
		{
			pPredictTestMean[k] /= (double)nd;
			PredictedY.push_back(pPredictTestMean[k]);
		}
	}
	
	delete []pAllfit;
	pAllfit = NULL;
	delete pResidual;
	pResidual = NULL;
	delete []pFittemp;
	pFittemp = NULL;
	delete []pPredictTrainMean;
	pPredictTrainMean = NULL;
	delete []pPredictTestMean;
	pPredictTestMean = NULL;
	delete []pPredictTesttemp;
	pPredictTesttemp = NULL;
	
}

void CBART::mf_CorrectedPeptidesIntensity(vector<CProtein>& proteins)
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dPeptideIntensityTemp;
	int i = 0;
	proteinIter = proteins.begin();
	for (; proteinIter != proteins.end(); proteinIter++)
	{
		i = 0;
		proteinIter->m_dProteinLFAQpep = 0.0;
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			proteinIter->m_bIfCalculateLFAQpep = true;
			if (proteinIter->m_iPeptidesNumber > UniquePeptidesCorrectThreshold)
			{			
				for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++)
				{
					dPeptideIntensityTemp = *peptidesIter;
					if (proteinIter->m_vPeptideQfactors.at(i) != 0.0)
						proteinIter->m_vPeptidesCorrectedIntensity[i] = dPeptideIntensityTemp / proteinIter->m_vPeptideQfactors.at(i);
					else
						proteinIter->m_vPeptidesCorrectedIntensity[i] = 0.0;
					proteinIter->m_dProteinLFAQpep += proteinIter->m_vPeptidesCorrectedIntensity[i];
					i++;
				}
			}
			else
			{ //not enough peptides for correction
				for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++)
				{
					dPeptideIntensityTemp = *peptidesIter;
					proteinIter->m_dProteinLFAQpep += dPeptideIntensityTemp;
				}
			}
			if (proteinIter->m_iPeptidesNumber > 0)
			{
				proteinIter->m_dProteinLFAQpep = proteinIter->m_dProteinLFAQpep / proteinIter->m_iPeptidesNumber;
			}
			else
			{
				proteinIter->m_dProteinLFAQpep = 0.0;
				cout << "\tWarning:\tProtein " << proteinIter->m_strProteinID << " do not have unique peptides.\n";
				cout << "\tSo we make its LFAQ 0.0.\n";
				flog.mf_Input("\tWarning:\tProtein " + proteinIter->m_strProteinID + " do not have unique peptides. \n");
				flog.mf_Input("So we make its LFAQ 0.0.\n");
			}
		}
	}
}

 //mf_ChooseBestParameter: select the parameters : pi.alpha pi.beta
//using ten-fold crossvalidation for parameters selection
//void CBART::mf_ChooseBestParameter(vector<CProtein> &proteins, CQuantificationParam trainparam)
//{
//	cout << "\tChoosing the optimal parameters\n";
//	flog.mf_Input("\tChoosing the optimal parameters\n");
//
//	size_t burn = 100;
//	size_t nd = 100;
//	size_t m = 200;
//	double lambda = 1;
//	double nu = 3;
//	double kfac = 2;
//	double alpha = 0.95;
//	double beta = 1;
//	vector<size_t> vecIndexPartation;
//	
//	string strParameterSelectByCVResultPath = trainparam.m_strQuantiResultPath + "\\ParameterOptimization.txt";
//	ofstream ParameterSelectByCVResult(strParameterSelectByCVResultPath);
//	ParameterSelectByCVResult << "lambda\tkfac\tnu\tm\t";
//	ParameterSelectByCVResult << "alpha\tbeta\tSquareSumOfRelativeError\n";
//
//	vector<double> vecTrainxForCV;
//	vector<double> vecTrainYForCV;
//	vector<double> vecTestXForCV;
//	vector<double> vecTestYForCV;
//	vector<double> vecPredictedY;
//	int iFold = 0;
//	double dBestSpearsonCorrelation = 0.0;
//	double dBestUPS2SpearsonCorr=0.0;
//	double dSpearsonCorrelationTemp;
//	double dUPS2PearsonCorr;
//	double dbestLamda;
//	double dBestNu;
//	double dBestalpha;
//	double dBestBeta;
//	size_t siBestm;
//	string strTemp;
//	double dMinREsq=1e15;
//	double dREsqTemp;
//
//	map<string, double> mapUPS2Id2Mols;
//	map<string, double> mapUPS2iBAQs;
//	map<string, double>::iterator Ups2MolsIter;
//	map<string, double>::iterator Ups2WMIter;
//	map<string, double>::iterator UPS2iBAQIter;
//
//	double meanReal = 0.0;
//	double meanTest = 0.0;
//	double sd1 = 0.0, sd2 = 0.0;
//	vector<string> vecUPS2IDIdentified;
//	vector<double> vecUPS2logPredictIdentified;
//	vector<double> vecUPS2logMolsIdentified;
//
//	string::size_type stPosition;
//	vector<CProtein>::iterator proteinIter;
//	vector<double>::iterator peptidesIter;
//	vector<double>::iterator PeptideQfactorsIter;
//	vector<double> CorrectedPeptidesIntensity;
//	vector<double> NativePeptidesIntensity;
//	double dPeptideIntensityTemp;
//	int index = 0;
//	int i = 0;
//	map<string, int>::iterator mapLowSpikedInUPS2Iter;
//	CProteinWorker proteinworker;
//	proteinworker.m_Params = trainparam;
//	if (trainparam.m_bIfCotainStandardProtein)
//	{
//		proteinworker.mf_LoadUPS2Mols(trainparam.m_strStandardProteinsFilePath, mapUPS2Id2Mols);
//	}
//	m_iNumberOfTrainProteins = mf_GetNumberOfTrainProteins(proteins);
//	mf_MakeIndexPartition(m_iNumberOfTrainProteins, vecIndexPartation);
//	//for (burn = 100; burn <= 500; burn = burn + 100)
//	//{
//	//	for (m = 100; m < 1000;m=m+100)
//		for (alpha = 0.8; alpha < 1; alpha = alpha + 0.05)
//			for (beta = 0.1; beta <= 2; beta = beta + 0.5)
//			{
//				if (!trainparam.m_bIfCotainStandardProtein)
//				{
//					iFold = 0;
//					vecTestYForCV.clear();
//					vecPredictedY.clear();
//					//for 10-fold CrossValidation
//					while (iFold <= 9)
//					{
//						mf_ConstructCVData(proteins, vecIndexPartation, iFold, vecTrainxForCV, vecTrainYForCV, vecTestXForCV, vecTestYForCV);
//						mf_TrainAndPredictForCV(proteins, vecTrainxForCV, vecTrainYForCV, vecTestXForCV, vecPredictedY, burn, nd, m, lambda, nu, kfac, alpha, beta);
//						iFold++;
//					}
//					dSpearsonCorrelationTemp = spearsonCorrelation(vecPredictedY, vecTestYForCV);
//
//					if (dSpearsonCorrelationTemp > dBestSpearsonCorrelation)
//					{
//						dBestSpearsonCorrelation = dSpearsonCorrelationTemp;
//						dbestLamda = lambda;
//						dBestNu = nu;
//						dBestalpha = alpha;
//						dBestBeta = beta;
//						siBestm = m;
//					}
//					ParameterSelectByCVResult << lambda << "\t" << kfac << "\t";
//					ParameterSelectByCVResult << nu << "\t" << m << "\t";
//					ParameterSelectByCVResult << alpha << "\t" << beta << "\t";
//					ParameterSelectByCVResult << dSpearsonCorrelationTemp << "\t";
//					ParameterSelectByCVResult << endl;
//				}
//				else
//				{
//					mf_ConstructXY(proteins);
//					mf_ConstructTestX(proteins);
//					mf_BartRegressionRun(proteins, burn, nd, m, lambda, nu, kfac, alpha, beta);
//					mf_CorrectedPeptidesIntensity(proteins);
//					proteinworker.mf_PredictMolsByStandAndLFAQ(proteins, dREsqTemp);
//					if (dREsqTemp < dMinREsq)
//					{
//							dMinREsq = dREsqTemp;
//							dbestLamda = lambda;
//							dBestNu = nu;
//							dBestalpha = alpha;
//							dBestBeta = beta;
//							siBestm = m;
//					}
//				}			
//
//				ParameterSelectByCVResult << lambda << "\t" << kfac << "\t";
//				ParameterSelectByCVResult << nu << "\t" << m << "\t";
//				ParameterSelectByCVResult << alpha << "\t" << beta << "\t";
//				ParameterSelectByCVResult << dREsqTemp << "\t";
//				ParameterSelectByCVResult << endl;
//			} //end for parameters 
//
//	cout << "\t\tThe best parameters of Bayes addictive regression trees is:"<<siBestm<<"\n";
//	cout << "\t\tdbestLamda = " << dbestLamda << endl;
//	cout << "\t\tdBestNu = " << dBestNu << endl;
//	cout << "\t\tdBestalpha = " << dBestalpha << endl;
//	cout << "\t\tdBestBeta = " << dBestBeta << endl;
//
//	flog.mf_Input("\t\tThe best parameters of Bayes addictive regression trees is:"+fDouble2String(siBestm)+"\n");
//	flog.mf_Input("\t\tdbestLamda = " + fDouble2String(dbestLamda)+"\n");
//	flog.mf_Input("\t\tdBestNu = " +fDouble2String(dBestNu)+"\n");
//	flog.mf_Input("\t\tdBestalpha = " +fDouble2String(dBestalpha)+"\n");
//	flog.mf_Input("\t\tdBestBeta = " +fDouble2String(dBestBeta)+"\n");
//
//	m_dBestAlpha = dBestalpha;
//	m_dBestBeta = dBestBeta;
//	m_iBestNumberOfTrees = siBestm;
//}

void CBART::mf_saveTestXY( string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "The testing set saved to " << path << endl;
	flog.mf_Input("The testing set saved to " +path+"\n");
	vector<CProtein>::iterator ProteinIter;


	int t = 0;

	for (size_t i = 0; i<m_iTestSampleNumber; i++)
	{
		for (size_t j = 0; j < m_iFeatureNumber; j++)
		{
			ofile << m_pTestx[j + i*m_iFeatureNumber] << "\t";
		}
		ofile << m_vTesty[i];
		ofile << "\n";
	}

	ofile.close();

}

void CBART::mf_saveRegressionResult( string path)
{
	ofstream fout(path.c_str(), ios::out);
	if (!fout)
	{
		cout << "Error:\tCannot open  " << path << endl;
		flog.mf_Input("Error:\tCannot open  " + path+"\n");
		flog.mf_Destroy();
		exit(1);
	}
	cout << "The BART regression result saved to " << path << endl;
	flog.mf_Input("The BART regression result saved to " + path+"\n");

	double RMSE = 0.0;
	fout << "Trainy\tPredictedy\n";
	for (size_t i = 0; i < m_iTrainSampleNumber; i++)
	{
		fout << m_vTrainy[i] << "\t" << m_vTrainPredicty[i] << endl;
		RMSE += (m_vTrainy[i] - m_vTrainPredicty[i])*(m_vTrainy[i] - m_vTrainPredicty[i]);
	}
	RMSE = sqrt(RMSE / m_iTrainSampleNumber);
	fout << "The RMSE of BART regression=" << RMSE << endl;

	double Correlation = spearsonCorrelation(m_vTrainy, m_vTrainPredicty);
	fout << "The Correlation between Prediction Y and Real Y(AutoTrain) = " << Correlation << endl;
	cout << "\tThe Correlation between Prediction Y and Real Y(AutoTrain) = " << Correlation << endl;
	flog.mf_Input("\tThe Correlation between Prediction Y and Real Y(AutoTrain) = " +fDouble2String(Correlation)+"\n");

	fout << "Testy\tPredictedy\n";
	for (size_t i = 0; i < m_iTestSampleNumber; i++)
	{
		fout << m_vTesty[i] << "\t" << m_vTestPredicty[i] << endl;
		RMSE += (m_vTesty[i] - m_vTestPredicty[i])*(m_vTesty[i] - m_vTestPredicty[i]);
	}
	RMSE = sqrt(RMSE / m_iTestSampleNumber);
	fout << "The RMSE of BART regression in test dataset=" << RMSE << endl;

	Correlation = spearsonCorrelation(m_vTesty, m_vTestPredicty);
	fout << "The Correlation between Prediction Y and Real Y(in test dataset) = " << Correlation << endl;
	cout << "\tThe Correlation between Prediction Y and Real Y(in test dataset) = " << Correlation << endl;
	flog.mf_Input("\tThe Correlation between Prediction Y and Real Y(in test dataset) = " +fDouble2String(Correlation)+"\n");
	fout.close();
}


void CBART::mf_saveTrainXY( string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "\tThe training set saved to " << path << endl;
	flog.mf_Input("\tThe training set saved to " + path+"\n");
	vector<CProtein>::iterator ProteinIter;

	int t = 0;

	for (size_t i = 0; i<m_iTrainSampleNumber; i++)
	{
		for (size_t j = 0; j < m_iFeatureNumber; j++)
		{
			ofile << m_vTrainx[j + i*m_iFeatureNumber] << "\t";
		}
		ofile << m_vTrainy[i];
		ofile << "\n";
	}
	ofile.close();
}

// Iteration for better result
bool CBART::mf_BARTIteration(vector<CProtein> &proteins, CQuantificationParam trainparam, size_t burn, size_t nd, \
	size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{
	// observe the changed trend of AutoCorrelation and peptides'CV with the iteration;
	ofstream pPeptidesCVFile("PeptidesCV.txt");
	if (!pPeptidesCVFile.is_open())
	{
		cout << "Error:\tCannot open peptidesCV.txt" << endl;
		flog.mf_Input("Error:\tCannot open peptidesCV.txt\n");
		flog.mf_Destroy();
		exit(1);
	}
	pPeptidesCVFile << "Iterations\tFirstQuantile\tMedian\tThirdQuantile\tAutoCorrelation\tPearsonCorrelationUPS2PredictedWithMols\n";
	double dFirstQuan=0.0;
	double dMedian=0.0;
	double dThirdQuan=0.0;

	mf_ConstructXY(proteins,"iteration");
	mf_BARTAutoRegressionRun(proteins, burn, nd, m, lambda, nu, kfac, alpha, beta);
	mf_CorrectedPeptidesIntensity(proteins);
	mf_CalculatePeptidesCVofProteins(proteins, dFirstQuan, dMedian, dThirdQuan);
	
	pPeptidesCVFile << 0 << "\t" << dFirstQuan << "\t" << dMedian << "\t" << dThirdQuan << "\t";
	pPeptidesCVFile << m_dAutoCorrelation << "\t" << m_dPearsonCorrelationUPS2PredictedWithMols << "\n";

	for (int i = 1; i < 20; i++)
	{
		mf_ConstructXY(proteins,"iteration");
		mf_BARTAutoRegressionRun(proteins, burn, nd, m, lambda, nu, kfac, alpha, beta);
		mf_CorrectedPeptidesIntensity(proteins);
		mf_CalculatePeptidesCVofProteins(proteins, dFirstQuan, dMedian, dThirdQuan);		
		pPeptidesCVFile << i << "\t" << dFirstQuan<<"\t"<<dMedian<<"\t"<<dThirdQuan <<"\t";
		pPeptidesCVFile << m_dAutoCorrelation << "\t" << m_dPearsonCorrelationUPS2PredictedWithMols << "\n";
	}
	return true;
}

bool CBART::mf_BARTAutoRegressionRun(vector<CProtein> &proteins, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{

	cout << "Bayes addictive regression trees\n";
	cout << "burn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
	cout << "lambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;
	flog.mf_Input("Bayes addictive regression trees\n");
	flog.mf_Input("burn,nd,number of trees: " + fSize_t2String(burn) +", " +fSize_t2String(nd) + ", " +fSize_t2String(m)+"\n");
	flog.mf_Input("lambda,nu,kfac: " +fDouble2String(lambda)+ ", " +fDouble2String(nu)+ ", " +fDouble2String(kfac)+"\n");

	double* pAllfit;//sum of fit of all trees
	double* pResidual; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	double* pFittemp; //fit of current tree
	double* pPredictTrainMean; //posterior mean of in-sample fit, sum draws,then divide

	//x cutpoints
	size_t nc = 100; //100 equally spaced cutpoints from min to max.
	makexinfo(m_iFeatureNumber, m_iTrainSampleNumber, &m_vTrainx[0], xi, nc);

	//trees
	tree cTreeTemp;
	cTreeTemp.setm(m_dMeany / m);
	m_vTrees.clear();
	for (size_t i = 0; i < m; i++)
	{
		m_vTrees.push_back(cTreeTemp);
	}

	pi.pb = .5; //prob of birth given  birth/death

	pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
	pi.beta = beta; //2 for bart means it is harder to build big trees.
	pi.tau = (m_dMaxy - m_dMiny) / (2 * kfac*sqrt((double)m));
	pi.sigma = m_dStdy;
	cout << "\talpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	cout << "\tsigma, tau: " << pi.sigma << ", " << pi.tau << endl;
	flog.mf_Input("\talpha, beta: " + fDouble2String(pi.alpha) + ", " + fDouble2String(pi.beta)+"\n");
	flog.mf_Input("\tsigma, tau: " +fDouble2String(pi.sigma)+ ", " +fDouble2String(pi.tau)+"\n");

	//--------------------------------------------------
	//dinfo
	pAllfit = new double[m_iTrainSampleNumber]; //sum of fit of all trees
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pAllfit[i] = m_dMeany;
	pResidual = new double[m_iTrainSampleNumber]; //y-(m_pAllfit-ftemp) = y-m_pAllfit+ftemp
	pFittemp = new double[m_iTrainSampleNumber]; //fit of current tree
	di.n = m_iTrainSampleNumber; di.p = m_iFeatureNumber; di.x = &m_vTrainx[0]; di.y = pResidual; //the y for each draw will be the residual 

	//--------------------------------------------------
	//storage for ouput in sample fit
	pPredictTrainMean = new double[m_iTrainSampleNumber]; //posterior mean of in-sample fit, sum draws,then divide
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pPredictTrainMean[i] = 0.0;

	//for sigma draw
	double rss;  //residual sum of squares
	double restemp; //a residual

	//--------------------------------------------------
	//mcmc
	//random number generation
	uint seed = 99;
	RNG gen(seed); //this one random number generator is used in all draws

	cout << "\tMCMC:\n";
	flog.mf_Input("\tMCMC:\n");

	clock_t tp;
	tp = clock();
	for (size_t i = 0; i<(nd + burn); i++)
	{
		if (i % 100 == 0)
		{
			cout << "i: " << i << endl;
			flog.mf_Input("i: " +fInt2String(i)+"\n");

		}
		//draw trees
		for (size_t j = 0; j<m; j++)
		{
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
			{
				pAllfit[k] = pAllfit[k] - pFittemp[k];
				pResidual[k] = m_vTrainy[k] - pAllfit[k];
			}
			bd(m_vTrees[j], xi, di, pi, gen);
			drmu(m_vTrees[j], xi, di, pi, gen);
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) pAllfit[k] += pFittemp[k];
		}
		//draw sigma
		rss = 0.0;
		for (size_t k = 0; k<m_iTrainSampleNumber; k++)
		{
			restemp = m_vTrainy[k] - pAllfit[k];
			rss += restemp*restemp;
		}
		pi.sigma = sqrt((nu*lambda + rss) / gen.chi_square(nu + m_iTrainSampleNumber));
		if (i >= burn)
		{
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
				pPredictTrainMean[k] += pAllfit[k];
		}
	}
	tp = clock() - tp;
	double thetime = (double)(tp) / (double)(CLOCKS_PER_SEC);
	cout << "time for loop: " << thetime << endl;
	flog.mf_Input("time for loop: " +fDouble2String(thetime)+"\n");

	std::ofstream timef("time.txt");
	timef << thetime << endl;

	m_vTrainPredicty.clear();
	for (size_t k = 0; k < m_iTrainSampleNumber; k++)
	{
		pPredictTrainMean[k] /= (double)nd;
		m_vTrainPredicty.push_back(pPredictTrainMean[k]);
	}

#ifdef Debug_Mode
	std::ofstream treef("trees.txt");
	//treef << xi << endl; //the cutpoints
	mf_xitofile(treef, xi);
	treef << m << endl;  //number of trees
	treef << m_iFeatureNumber << endl;  //dimension of x's
	for (size_t j = 0; j<m; j++) treef << m_vTrees[j] << endl;  //all the trees
#endif

	m_dAutoCorrelation = spearsonCorrelation(m_vTrainy, pPredictTrainMean);
	cout << "\tThe Correlation between Prediction Y and Real Y is " << m_dAutoCorrelation << endl;
	flog.mf_Input("\tThe Correlation between Prediction Y and Real Y is " +fDouble2String(m_dAutoCorrelation)+"\n");

	Discretization(pPredictTrainMean, m_iTrainSampleNumber,10);
	//Assertion(pPredictTestMean, m_iTestSampleNumber);

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	int index = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfInTrainSet==true)
		{
			proteinIter->m_vPeptideQfactors.clear();
			for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) 
			{
				proteinIter->m_vPeptideQfactors.push_back(pPredictTrainMean[index]);
				index++;
			}
		}

	}

	delete pAllfit;
	pAllfit = NULL;
	delete pResidual;
	pResidual = NULL;
	delete pFittemp;
	pFittemp = NULL;
	delete pPredictTrainMean;
	pPredictTrainMean = NULL;

	return true;
}
//calculate the peptides CV of every protein, then get the meadian and First and third quartile of all CVs;
bool CBART::mf_CalculatePeptidesCVofProteins(vector<CProtein> proteins, double &FirstQ, double &median, double &ThirdQ)
{

	vector<double> vecPeptidesCVs;
	vector<double> vecPeptidesIntensityOfOneProtein;
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	vector<double>::iterator PeptideQfactorsIter;

	double dPeptideIntensityTemp;
	double dCvTemp;

	int i = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		i = 0;
		if (proteinIter->m_bIfInTrainSet==true)
		{
			for (peptidesIter = proteinIter->m_vPeptidesCorrectedIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesCorrectedIntensity.end(); peptidesIter++) 
			{
				dPeptideIntensityTemp = *peptidesIter;
				vecPeptidesIntensityOfOneProtein.push_back(dPeptideIntensityTemp);
			}
			dCvTemp = CalculatelogCV(vecPeptidesIntensityOfOneProtein);
			vecPeptidesIntensityOfOneProtein.clear();
			vecPeptidesCVs.push_back(dCvTemp);
		}
		
	}

	GetQuantiles(vecPeptidesCVs, FirstQ, median, ThirdQ);
	return true;

}
