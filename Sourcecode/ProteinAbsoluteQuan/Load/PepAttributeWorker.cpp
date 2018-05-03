#include "StdAfx.h"
#include "PepAttributeWorker.h"
#include "Constants.h"
#include<hash_map>
#include<algorithm>
#include<cmath>
//using namespace std;
PepAttributeWorker::PepAttributeWorker()
{
	m_strAAindexFilePath = "aaindex1.formatted";
}


PepAttributeWorker::~PepAttributeWorker(void)
{
}
vector<string> PepAttributeWorker::mf_LoadAAindex()
{
	vector<string > vStrTemp;
	//cout << "Loading AAIndex!\n";

	FILE * pFile;
	pFile = fopen(m_strAAindexFilePath.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Error when opening \"" << m_strAAindexFilePath << endl;
	}
	//cout<<"open "<<strProteinFilePath<<endl;
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	fgets(Buffer, BUFFERLENGTH, pFile);//首行
	fgets(Buffer, BUFFERLENGTH, pFile);
	int row = 0, colume = 0;
	while (!feof(pFile))
	{
		//cout << "row " << row << endl;
		pstr = Buffer;
		colume = 0;
		//cout << pstr << endl;
		pstr1 = strstr(pstr, "\t");
		*pstr1 = '\0';
		vStrTemp.push_back(pstr);
		pstr = pstr1 + 1;

		//cout << pstr << endl;

		pstr1 = strstr(pstr, "\t");
		pstr = pstr1 + 1;
		//cout << pstr << endl;

		pstr1 = strstr(pstr, "\t");
		while (pstr1 != NULL)
		{
			*pstr1 = '\0';
			//cout << pstr << endl;
			if (pstr != "NA")
			{
				m_arrayAAindexAttributeValue[row][colume] = atof(pstr);
				m_arrayIfIndexAttributeExist[row][colume] = true;
			}
			else
			{
				m_arrayAAindexAttributeValue[row][colume] = 0.0;
				m_arrayIfIndexAttributeExist[row][colume] = false;
			}

			pstr = pstr1 + 1;
			pstr1 = strstr(pstr, "\t");
			colume++;
		}
		pstr1 = strstr(pstr, "\n");
		*pstr1 = '\0';
		if (pstr != "NA")                   //最后一列
		{
			m_arrayAAindexAttributeValue[row][colume] = atof(pstr);
			m_arrayIfIndexAttributeExist[row][colume] = true;
		}
		else
		{
			m_arrayAAindexAttributeValue[row][colume] = 0.0;
			m_arrayIfIndexAttributeExist[row][colume] = false;
		}
		row++;
		fgets(Buffer, BUFFERLENGTH, pFile);

	}
	//cout << "Have Loaded aaindexvalue!" << endl;
	return vStrTemp;
}
void PepAttributeWorker::mf_showAAindexValue(string path)
{
	ofstream ofile(path.c_str());
	for (int i = 0; i < 544; i++)
	{
		for (int j = 0; j < 20; j++)
			ofile << m_arrayAAindexAttributeValue[i][j] << "\t";
		ofile << endl;
	}


	ofile.close();
}

/*
Adds the sequence length attribute to each APEXPeptide within the APEXProtein
@param protein protein containing the APEXPeptide list for evaluation
*/
int PepAttributeWorker::mf_getLength(string seqChars)
{
	return (int)seqChars.length();
}

/*
Returns the mass of the amino acid frequency
@param seqChars input sequence
@return returns the mass of the sequence
*/
double PepAttributeWorker::mf_getMass(string seqChars)
{
	double mass = m_constants.m_dMASS_INITIAL;
	int aaIndex;
	for (aaIndex = 0; aaIndex<(int)seqChars.length(); aaIndex++)
	{
		if (seqChars[aaIndex] == 'A')
			mass += m_constants.m_dMASS_A;
		else if (seqChars[aaIndex] == 'C')
			mass += m_constants.m_dMASS_C;
		else if (seqChars[aaIndex] == 'D')
			mass += m_constants.m_dMASS_D;
		else if (seqChars[aaIndex] == 'E')
			mass += m_constants.m_dMASS_E;
		else if (seqChars[aaIndex] == 'F')
			mass += m_constants.m_dMASS_F;
		else if (seqChars[aaIndex] == 'G')
			mass += m_constants.m_dMASS_G;
		else if (seqChars[aaIndex] == 'H')
			mass += m_constants.m_dMASS_H;
		else if (seqChars[aaIndex] == 'I')
			mass += m_constants.m_dMASS_I;
		else if (seqChars[aaIndex] == 'K')
			mass += m_constants.m_dMASS_K;
		else if (seqChars[aaIndex] == 'L')
			mass += m_constants.m_dMASS_L;
		else if (seqChars[aaIndex] == 'M')
			mass += m_constants.m_dMASS_M;
		else if (seqChars[aaIndex] == 'N')
			mass += m_constants.m_dMASS_N;
		else if (seqChars[aaIndex] == 'P')
			mass += m_constants.m_dMASS_P;
		else if (seqChars[aaIndex] == 'Q')
			mass += m_constants.m_dMASS_Q;
		else if (seqChars[aaIndex] == 'R')
			mass += m_constants.m_dMASS_R;
		else if (seqChars[aaIndex] == 'S')
			mass += m_constants.m_dMASS_S;
		else if (seqChars[aaIndex] == 'T')
			mass += m_constants.m_dMASS_T;
		else if (seqChars[aaIndex] == 'V')
			mass += m_constants.m_dMASS_V;
		else if (seqChars[aaIndex] == 'W')
			mass += m_constants.m_dMASS_W;
		else if (seqChars[aaIndex] == 'Y')
			mass += m_constants.m_dMASS_Y;
	}
	return mass / 1000;
}
/**
* Returns the freqency of amino acid residues in the provided sequence.
* The order of the requencies is dictated by the String of aaSymbols
* @param seq sequence to analyze
* @param aaSymbols ordered set of aaSymbols
* @return a Vector of Float objects corresponding to the aa frequncies
* ordered by the aaSymbols String.
*/
vector<double> PepAttributeWorker::mf_getAAFrequencies(string seqChars, string aaSymbols)
{
	hash_map<string, double>  aaCounts;
	typedef pair<string, double> myPair;
	string aa;
	double count;
	vector<double> frequencies;

	for (int loc = 0; loc<(int)seqChars.length(); loc++)
	{
		aa = seqChars[loc];
		if (aaCounts.find(aa) != aaCounts.end())//can find it
		{
			count = (double)(aaCounts.at(aa)) + (double)1.;
			aaCounts[aa] = count;
			//			aaCounts.insert(myPair(aa,count));
		}
		else
		{
			aaCounts[aa] = (double)1.;
		}
	}

	for (int loc = 0; loc<(int)aaSymbols.length(); loc++)
	{
		aa = aaSymbols[loc];
		if (aaCounts.find(aa) != aaCounts.end())//can find it
		{
			frequencies.push_back((double)aaCounts.at(aa) / (double)seqChars.length());
		}
		else//cannot find it
		{
			frequencies.push_back((double)0);
		}
	}
	return frequencies;
}

int PepAttributeWorker::mf_GetMissingCleavageNumber(string Sequence)
{
	int count = 0;
	for (int i = 0; i < Sequence.size() - 1; i++)
	{
		if (((Sequence.at(i) == 'K') || (Sequence.at(i) == 'R')) && (Sequence.at(i+1)!='P'))  //add "P" by gzhq 20150723
			count++;
	}
	return count;
}

vector<int> PepAttributeWorker::mf_Sequence2Index(string sequence)
{
	vector<int > indexTemp;
	for (int i = 0; i < sequence.size(); i++)
	{
		switch (sequence[i])
		{
		case 'A':
			indexTemp.push_back(0);
			break;
		case 'C':
			indexTemp.push_back(1);
			break;
		case 'D':
			indexTemp.push_back(2);
			break;
		case'E':
			indexTemp.push_back(3);
			break;
		case 'F':
			indexTemp.push_back(4);
			break;
		case 'G':
			indexTemp.push_back(5);
			break;
		case 'H':
			indexTemp.push_back(6);
			break;
		case 'I':
			indexTemp.push_back(7);
			break;
		case 'K':
			indexTemp.push_back(8);
			break;
		case 'L':
			indexTemp.push_back(9);
			break;
		case 'M':
			indexTemp.push_back(10);
			break;
		case 'N':
			indexTemp.push_back(11);
			break;
		case 'P':
			indexTemp.push_back(12);
			break;
		case 'Q':
			indexTemp.push_back(13);
			break;
		case 'R':
			indexTemp.push_back(14);
			break;
		case 'S':
			indexTemp.push_back(15);
			break;
		case 'T':
			indexTemp.push_back(16);
			break;
		case 'V':
			indexTemp.push_back(17);
			break;
		case 'W':
			indexTemp.push_back(18);
			break;
		case 'Y':
			indexTemp.push_back(19);
			break;

		}
	}
	return indexTemp;
}
/**
* Computes net charge for a given sequence at a particular pH
* @param seqChars The amino acid sequence of the peptide
* @param pH The pH at which to evaluate net charge
* @return net charge at the given pH
*/

double PepAttributeWorker::mf_computeNetCharge(string seqChars, double pH)
{
	double ten2PKa;
	double ten2pH = pow((double)10, pH);
	//		double ten2pH = 
	double pKa;
	double netZ = (double)0;
	double sign = 1;
	bool contributes;

	ten2PKa = pow(10, m_constants.m_dPKa_N_TERM);
	//handle -COOH and -NH2			
	netZ = ten2PKa / (ten2pH + ten2PKa);
	ten2PKa = pow(10, m_constants.m_dPKa_C_TERM);
	netZ -= ten2pH / (ten2pH + ten2PKa);

	for (int aa = 0; aa < (int)seqChars.length(); aa++) {
		contributes = false;
		if (seqChars[aa] == 'R') {
			contributes = true;
			sign = 1;
			ten2PKa = pow(10, m_constants.m_dPKa_R);
		}
		else if (seqChars[aa] == 'K') {
			contributes = true;
			sign = 1;
			ten2PKa = pow(10, m_constants.m_dPKa_K);
		}
		else if (seqChars[aa] == 'H') {
			contributes = true;
			sign = 1;
			ten2PKa = pow(10, m_constants.m_dPKa_H);
		}
		else if (seqChars[aa] == 'C') {
			contributes = true;
			sign = -1;
			ten2PKa = pow(10, m_constants.m_dPKa_C);
		}
		else if (seqChars[aa] == 'D') {
			contributes = true;
			sign = -1;
			ten2PKa = pow(10, m_constants.m_dPKa_D);
		}
		else if (seqChars[aa] == 'E') {
			contributes = true;
			sign = -1;
			ten2PKa = pow(10, m_constants.m_dPKa_E);
		}
		else if (seqChars[aa] == 'Y') {
			contributes = true;
			sign = -1;
			ten2PKa = pow(10, m_constants.m_dPKa_Y);
		}

		if (contributes) {
			if (sign>0)
				netZ += sign*(ten2PKa / (ten2pH + ten2PKa));
			else
				netZ += sign*(ten2pH / (ten2pH + ten2PKa));
		}
	}
	return netZ;
}


int PepAttributeWorker::mf_getSumNeutralResidues(string seqChars)
{
	int neutralCount = 0;	//neutral residue count
	for (int aaIndex = 0; aaIndex < (int)seqChars.length(); aaIndex++)
	{
		if (seqChars[aaIndex] == 'E' || seqChars[aaIndex] == 'D' || seqChars[aaIndex] == 'R' || seqChars[aaIndex] == 'K' || seqChars[aaIndex] == 'H')// include H
			continue;
		neutralCount++;
	}
	return neutralCount;
}

int PepAttributeWorker::mf_getSumBasicResidues(string seqChars)
{
	int basicCount = 0;//basic residue count  R K H
	for (int aaIndex = 0; aaIndex<(int)seqChars.length(); aaIndex++)
	{
		if (seqChars[aaIndex] == 'R' || seqChars[aaIndex] == 'K' || seqChars[aaIndex] == 'H')
			basicCount++;
	}
	return basicCount;
}

//below 2 methods are added by Lei Wang on 12/09/2013, revised @ 16:06 Jan 7, 2014
//according to http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
double PepAttributeWorker::mf_computeProportionLargeSizedResidues(string seqChars)
{
	int count = 0;
	for (int index = 0; index<(int)seqChars.length(); index++)
	{
		//V, H, E, Q belong to the Medium Sized Residues Group
		if (seqChars[index] == 'M' || seqChars[index] == 'I' || seqChars[index] == 'L' || seqChars[index] == 'K' || seqChars[index] == 'R' || seqChars[index] == 'F' || seqChars[index] == 'Y' || seqChars[index] == 'W')
			count++;
	}
	return (double)count / (double)seqChars.length();
}

double PepAttributeWorker::mf_computeProportionSmallSizedResidues(string seqChars)
{
	int count = 0;
	for (int index = 0; index<(int)seqChars.length(); index++)
	{
		//V, H, E, Q belong to the Medium Sized Residues Group
		if (seqChars[index] == 'C' || seqChars[index] == 'P' || seqChars[index] == 'T' || seqChars[index] == 'D' || seqChars[index] == 'N' || seqChars[index] == 'A' || seqChars[index] == 'G' || seqChars[index] == 'S')
			count++;
	}
	return (double)count / (double)seqChars.length();

}
/*----------------------- below are properties from Haixu Tang--------------------------*/

double PepAttributeWorker::mf_getMass2LengthRatio(string seqChars)
{
	return mf_getMass(seqChars) / (double)mf_getLength(seqChars);
}

double PepAttributeWorker::mf_computeSequenceComplexity(vector<double> aaFrequencies)
{
	vector<double>::iterator iter;
	double entropy = 0;// the entropy on the informatics level
	double pi;//the temp variable, indicates the probability or the frequency
	for (iter = aaFrequencies.begin(); iter != aaFrequencies.end(); iter++)
	{
		//define that 0log=0
		if ((*iter) == 0)
			continue;
		pi = (*iter);
		entropy += pi * log(pi) / log(2.0);
	}
	if (entropy == 0)// to prevent the -0 from happening
		return 0;
	entropy = -entropy;
	return entropy*1000;  //change *1000 by gzhq 20150723
}

double PepAttributeWorker::mf_getVihinenFlexibility(string seqChars)
{
	double result = 0.0;//this is the one that should be returned to whichever calls the function.

	string const AminoAcids = "ACDEFGHIKLMNPQRSTVWYXBZ";
	double const FlexScales[] = { 0.984, 0.906, 1.068, 1.094, 0.915, 1.031, 0.950, 0.927, 1.102, 0.935, 0.952, 1.048, 1.049, 1.037, 1.008, 1.046, 0.997, 0.931, 0.904, 0.929, 0.9906, 1.068, 1.094 };
	double const window[] = { 0.25, 0.4375, 0.625, 0.8125, 1, 0.8125, 0.625, 0.4375, 0.25 };
	double const windowSum = 5.25;//this is the sum of all elements in array 'window'
	double windSum;//this is the sum  of all elements in array 'wind'
	double proportion;//this is windowSum/windSum.
	int sequenceLength;//this is the length of the sequence, i.e. seqChars
	int strIndex, index;//this is a index of a string
	int begin, end;//this is the start/end index of the sub of window.
	int windowSize = 9;
	//	vector<double> prediction;//this is the result that should be returned.
	vector<double> subFlexScales;//this contains the sub of the FlexScales

	//	prediction.clear();
	sequenceLength = seqChars.length();
	for (strIndex = 0; strIndex < sequenceLength; strIndex++)
	{
		//cout<<"strIndex: "<<strIndex<<endl;
		subFlexScales.clear();
		int Max = max(0, strIndex - 4), Min = min(strIndex + 4, sequenceLength - 1);
		for (index = Max; index <= Min; index++)
			subFlexScales.push_back(FlexScales[AminoAcids.find(seqChars[index])]);

		if (strIndex == 0)
		{
			begin = 4;
			end = windowSize - 1;
		}
		else if (strIndex == 1)
		{
			begin = 3;
			end = windowSize - 1;
		}

		else if (strIndex == 2)
		{
			begin = 2;
			end = windowSize - 1;
		}

		else if (strIndex == 3)
		{
			begin = 1;
			end = windowSize - 1;
		}

		else if (strIndex == sequenceLength - 4)
		{
			begin = 0;
			end = windowSize - 2;
		}
		else if (strIndex == sequenceLength - 3)
		{
			begin = 0;
			end = windowSize - 3;
		}
		else if (strIndex == sequenceLength - 2)
		{
			begin = 0;
			end = windowSize - 4;
		}
		else if (strIndex == sequenceLength - 1)
		{
			begin = 0;
			end = windowSize - 5;
		}
		else
		{
			begin = 0;
			end = windowSize - 1;
		}

		if ((end - begin + 1) != subFlexScales.size())
		{
			cout << "Vector wind's size is not equal to Vector subFlexScales's size" << endl;
			cout << "Vector wind size: " << end - begin + 1 << endl;
			cout << "Vector subFlexScales size: " << subFlexScales.size() << endl;
			system("pause");
		}
		windSum = 0;
		for (index = begin; index <= end; index++)
			windSum += window[index];
		proportion = windowSum / windSum;

		double temp = 0;
		for (index = 0; index<(int)subFlexScales.size(); index++)
		{
			temp += subFlexScales.at(index) * window[index + begin] * proportion;
		}
		//prediction.push_back(temp);
		//cout<<temp<<endl; //this is to test the intermideate result of the function
		result += temp;
	}//for
	return result / seqChars.length();
}
double PepAttributeWorker::mf_getHydrophobicMoment(string seqChars, double window, double angle)
{
	/*this function calculates the hydrophobic moment for a protein sequence,
	Input:
	protein sequence, i.e. seqChars
	window: window size (a number)
	angle: rotation angle in degrees (100 degrees for helix)
	Output:
	the mean of the result from vector of hydrophobic moments
	*/
	double result = 0;//this is the one that should be returned when this function is called
	string const AminoAcids = "ACDEFGHIKLMNPQRSTVWY";
	//vector<double> hm;//this is the only result that should be returned.
	vector<int> aaIndex;//this is a vector of indicies corresponding to particular aa in the whole sequence.
	double hydrophobicIndicies[] = { 0.25, 0.04, -0.72, -0.62, 0.61, 0.16, -0.40, 0.73, -1.10, 0.53, 0.26, -0.64, -0.07, -0.69, -1.76, -0.26, -0.18, 0.54, 0.37, 0.02 };
	int halfWindow, sequenceLength, index, index2;
	int Max, Min;
	double hmSin, hmCos;
	halfWindow = floor(window / 2.);
	angle = angle * m_constants.m_dPi / 180;//convert an angle into radians
	sequenceLength = seqChars.length();
	for (index = 0; index<sequenceLength; index++)
		aaIndex.push_back(AminoAcids.find(seqChars[index]));

	//calculate Hydrophobic Indicies over the input sequence.
	for (index = 0; index<sequenceLength; index++)
	{
		hmSin = hmCos = 0;
		Max = max(0, index - halfWindow);
		Min = min(index + halfWindow, sequenceLength - 1);
		for (index2 = Max; index2 <= Min; index2++)
		{
			hmSin += sin(angle * index2) * hydrophobicIndicies[aaIndex.at(index2)];
			hmCos += cos(angle * index2) * hydrophobicIndicies[aaIndex.at(index2)];
		}

		//hydrophobic moment is normalized with the window length
		double temp = Min - Max + 1;
		double temp2;
		if (temp< window)
		{
			temp2 = sqrt(hmSin*hmSin + hmCos*hmCos) / temp;
			//hm.push_back(temp2);
			result += temp2;
			//			cout<<temp2<<endl;
		}
		else
		{
			temp2 = sqrt(hmSin*hmSin + hmCos*hmCos) / window;
			//hm.push_back(temp2);
			result += temp2;
			//			cout<<temp2<<endl;
		}
	}
	return result / seqChars.length();
}

//these functions are the ones that I write to simulate the corresponding functions in Matlab used by the functions above
/*
Declaration: Actually, this filter() function is from the website: http://mechatronics.ece.usu.edu/yqchen/filter.c/FILTER.C by YangQuan Chen
the inputs:
ord means order, i.e. (size of array a) - 1.
a means array a
b means array b.
np is the size of array x
x is the variables
y is the result

Comment:
１. this filter() assumes that  the input parameters have the requirement: np >= order +１
2.  the original filter() function has an error in the sentence: 　		for (i=ord+1;i<np+1;i++), actually it should be like this:
for (i=ord+1;i<np;i++), if np +1, then access violation for the arrays x[] and y[] whose size length is np, i.e. index range: 0 - (np-1)
*/

void PepAttributeWorker::mf_filter(int ord, double *a, double *b, int np, double *x, double *y)
{
	int i, j;
	y[0] = b[0] * x[0];
	for (i = 1; i<ord + 1; i++)
	{
		y[i] = 0.0;
		for (j = 0; j<i + 1; j++)
			y[i] = y[i] + b[j] * x[i - j];
		for (j = 0; j<i; j++)
			y[i] = y[i] - a[j + 1] * y[i - j - 1];
	}
	/* end of initial part */
	for (i = ord + 1; i<np; i++)
	{
		y[i] = 0.0;
		for (j = 0; j<ord + 1; j++)
			y[i] = y[i] + b[j] * x[i - j];
		for (j = 0; j<ord; j++)
			y[i] = y[i] - a[j + 1] * y[i - j - 1];
	}
	/* end of filter */
}

void PepAttributeWorker::mf_filter_twosided(double **matrix_protein, int matrixRow, int matrixColumn, double * filt_onesided, int filtOnesidedSize)
{
	int index, indexMatirxColumn;
	double *denominator = new double[filtOnesidedSize];

	/*
	WhH DO WE NEED THE BIGGER ONE OF filterOnesidedSize and matrixRow?
	to make sure that the size of array one[] and scale[] is the bigger one between filtOnesidedSize and matrixRow.
	we only need matrixRow-sized one[] and scale[], but we have to allocate at least (ord+1)-sized memory to avoid the access violation in the filter() function.
	*/
	int biggerLength = (filtOnesidedSize > matrixRow ? filtOnesidedSize : matrixRow);
	double *one = new double[biggerLength];
	double *scale = new double[biggerLength];

	for (index = 0; index<filtOnesidedSize; index++)
		denominator[index] = 0;
	denominator[0] = 1;

	for (index = 0; index<matrixRow; index++)
		one[index] = 1;

	//commented @ 14:41 Feb 14, 2014
	mf_filter((filtOnesidedSize - 1), denominator, filt_onesided, matrixRow, one, scale);


	//scale = (scale1 + scale2)/2
	for (index = 0; index<matrixRow / 2; index++)
	{
		scale[index] = (scale[index] + scale[matrixRow - 1 - index]) / 2.;
		scale[matrixRow - 1 - index] = scale[index];
	}

	double *data = new double[biggerLength];
	double *fdata = new double[biggerLength];
	double *fdata2 = new double[biggerLength];
	for (indexMatirxColumn = 0; indexMatirxColumn < matrixColumn; indexMatirxColumn++)
	{

		for (index = 0; index<matrixRow; index++)
			data[index] = matrix_protein[index][indexMatirxColumn];

		mf_filter((filtOnesidedSize - 1), denominator, filt_onesided, matrixRow, data, fdata);

		//now flip the data
		for (index = 0; index<matrixRow / 2; index++)
			swap(data[index], data[matrixRow - 1 - index]);

		mf_filter((filtOnesidedSize - 1), denominator, filt_onesided, matrixRow, data, fdata2);

		//scale = (fdata1 + fdata2)/2
		for (index = 0; index<matrixRow; index++)
			fdata[index] = (fdata[index] + fdata2[matrixRow - 1 - index]) / 2.;

		for (index = 0; index<matrixRow; index++)
			matrix_protein[index][indexMatirxColumn] = fdata[index] / scale[index];
	}

	//now delete all unnecessary arrays to release memory.
	delete[] fdata2;
	delete[] fdata;
	delete[] data;
	delete[] denominator;
	delete[] one;
	delete[] scale;
}

vector<double> PepAttributeWorker::mf_getVL2Disorder(string seqChars)
{
	int length = seqChars.length();
	int index, indexColumn;
	double **pfilt;
	pfilt = new double *[length];
	for (index = 0; index<length; index++)
		pfilt[index] = new double[4];
	//initialize pfilt[][]
	for (index = 0; index<length; index++)
		for (indexColumn = 0; indexColumn<4; indexColumn++)
			pfilt[index][indexColumn] = 0;

	mf_predictLinear(seqChars, pfilt, length, 4);

	double result = 0;
	vector<double> vec;
	vec.clear();
	for (indexColumn = 0; indexColumn<4; indexColumn++)
	{
		result = 0;
		for (index = 0; index<length; index++)
		{
			result += pfilt[index][indexColumn];
			//			cout<<pfilt[index][indexColumn]<<" ";
		}
		//		cout<<endl;
		result /= seqChars.size();
		vec.push_back(result);
	}

	//delete arrays
	for (index = 0; index<length; index++)
		delete[] pfilt[index];
	delete[] pfilt;

	return vec;
}

void PepAttributeWorker::mf_predictLinear(string seqChars, double **pfilt, int pfiltRow, int pfiltColumn)
{
	double modelsBeta[4][21] =
	{
		-25.360122772689530, -0.696903041342054, 1.759821058633487, 0.088936843652238, 1.816910774835402, 2.580379123489578, 1.199846997885154, -0.095699661003024, 1.098553523477732, 2.636616803052733, -1.665418976414453, 1.352861445754487, -0.580508778536282, 2.745993680647693, -1.085823401526897, 0.740717492440777, 0.533128494727820, -4.586324163702638, 0.739585357799028, 24.425480517325200, 0.230030178184243,
		-14.300860682065565, 1.434701221055386, -1.028707033124128, -3.013912732313829, -2.376741271298584, 1.273834113496491, -0.964875311764001, -1.094035314236643, -1.044061215564152, 2.773847613995727, -0.553544312680591, -0.152639027401581, 0.678673838499450, 1.016487071159430, -0.165345043318697, 0.014607300630189, 0.925090573506023, 2.515930202256339, 0.974841031279123, 15.431464279382766, -0.221008947365943,
		-30.337058421955790, 0.854440689066855, 5.339030779290762, 1.010144293991523, 3.381492895434795, 1.206914784027385, 2.371572902795911, -0.850297513246717, 3.987712522463211, 4.790343477250348, 1.081410650555920, -0.169321385683294, 1.554378338506579, 3.520055521634506, 2.432348578718518, 1.738359018604442, 3.223434358612503, 1.679678304063250, -1.456869987959105, 30.284784500965750, -0.354270870941984,
		-25.970718651350940, 0.683574136190159, 3.082311205666289, -0.370673832896873, 1.626657146800731, 2.486017822353065, 1.360178852772070, -0.732837806618001, 1.769823119646163, 4.086509434704226, -0.328578025097438, 0.252295819229421, 0.431583241825593, 2.009388278166486, 0.456480475170492, 0.666658611714632, 1.805648436817136, -0.545838679329645, -0.127664204931993, 26.489721415658103, -0.257148340210382
	};

	int index, indexColumn;
	const int W_IN = 21, W_OUT = 21;
	double * filt_onesided = new double[W_IN];
	for (index = 0; index<W_IN; index++)
		filt_onesided[index] = 1.;
	filt_onesided[0] = 0.5;

	//protein_flat_filter()
	int sequenceLength = seqChars.length();
	double **data;//[length x 26 ]
	data = new double *[sequenceLength];
	for (index = 0; index<sequenceLength; index++)
		data[index] = new double[26];

	//initialize data[][]
	for (index = 0; index<sequenceLength; index++)
		for (indexColumn = 0; indexColumn<26; indexColumn++)
			data[index][indexColumn] = 0;

	mf_proteinFlatFilter(data, sequenceLength, 26, seqChars, filt_onesided, W_IN);//to change 'data'

	//commented @ 18:12 Feb 15, 2014
	////test data[][]
	//cout<<"This is to test data[][] after protein_flat_filter():"<<endl;
	//for(index=0; index< sequenceLength; index++)
	//{
	//	cout<<"Line "<<index+1<<": ";
	//	for(indexColumn=0; indexColumn<26; indexColumn++)
	//		cout<<data[index][indexColumn]<<" ";
	//	cout<<endl;
	//}		
	//cout<<"This is the end of the data[][]"<<endl;
	//system("pause");

	//to calculate entropy
	double entropy;
	double temp;
	for (index = 0; index<sequenceLength; index++)
	{
		entropy = 0;
		for (indexColumn = 0; indexColumn<20; indexColumn++)
		{
			temp = data[index][indexColumn];
			if (temp == 0)
				continue;
			entropy += temp * log(1. / temp) / log(2.);
		}
		data[index][21] = entropy;

		//This is to test data[i][21]
		//cout<<"Line "<<index +1 <<": "<<data[index][21]<<endl;
	}

	//system("pause");

	for (index = 0; index<W_OUT; index++)
		filt_onesided[index] = 1;
	filt_onesided[0] = 0.5;

	int lmod = 4;//length(mdels)
	int Attributes[20] = { 0, 1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };
	int attIndex;//index of the attributes
	int lenIndex;// index of the sequenceLength
	double sum;
	double *p = new double[sequenceLength];//correspondingly, p(: ; j)
	for (index = 0; index< lmod; index++)
	{
		//		cout<<"Line "<<index<<": "<<endl;
		for (lenIndex = 0; lenIndex < sequenceLength; lenIndex++)
		{
			sum = modelsBeta[index][0];
			for (attIndex = 0; attIndex<20; attIndex++)
				sum += data[lenIndex][Attributes[attIndex]] * modelsBeta[index][attIndex + 1];
			p[lenIndex] = sum;

			//This is to test p[]
			//cout<<lenIndex<<": "<<p[lenIndex]<<endl;
		}

		mf_filter_prot(p, pfilt, pfiltRow, index, filt_onesided, W_OUT);
	}

	//delete arrays
	delete[] p;
	delete[] filt_onesided;
	for (index = 0; index<sequenceLength; index++)
		delete[] data[index];
	delete[] data;
}
void PepAttributeWorker::mf_filter_prot(double *pred, double **pfilt, int pfiltRow, int pfiltColumn, double *filt_onesided, int filtOnesidedSize)
{
	/*
	Calculating what is the number of proteins in the database with chance of long disorder -- quoted from the original Matlab source code of P.P @ IUB
	*/
	int biggerLength = (filtOnesidedSize > pfiltRow ? filtOnesidedSize : pfiltRow);
	double *p = new double[biggerLength];
	int index;
	for (index = 0; index<pfiltRow; index++)
		p[index] = pred[index];
	mf_filter_error(pfilt, pfiltRow, pfiltColumn, p, pfiltRow, filt_onesided, filtOnesidedSize);

	//delete array
	delete[] p;
}
void PepAttributeWorker::mf_filter_error(double **pfilt, int pfiltRow, int pfiltColumn, double *data, int dataLength, double *filt_onesided, int filtOnesidedSize)
{
	/*
	Quoted from the original Matlab source code 26 Jun 1999
	--This function performs proper two-sided filtering of error to aleviate edge effects
	--True filter size is twice the size of one-sided filter!!
	*/
	int index;
	double *denominator = new double[filtOnesidedSize];
	int biggerLength = (filtOnesidedSize > dataLength ? filtOnesidedSize : dataLength);
	double * one = new double[biggerLength];
	double *scale = new double[biggerLength];
	for (index = 0; index<dataLength; index++)
		one[index] = 1;

	for (index = 0; index<filtOnesidedSize; index++)
		denominator[index] = 0;
	denominator[0] = 1;

	mf_filter((filtOnesidedSize - 1), denominator, filt_onesided, dataLength, one, scale);

	//scale = (scale1 + scale2)/2
	for (index = 0; index<dataLength / 2; index++)
	{
		scale[index] = (scale[index] + scale[dataLength - 1 - index]) / 2.;
		scale[dataLength - 1 - index] = scale[index];
	}

	double *fdata = new double[biggerLength];
	double *fdata2 = new double[biggerLength];

	mf_filter((filtOnesidedSize - 1), denominator, filt_onesided, dataLength, data, fdata);

	//now flip the data
	for (index = 0; index<dataLength / 2; index++)
		swap(data[index], data[dataLength - 1 - index]);

	mf_filter((filtOnesidedSize - 1), denominator, filt_onesided, dataLength, data, fdata2);

	//scale = (fdata1 + fdata2)/2
	for (index = 0; index<dataLength; index++)
		fdata[index] = (fdata[index] + fdata2[dataLength - 1 - index]) / 2.;

	for (index = 0; index<pfiltRow; index++)
		pfilt[index][pfiltColumn] = fdata[index] / scale[index];

	//delete arrays
	delete[] denominator;
	delete[] one;
	delete[]scale;
	delete[] fdata;
	delete[] fdata2;
}

void PepAttributeWorker::mf_proteinFlatFilter(double **matrix_protein, int matrixRow, int matrixColum, string sequence, double * filt_onesided, int filtOnesidedSize)
{
	double FHC[4][20] =
	{ 65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 80, 81, 82, 83, 84, 86, 87, 89,
	0.984000, 0.906000, 1.068000, 1.094000, 0.915000, 1.031000, 0.950000, 0.927000, 1.102000, 0.935000, 0.952000, 1.048000, 1.049000, 1.037000, 1.008000, 1.046000, 0.997000, 0.931000, 0.904000, 0.929000,
	1.800000, 2.500000, -3.500000, -3.500000, 2.800000, -0.400000, -3.200000, 4.500000, -3.900000, 3.800000, 1.900000, -3.500000, -1.600000, -3.500000, -4.500000, -0.800000, -0.700000, 4.200000, -0.900000, -1.300000,
	4.970000, 7.490000, 3.990000, 3.970000, 5.990000, 4.720000, 6.110000, 5.920000, 4.470000, 6.000000, 4.890000, 5.020000, 3.710000, 4.820000, 4.750000, 4.310000, 5.000000, 5.600000, 5.810000, 6.130000
	};

	int index, indexSequenceASCII;
	int sequenceLength = sequence.length();
	vector<int> sequenceASCII;
	vector<int> sequenceUniqueASCII;
	vector<int>::iterator  iter;
	for (index = 0; index<sequenceLength; index++)
	{
		int temp = (int)sequence[index];
		sequenceASCII.push_back(temp);
		iter = std::find(sequenceUniqueASCII.begin(), sequenceUniqueASCII.end(), temp);
		if (iter == sequenceUniqueASCII.end())//if not found
			sequenceUniqueASCII.push_back(temp);
	}

	iter = sequenceUniqueASCII.begin();
	for (; iter != sequenceUniqueASCII.end(); iter++)
	{
		for (index = 0; index<20; index++)
			if (FHC[0][index] == (*iter))
				break;

		for (indexSequenceASCII = 0; indexSequenceASCII < (int)sequenceASCII.size(); indexSequenceASCII++)
		{
			if (sequenceASCII.at(indexSequenceASCII) != (*iter))
				continue;
			matrix_protein[indexSequenceASCII][index] = 1;
			matrix_protein[indexSequenceASCII][20] = FHC[1][index];
			matrix_protein[indexSequenceASCII][21] = FHC[2][index];
			matrix_protein[indexSequenceASCII][22] = FHC[3][index];
		}
	}

	for (index = 0; index< matrixRow; index++)
	{
		matrix_protein[index][23] = 0;
		matrix_protein[index][24] = 1;
		matrix_protein[index][25] = (index + 1);//index or index+1 ?
	}

	//cout<<"the address of [] filt_onesided in function() proteinFlatFilter:"<<endl;
	//cout<<filt_onesided<<endl;

	//filter_two sided() function 
	mf_filter_twosided(matrix_protein, sequenceLength, 23, filt_onesided, filtOnesidedSize);

	////test codes
	//cout<<"This is to test protein_flat_filter():"<<endl;
	//int i;
	//for(index=0; index<matrixRow; index++)
	//{
	//	for(i=0; i<matrixColum;i++)
	//		cout<<matrix_protein[index][i]<<" ";
	//	cout<<endl;
	//}
	//system("pause");
	//	cout<<"This is the end of the method -- proteinFlatFilter"<<endl;
}

double PepAttributeWorker::mf_getBFactorPrediction(string seqChars, int Win, int Wout)
{
	double result = 0;//this is the result that should be returned when this function is called
	//load linpred_ns; load FHC
	double const meanv[20] =
	{
		0.082968313140730, 0.011736564150357, 0.067809878844363, 0.064588381484934, 0.033886300093196, 0.088146940043496, 0.022098477788133, 0.045552966759862, 0.064153463808637, 0.074777881329608, 0.020540540540540, 0.051916744330536, 0.053502640571605, 0.039798073936004, 0.044307238272754, 0.067135756446103, 0.061397949673812, 0.061351351351352, 1.006507096924509, 2.053169508147902
	};
	double const stdv[20] =
	{
		0.130331719672796, 0.049309161430501, 0.112881088280996, 0.110767901251959, 0.081279119138629, 0.128216021352435, 0.066946723059315, 0.094605357612097, 0.112588332693937, 0.119775107696711, 0.065689847423736, 0.101495390561351, 0.101360742817996, 0.089508175324723, 0.093472444265776, 0.116208004856971, 0.110493688968059, 0.109700608158085, 0.026172472801445, 0.309130677346108
	};
	double const resultsBeta[21] =
	{
		-0.127150878905927, 0.051137322783183, -0.284651567902609, -0.430461262094386, -0.021608943274078, -0.154792867992833, 0.049308110391981, -0.104625789276418, -0.467049531281995, 0.009179581249833, -0.005886075350205, -0.219643270163331, -0.176474696960102, -0.220613294186365, -0.150917721214844, -0.253207400169961, -0.121781926993112, -0.056917028712712, 1.137733284555524, -0.185855943167836, 0.0003732630317847297
	};
	double FHC[4][20] =
	{ 65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 80, 81, 82, 83, 84, 86, 87, 89,
	0.984000, 0.906000, 1.068000, 1.094000, 0.915000, 1.031000, 0.950000, 0.927000, 1.102000, 0.935000, 0.952000, 1.048000, 1.049000, 1.037000, 1.008000, 1.046000, 0.997000, 0.931000, 0.904000, 0.929000,
	1.800000, 2.500000, -3.500000, -3.500000, 2.800000, -0.400000, -3.200000, 4.500000, -3.900000, 3.800000, 1.900000, -3.500000, -1.600000, -3.500000, -4.500000, -0.800000, -0.700000, 4.200000, -0.900000, -1.300000,
	4.970000, 7.490000, 3.990000, 3.970000, 5.990000, 4.720000, 6.110000, 5.920000, 4.470000, 6.000000, 4.890000, 5.020000, 3.710000, 4.820000, 4.750000, 4.310000, 5.000000, 5.600000, 5.810000, 6.130000
	};

	int index, indexColumn;
	int sequenceLength = seqChars.length();
	double **xtemp;
	xtemp = new double *[sequenceLength];
	for (index = 0; index<sequenceLength; index++)
		xtemp[index] = new double[26];

	//initialize x[][]
	for (index = 0; index<sequenceLength; index++)
		for (indexColumn = 0; indexColumn<26; indexColumn++)
			xtemp[index][indexColumn] = 0;

	mf_makeAttribute(xtemp, sequenceLength, 26, seqChars, Win);

	////	test the result xtemp vector of makeAttribute()
	//	cout<<"xtemp vector from makeAttribute(): "<<endl;
	//	for(index=0;index<sequenceLength; index++)
	//	{
	//		for(indexColumn=0; indexColumn<26;indexColumn++)
	//			cout<<x[index][indexColumn]<<" ";
	//		cout<<endl;
	//	}
	//	system("pause");

	//cout<<"This is to test normalized x vector: "<<endl;
	int indexMS;//index of meanv[] and stdv[]
	for (index = 0; index<sequenceLength; index++)//normalize
	{
		indexMS = 0;
		for (indexColumn = 0; indexColumn< 24; indexColumn++)
		{
			if (indexColumn == 18 || indexColumn == 19 || indexColumn == 21 || indexColumn == 22)
				continue;
			xtemp[index][indexColumn] = (xtemp[index][indexColumn] - meanv[indexMS]) / stdv[indexMS++];
			//	cout<<x[index][indexColumn]<<" ";
		}
		//		cout<<endl;
	}

	//test raw_prediction[]
	//	cout<<"This is to test raw_prediction[]: "<<endl;
	double * raw_prediction = new double[sequenceLength];
	int indexResultsBeta;
	double sum;
	for (index = 0; index<sequenceLength; index++)
	{
		indexResultsBeta = 0;
		sum = 0;
		for (indexColumn = 0; indexColumn< 24; indexColumn++)
		{
			if (indexColumn == 18 || indexColumn == 19 || indexColumn == 21 || indexColumn == 22)
				continue;
			sum += xtemp[index][indexColumn] * resultsBeta[indexResultsBeta++];
		}
		sum += resultsBeta[indexResultsBeta];
		raw_prediction[index] = 1. / (1 + exp(-sum));
	}//for

	//moving average
	double * prediction = new double[sequenceLength];
	for (index = 0; index<sequenceLength; index++)
		prediction[index] = 0;

	/*
	Actually, after the moving_average() function, prediction[] is filled with the same result which is also the output of
	the corresponding Matlab predictBfactors() function.
	*/

	mf_moving_average(raw_prediction, Wout, prediction, sequenceLength);

	//cout<<"prediction vector: "<<endl;
	for (index = 0; index<sequenceLength; index++)
	{
		result += prediction[index];
		//cout<<prediction[index]<<endl;
	}
	//delete array
	delete[] prediction;
	delete[] raw_prediction;
	for (index = 0; index<sequenceLength; index++)
		delete[]xtemp[index];
	delete[]xtemp;

	return result / seqChars.length();
}
void  PepAttributeWorker::mf_moving_average(double *x, int window, double * y, int length)//length(x) == length(y)
{
	/*
	According to moving_average.m by Predrag Radivojac, if a sequence length is comparable to the window size (or shorter)
	we do quadratic loop.
	*/

	//the assumption is that both window and the sequence x are short

	int index, xIndex;
	double sum;
	int maximum, minimum;
	if (window >= length)
	{
		for (index = 0; index<length; index++)
		{
			//y[index] = 0;
			maximum = max(index - (int)floor(window / 2.), 0);
			minimum = min(index + (int)floor(window / 2.), length - 1);
			sum = 0;
			for (xIndex = maximum; xIndex <= minimum; xIndex++)
				sum += x[xIndex];
			y[index] = sum / (minimum - maximum + 1);
		}
		return;
	}

	//otherwise, we do linear loop
	int one_sided_window = floor(window / 2.) + 1;//it should be set up to zero vector, however, i iintentionally not do it.
	double t = 0;//t is a running sum
	for (index = 0; index < one_sided_window; index++)
		t += x[index];

	//first, we expand the window at the left end
	for (index = 0; index<one_sided_window - 1; index++)
	{
		y[index] = t / (one_sided_window + index);
		t += x[one_sided_window + index];
		//		cout<<y[index]<<endl;
	}

	//second, we work with the full-sized window ('middle' of the sequence x)
	for (index = one_sided_window - 1; index < length - one_sided_window; index++)
	{
		y[index] = t / window;
		t = t - x[index - one_sided_window + 1] + x[one_sided_window + index];
	}

	//third, we collapse the window at the right end
	int j = 0;
	for (index = length - one_sided_window; index < length; index++)
	{
		y[index] = t / (window - j);
		t = t - x[index - one_sided_window + 1];
		j++;
	}
}

void  PepAttributeWorker::mf_makeAttribute(double **data, int dataRow, int dataColumn, string seqChars, int W_IN)
{
	double * filt_onesided = new double[W_IN];
	int index, indexColumn;

	for (index = 0; index<W_IN; index++)
		filt_onesided[index] = 1.;
	filt_onesided[0] = 0.5;

	mf_proteinFlatFilter(data, dataRow, dataColumn, seqChars, filt_onesided, W_IN);//to change 'data'

	//to calculate entropy
	double entropy;
	for (index = 0; index<dataRow; index++)
	{
		entropy = 0;
		for (indexColumn = 0; indexColumn<20; indexColumn++)
		{
			double temp = data[index][indexColumn];
			if (temp == 0)
				continue;
			entropy += temp * log(1. / temp) / log(2.);
		}
		data[index][23] = entropy;
	}

	//added @ 15:29 Feb 14, 2014
	delete[] filt_onesided;
	////test codes
	//cout<<"This is to test makeAttribute():"<<endl;
	//int i;
	//for(index=0; index<dataRow; index++)
	//{
	//	for(i=0; i<dataColumn;i++)
	//		cout<<data[index][i]<<" ";
	//	cout<<endl;
	//}
	//system("pause");

}
/*--------------Supported Vector Machine---------------------------*/
int PepAttributeWorker::mf_getNumberNonPolarHydrophobicResidues(string seqChars)
{
	//nonpolar == hydrophobic
	int sum = 0;
	for (int aaIndex = 0; aaIndex < (int)seqChars.length(); aaIndex++)
	{
		//F A L M I W P V
		if (seqChars[aaIndex] == 'F')
			sum++;
		else if (seqChars[aaIndex] == 'A')
			sum++;
		else if (seqChars[aaIndex] == 'L')
			sum++;
		else if (seqChars[aaIndex] == 'M')
			sum++;
		else if (seqChars[aaIndex] == 'I')
			sum++;
		else if (seqChars[aaIndex] == 'W')
			sum++;
		else if (seqChars[aaIndex] == 'P')
			sum++;
		else if (seqChars[aaIndex] == 'V')
			sum++;
	}
	return sum;
}

int PepAttributeWorker::mf_getNumberPolarHydrophilicResidues(string seqChars)
{
	//polar = hydrophilic
	int sum = 0;
	for (int aaIndex = 0; aaIndex<(int)seqChars.length(); aaIndex++)
	{
		//acid residues
		//D E
		if (seqChars[aaIndex] == 'D')
			sum++;
		else if (seqChars[aaIndex] == 'E')
			sum++;

		//basic residues
		// R K H
		else if (seqChars[aaIndex] == 'R')
			sum++;
		else if (seqChars[aaIndex] == 'K')
			sum++;
		else if (seqChars[aaIndex] == 'H')
			sum++;

		//polar uncharged residues
		//C G Q N S Y T
		else if (seqChars[aaIndex] == 'C')
			sum++;
		else if (seqChars[aaIndex] == 'G')
			sum++;
		else if (seqChars[aaIndex] == 'Q')
			sum++;
		else if (seqChars[aaIndex] == 'N')
			sum++;
		else if (seqChars[aaIndex] == 'S')
			sum++;
		else if (seqChars[aaIndex] == 'Y')
			sum++;
		else if (seqChars[aaIndex] == 'T')
			sum++;
	}
	return sum;
}

int PepAttributeWorker::mf_getNumberUnchargedPolarHydrophilicResidues(string seqChars)
{
	//uncharged polar hydrophilic residues
	//C G Q N S Y T
	int sum = 0;
	for (int aaIndex = 0; aaIndex<(int)seqChars.length(); aaIndex++)
	{
		if (seqChars[aaIndex] == 'C')
			sum++;
		else if (seqChars[aaIndex] == 'G')
			sum++;
		else if (seqChars[aaIndex] == 'Q')
			sum++;
		else if (seqChars[aaIndex] == 'N')
			sum++;
		else if (seqChars[aaIndex] == 'S')
			sum++;
		else if (seqChars[aaIndex] == 'Y')
			sum++;
		else if (seqChars[aaIndex] == 'T')
			sum++;
	}
	return sum;
}

int PepAttributeWorker::mf_getNumberChargedPolarHydrophilicResidues(string seqChars)
{
	//charged polar hydrophilic residues.   Acidic residues +  Basic residues
	//D E R K H
	int sum = 0;
	for (int aaIndex = 0; aaIndex<(int)seqChars.length(); aaIndex++)
	{
		if (seqChars[aaIndex] == 'D')
			sum++;
		else if (seqChars[aaIndex] == 'E')
			sum++;
		else if (seqChars[aaIndex] == 'R')
			sum++;
		else if (seqChars[aaIndex] == 'K')
			sum++;
		else if (seqChars[aaIndex] == 'H')
			sum++;
	}
	return sum;
}


vector<double>  PepAttributeWorker::mf_GetAttributeFromSequence(string Sequence, string SequenceWithFlankingRegion)
{
	vector<double> Attributes;
	//Attributes Computed without FlankingRegion
	vector<double> aaFrequencies;
	vector<double>::iterator iterator;


	Attributes.push_back(((double)mf_getLength(Sequence)));
	Attributes.push_back(((double)mf_GetMissingCleavageNumber(Sequence)));
	Attributes.push_back(mf_getMass(Sequence));

	aaFrequencies = mf_getAAFrequencies(Sequence, "ACDEFGHIKLMNPQRSTVWY");
	iterator = aaFrequencies.begin();
	for (; iterator != aaFrequencies.end(); iterator++)
		Attributes.push_back((*iterator));

	//从sequence算544种特征

	vector<int> IndexOfSequence;
	IndexOfSequence = mf_Sequence2Index(Sequence);
	double FeatureValueTemp = 0.0;
	vector<int>::iterator iter;
	int count = 0;
	for (int row = 0; row < 544; row++)
	{
		count = 0;
		FeatureValueTemp = 0.0;
		for (iter = IndexOfSequence.begin(); iter != IndexOfSequence.end(); iter++)
		{
			if (m_arrayIfIndexAttributeExist[row][*iter])
			{
				FeatureValueTemp += m_arrayAAindexAttributeValue[row][*iter];
				count++;
			}

		}

		FeatureValueTemp = FeatureValueTemp / count;
		Attributes.push_back(FeatureValueTemp);
	}

	Attributes.push_back(mf_computeNetCharge(Sequence, 7));

	/*----------------------From  SVM Paper, in which all of these 4 functions have been implemented.------------*/
	Attributes.push_back((double)mf_getNumberNonPolarHydrophobicResidues(Sequence));
	Attributes.push_back((double)mf_getNumberPolarHydrophilicResidues(Sequence));
	Attributes.push_back((double)mf_getNumberUnchargedPolarHydrophilicResidues(Sequence));
	Attributes.push_back((double)mf_getNumberChargedPolarHydrophilicResidues(Sequence));

	Attributes.push_back((double)mf_getSumBasicResidues(Sequence));
	Attributes.push_back((double)mf_getSumNeutralResidues(Sequence));
	Attributes.push_back(mf_computeProportionLargeSizedResidues(Sequence));
	Attributes.push_back(mf_computeProportionSmallSizedResidues(Sequence));

	/*------------Attributes From Haixu Tang------------*/
	Attributes.push_back(mf_getMass2LengthRatio(Sequence));
	Attributes.push_back(mf_computeSequenceComplexity(aaFrequencies));

	//Attributes Computed with FlankingRegion
	//vector<double> Attributes;//this vector contains all attriibutes of the peptide sequence, i.e. SequenceWithFlankingRegion
	vector<double> VL2;
	vector<double>::iterator VL2iter;
	Attributes.push_back(mf_getVihinenFlexibility(SequenceWithFlankingRegion));
	Attributes.push_back(mf_getHydrophobicMoment(SequenceWithFlankingRegion, 11, 100));
	Attributes.push_back(mf_getHydrophobicMoment(SequenceWithFlankingRegion, 11, 160));
	Attributes.push_back(mf_getHydrophobicMoment(SequenceWithFlankingRegion, 11, 120));
	Attributes.push_back(mf_getBFactorPrediction(SequenceWithFlankingRegion, 3, 5));

	VL2 = mf_getVL2Disorder(SequenceWithFlankingRegion);
	VL2iter = VL2.begin();
	for (; VL2iter != VL2.end(); VL2iter++)
		Attributes.push_back((*VL2iter));

	return Attributes;
}

ostream& operator <<(ostream & output, CPeptideAttribute peptideAttribute)
{
	vector<double>::iterator peptideAttrIter = peptideAttribute.m_vecAttributes.begin();
	for (; peptideAttrIter!= peptideAttribute.m_vecAttributes.end()-1; peptideAttrIter++)
	{

		output << *peptideAttrIter << ";";
	}
	output << *peptideAttrIter;
	return output;
}