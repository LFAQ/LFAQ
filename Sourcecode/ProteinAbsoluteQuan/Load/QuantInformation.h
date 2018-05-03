#include"loadParam.h"
using namespace std;


class CQuantInformation
{
public:
	map<string, int> m_mapPeptideSequenceAndAppearNumber;
	map<string, string> m_mapExperimentNameAndPeptideIntensityName;
	map<string, vector<string>> mapProteinAndPeptideIDs;
	map<string, map<string, double>>  mapPeptideAndExperimenIntensity;
	map<string, string> mapPeptideIdAndSequence;
	string strProteinAccession;
	vector<string> vecPeptideIds;
	vector<string> vecExperimentNames;
	string StrPeptideId;
};