#include<fstream>
#include<string>
#include<iostream>
using namespace std;

class CLog
{
public:
	CLog()
	{
		m_bIfReady = false;
	}
	bool mf_Init(string strResultPath,ios_base::openmode=ios::out);
	void mf_Destroy();
	void mf_Input(string str);
	bool m_bIfReady;
private:
	ofstream fstr; 
};
extern CLog flog;
