#include"stdafx.h"
CLog flog;
bool CLog::mf_Init(string strResultPath,ios_base::openmode mode)
{
	string strLogPath = strResultPath + "\\log.txt";
	fstr.open(strLogPath.c_str(), mode);
	if (!fstr)
	{
		m_bIfReady = true;
		return false;
	}
	else
	{
		m_bIfReady = true;
		return true;
	}
}
void CLog::mf_Input(string str)
{
	if (!fstr.is_open())
	{
		cout << "Please check the log file is opening." << endl;	
		exit(-4);
	}
	else
	{
		fstr << str;
	}

}
void CLog::mf_Destroy()
{
	m_bIfReady = false;
	fstr.close();
}