// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include"CLog.h"
#include <tchar.h>
#include<time.h>

using namespace std;

enum DataType{ MaxQuantTpye = 0,mzQuantMLType=1, PeakViewType=3};
const int BUFFERLENGTH = 100000;

// TODO: reference additional headers your program requires here