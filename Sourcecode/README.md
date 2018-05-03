# For developers
## Download
The source codes of LFAQ can be downloaded for free from the github.
## Pre-Requisities(Software)
* Have a <a href="https://www.visualstudio.com">Visual Studio</a> (version 2013 or above) installed.
* The <a href="http://www.microsoft.com/net/download/framework">.net framework</a> (version 4.5 or above) should be installed in advance. 
* DownLoad and install the <a href="http://www.boost.org/"> Boost C++ library</a>（version 1.57.0 or above） in advance for compiling the project.
* <a href="http://xerces.apache.org/xerces-c/">Xerces-C++</a> (version 3.1.4 or above) is needed for reading XML data. To use Xerces-C++ in your project, you should build it with your compiler at first (See <a href="http://xerces.apache.org/xerces-c/build-3.html">details</a>).   
* To use the function of GUI visulization, the <a href=" https://www.r-project.org/">R</a> (vesion 3.2.5 or above) is needed.

## Compiling the solution
1. Open the sln file in .\LFAQ\Sourcecode\ProteinAbsoluteQuan by visual studio 2013 (or above).
2. Before compiling the project, the following properties must be configured at first:
	1. Load project
	  * Include Directories:
	 	 *  Add ..\ProteinAbsoluteQuan;
	  	 *  Add the "Boost include path"; 
	 	 *  Add the "Xerces inlcude path";
	 * Additional Include Directories
	 	 * Add the "Xerces lib path";
	 * Additional Dependencies
	 	 * Add "xerces-c\_3\_1D.lib" or "xerces-c\_3\_1.lib" depending on whether you choose "release" mode or "debug" mode.
	 * Copy DLL file of Xerces to the debug and release directory. 
	2. ProteinAbsoluteQuan project
	  * Include Directories:
	  	 *  Add the "Boost include path"; 

## Contact
May you have a pleasure journey with LFAQ. If you have any problems or advices, please let us know.  
Email:  
Cheng Chang(<a href="mailto:1987ccpacer@163.com">1987ccpacer@163.com</a>) or   
  ZhiQiang Gao (<a href="mailto:gaozhiqiang13@mails.ucas.ac.cn">gaozhiqiang13@mails.ucas.ac.cn</a>)  