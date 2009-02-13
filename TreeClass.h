//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef TreeClass_h
#define TreeClass_h
// Include de C++
using namespace std;

#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sstream>
#include <iomanip>

// Include de ROOT
// #include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

// Include Program
#include "ModuleClass.h"
#include "EventClass.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "FuncLib.h"
// Include log4cpp
#include "log4cpp/Category.hh"

// TODO
// Change InFileName from a vector to a list

class TreeClass {
	public :
		TreeClass()	{};
		~TreeClass()	{};
		void SaveHistos()	{Hist.Save();};

		bool InitAll( string Treename, ConfigFile &Conf, log4cpp::Category *TmpLog);
		bool InitTree( string Treename, ConfigFile &Conf);
		void CloseTree() { TreeFile->Close();};
		void LoopTree();
		void Register( ModuleClass *Class);
		bool ReadMofiaLog( string Filename, ConfigFile &Conf);
		void StoreExtraValues();
	private :
		TTree*	Tree;
		TFile*	TreeFile;
		vector<ModuleClass*> ModList;
		EventClass Evt;
		log4cpp::Category *Log;
		HistogramFactory Hist;

		int FirstEntry;
		int LastEntry;

		bool FirstInitTree;
		bool CheckLog;
		vector<string> nthrown;
		vector<string> reqnum;
		vector<string> MofiaLogs;
};

inline void TreeClass::Register( ModuleClass *Class)
{
	ModList.push_back(Class);
}


#endif
