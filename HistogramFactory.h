//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HistogramFactory_h
#define HistogramFactory_h
// Include de C++
using namespace std;

#include <iostream>
#include <string>
#include <map>
#include <vector>

// Include de ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

// Include Program
#include "log4cpp/Category.hh"
#include "EventClass.h"


class HistogramFactory {
	public :
		HistogramFactory()	{};
		~HistogramFactory()	{};
		void	DefineTH1D(string Path, string Name, string Title, int xBins, double xMin, double xMax);
		void	DefineTH2D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax);
		void	DefineTH1D_varwidth(string Path, string Name, string Title, vector<double> BinVect);
		void	AddCut(string InCut);
		void	DefineEventsBeforeCut(const char* Path, const char* Title);
		void	DefineEventsRejectedByCut(const char* Path, const char* Title);
		void	DefineNbExtraTracksPerCut(const char* Path, const char* Title, int yBins, double yMin, double yMax);
		void	DoDirectory(string FilePath);
		void	Fill(string Name, double Val);
		void	Fill(string Name, double Val1, double Val2);
		void	Fill(string Name, double Val1, double Val2, double Val3);
		bool	Init(log4cpp::Category *TmpLog, const char* Filename);
		void	Save();
		void	Store(string Name, string Path, string TmpType);
		void	CutApplied(string LastCut);
		void	NbCandidateTracks(string ComingCut, EventClass &E);

	private :
		map<string, TH1D*> H1D;
		map<string, TH2D*> H2D;
		map<string, string> Type;
		map<string, string> Dir;
		map<string, int> Cuts;
		vector<string> OrderedCuts;
		vector<string> OrderedHists;
		TFile *File;
		log4cpp::Category *Log;
};


#endif
