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
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH3D.h"
#include "TAxis.h"

// Include Program
#include "log4cpp/Category.hh"
#include "EventClass.h"
#include "RooUnfoldDummy.h"

class HistogramFactory {
  struct StoredObjectData {
    TObject *obj;
    std::string path;
    std::string name;
    int writeOpt;
    StoredObjectData(TObject *o, const std::string& p, const std::string& n, int opt)
      : obj(o), path(p), name(n), writeOpt(opt)
    {}
  };

public :
  HistogramFactory() : File(0) {}
  ~HistogramFactory();

		bool	Init(log4cpp::Category *TmpLog, const char* Filename);

		// =============================== definition section ==================================
		TH1D*	DefineTH1D(string Path, string Name, string Title, int xBins, double xMin, double xMax);
		TProfile* DefineTProfile(string Path, string Name, string Title, int xBins, double xMin, double xMax, const std::string& options="");
		TH2D*	DefineTH2D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax);
		TProfile2D* DefineTProfile2D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax, const std::string& options="");

		TH1D*	DefineTH1D_varwidth(string Path, string Name, string Title, vector<double> BinVect);
		TH2D*	DefineTH2D_Yvarwidth(string Path, string Name, string Title, int xBins, double xMin, double xMax, vector<double> yBinVect);
		TH3D*	DefineTH3D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax, int zBins, double zMin, double zMax);
		void	DefineArrayOfStr(string Path, string Name);
		void	DefineObjString(string Path, string Name);
		void	AddCut(string InCut);
		void	DefineEventsBeforeCut(const char* Path, const char* Title);
		void	DefineEventsRejectedByCut(const char* Path, const char* Title);
		void	DefineNbExtraTracksPerCut(const char* Path, const char* Title, int yBins, double yMin, double yMax);
		// ================================= Fill section ====================================== 
		void	Fill(string Name, double Val);
		void	Fill(string Name, double Val1, double Val2);
		void	Fill(string Name, double Val1, double Val2, double Val3);
		void	ListToTObjArr( string GlobArr, vector<string> List);
		void	AddStrToArray( string Name, string InStr);
		void	AddObjToArray( string Name, TObject *Obj);
		void	FillObjString(string Name, string InStr);
		void	CutApplied(string LastCut);
		void	NbCandidateTracks(string ComingCut, EventClass &E);

		// =========================== Miscellaneous section ===================================
		void	Sumw2(string Name);

		// ============================= File saving section ===================================
		void	Save();
		void	DoDirectory(string FilePath);
		void	Store(TObject* obj, string Name, string Path, string TmpType);
		void	Store(TObject* obj, string Name, string Path, int writeOpt = 0);
		void	Store(RooUnfoldDummy* obj, string Name, string Path, int writeOpt = 0) {}


	private :
		map<string, TH1D*> H1D;
		map<string, TH2D*> H2D;
		map<string, TH3D*> H3D;
		map<string, TObjArray*> ObjArray;
		map<string, TObjString*> ObjString;
		map<string, string> Type;
		map<string, string> Dir;
		map<string, int> Cuts;
		vector<string> OrderedCuts;
		vector<StoredObjectData> OrderedHists;
		TFile *File;
		log4cpp::Category *Log;
};


#endif
