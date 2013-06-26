#include "HistogramFactory.h"

bool HistogramFactory::Init(log4cpp::Category *TmpLog, const char* Filename)
{
	Log		= TmpLog;
	Log->info("Storing the histograms in %s", Filename);
	File	= new TFile(Filename, "RECREATE");
	if ( File == NULL )
	{
		Log->error("Could not create the ROOT file %s for the histograms.", Filename);
		return false;
	}
	return true;
}

// =============================== definition section ==================================

TH1D* HistogramFactory::DefineTH1D(string Path, string Name, string Title, int xBins, double xMin, double xMax)
{
	H1D[Name]	= new TH1D(Name.c_str(), Title.c_str(), xBins, xMin, xMax);
	Store(H1D[Name], Name, Path, "TH1D");
        return H1D[Name];
}

TH1D* HistogramFactory::DefineTH1D_varwidth(string Path, string Name, string Title, vector<double> BinVect)
{
	double *bins = new double[BinVect.size()];
	for ( int i = 0; i < BinVect.size(); i++)
	{
		bins[i]	= BinVect[i];
	}	
	H1D[Name]	= new TH1D(Name.c_str(), Title.c_str(), BinVect.size()-1, bins);
	H1D[Name]->Sumw2();			// use sqrt(sum of weights) instead of sqrt(entries)
	H1D[Name]->SetOption("e1");	// show horizontal error bars
	Store(H1D[Name], Name, Path, "TH1D");
        return H1D[Name];
}

TH2D* HistogramFactory::DefineTH2D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax)
{
	H2D[Name]	= new TH2D(Name.c_str(), Title.c_str(), xBins, xMin, xMax, yBins, yMin, yMax);
	Store(H2D[Name], Name, Path, "TH2D");
        return H2D[Name];
}

TProfile2D* HistogramFactory::DefineTProfile2D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax, const std::string& options)
{
  TProfile2D *res = new TProfile2D(Name.c_str(), Title.c_str(), xBins, xMin, xMax, yBins, yMin, yMax, options.c_str());
  Store(H2D[Name] = res, Name, Path, "TH2D");
  return res;
}

TH2D* HistogramFactory::DefineTH2D_Yvarwidth(string Path, string Name, string Title, int xBins, double xMin, double xMax, vector<double> yBinVect)
{
	double *ybins = new double[yBinVect.size()];
	for ( int i = 0; i < yBinVect.size(); i++)
	{
		ybins[i]	= yBinVect[i];
	}	
	H2D[Name]	= new TH2D(Name.c_str(), Title.c_str(), xBins, xMin, xMax, yBinVect.size()-1, ybins);
	H2D[Name]->Sumw2();			// use sqrt(sum of weights) instead of sqrt(entries)
	Store(H2D[Name], Name, Path, "TH2D");
        return H2D[Name];
}

TH3D* HistogramFactory::DefineTH3D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax, int zBins, double zMin, double zMax)
{
	H3D[Name]	= new TH3D(Name.c_str(), Title.c_str(), xBins, xMin, xMax, yBins, yMin, yMax, zBins, zMin, zMax);
	Store(H3D[Name], Name, Path, "TH3D");
        return H3D[Name];
}

void HistogramFactory::DefineArrayOfStr(string Path, string Name)
{
	ObjArray[Name]	= new TObjArray();
	ObjArray[Name]->SetName(Name.c_str());
	Store(ObjArray[Name], Name, Path, "TObjArray");
}

void HistogramFactory::DefineObjString( string Path, string Name)
{
	ObjString[Name] = new TObjString();
	// ObjString[Name]->SetName(Name.c_str());
	Store(ObjString[Name], Name, Path, "TObjString");
}



void HistogramFactory::AddCut(string InCut)
{
	Cuts[InCut]	= OrderedCuts.size();
	OrderedCuts.push_back(InCut);
}

void HistogramFactory::DefineEventsBeforeCut(const char* Path, const char* Title)
{
	DefineTH1D( "",	"EventsBeforeCut",Title,Cuts.size(), -0.5,Cuts.size()-0.5);
	for (int b = 0; b < OrderedCuts.size(); b++)
		H1D["EventsBeforeCut"]->GetXaxis()->SetBinLabel( b+1, OrderedCuts[b].c_str() );
}

void HistogramFactory::DefineEventsRejectedByCut(const char* Path, const char* Title)
{
	DefineTH1D( "",	"EventsRejectedByCut",Title,Cuts.size(), -0.5,Cuts.size()-0.5);
	for (int b = 0; b < OrderedCuts.size(); b++)
		H1D["EventsRejectedByCut"]->GetXaxis()->SetBinLabel( b+1, OrderedCuts[b].c_str() );
}

void HistogramFactory::DefineNbExtraTracksPerCut(const char* Path, const char* Title, int yBins, double yMin, double yMax)
{
	DefineTH2D( "",	"NbExtraTracksPerCut",Title,Cuts.size(), -0.5,Cuts.size()-0.5, yBins, yMin, yMax);
	for (int b = 0; b < OrderedCuts.size(); b++)
		H2D["NbExtraTracksPerCut"]->GetXaxis()->SetBinLabel( b+1, OrderedCuts[b].c_str() );
}

void HistogramFactory::Store(TObject *obj, string Name, string Path, string TmpType)
{
  TH1 *hh = dynamic_cast<TH1*>(obj);
  if(hh) {
    hh->SetDirectory(0); // avoid name clashes
  }
  Dir[Name]	= Path;
  Type[Name]	= TmpType;
  File->cd();
  int opt = ((TmpType == "TObjArray") || (TmpType == "TObjString")) ? TObject::kSingleKey: 0;
  OrderedHists.push_back(StoredObjectData(obj, Path, Name, opt));
}


// ================================= Fill section ======================================

void HistogramFactory::Fill( string Name, double Val)
{
		H1D[Name]->Fill(Val);
}

void HistogramFactory::Fill( string Name, double Val1, double Val2)
{
	if ( Type[Name] == "TH2D" )
		H2D[Name]->Fill(Val1, Val2);
	else
		H1D[Name]->Fill(Val1, Val2);
}

void HistogramFactory::Fill( string Name, double Val1, double Val2, double Val3)
{
	if ( Type[Name] == "TH2D" )
		H2D[Name]->Fill(Val1, Val2, Val3);
	else
		H3D[Name]->Fill(Val1, Val2, Val3);
}

void HistogramFactory::ListToTObjArr( string GlobArr, vector<string> List)
{
	TObjArray *TmpObjA	= new TObjArray();
	for(vector<string>::iterator N = List.begin(); N != List.end(); N++)
	{
		TmpObjA->Add(new TObjString((*N).c_str()));
	}
	ObjArray[GlobArr]->Add(TmpObjA);
}

void HistogramFactory::AddStrToArray( string Name, string InStr)
{
	ObjArray[Name]->Add(new TObjString(InStr.c_str()));
}

void HistogramFactory::AddObjToArray( string Name, TObject *Obj)
{
	ObjArray[Name]->Add(Obj);
}

void HistogramFactory::FillObjString( string Name, string InStr)
{
	ObjString[Name]->SetString(InStr.c_str());
}


void	HistogramFactory::CutApplied(string LastCut)
{ 
	if( H1D.find("EventsBeforeCut") != H1D.end() ) 
		for ( int i = 0; i < Cuts[LastCut]+1; i++)
		{
			Fill("EventsBeforeCut",i);
		}
	if( H1D.find("EventsRejectedByCut") != H1D.end() ) 
	{
		Fill("EventsRejectedByCut",Cuts[LastCut]);
	}
}

void	HistogramFactory::NbCandidateTracks(string ComingCut, EventClass &E)
{
	if( H2D.find("NbExtraTracksPerCut") != H2D.end() ) 
	{
		Fill("NbExtraTracksPerCut",Cuts[ComingCut], E.seltrack.size()-1);
	}
}

// =========================== Miscellaneous section ===================================

void HistogramFactory::Sumw2( string Name)
{
	if ( Type[Name] == "TH1D" )
		H1D[Name]->Sumw2();
	else if ( Type[Name] == "TH2D" )
		H2D[Name]->Sumw2();
	else if ( Type[Name] == "TH3D" )
		H3D[Name]->Sumw2();
}


// ============================= File saving section ===================================
namespace {
  TDirectory *mkdir_p(TDirectory *parent, const std::string& path) {
    if(path.empty()) {
      return parent;
    }
    else {
      const std::string::size_type slashpos = path.find('/');
      const std::string mydir(path.substr(0, slashpos));
      const std::string subpath(((slashpos != std::string::npos) && (slashpos+1 < path.size())) ?
                                path.substr(slashpos + 1) : "" );

      TDirectory *my = parent->GetDirectory(mydir.c_str(), false);
      if(!my) {
        // The directory doesn't exist. Create it.
        my = parent->mkdir(mydir.c_str());
      }
      return mkdir_p(my, subpath);
    }
  }
}

void HistogramFactory::DoDirectory(string FilePath)
{
  File->cd();
  mkdir_p(File, FilePath);
  File->cd(FilePath.c_str());
}


void HistogramFactory::Save()
{
  int H = H1D.size() + H2D.size();
  if ( H > 0)
    {
      for(vector<StoredObjectData>::iterator N=OrderedHists.begin(); N != OrderedHists.end(); N++)
        {
          Log->info("Storing %s %s/%s with opt = %d",Type[N->name].c_str(), (N->path).c_str(), (N->name).c_str(), N->writeOpt);
          DoDirectory(N->path);
          N->obj->Write(N->name.c_str(), N->writeOpt);
          File->cd();
        }
    }
}

HistogramFactory::~HistogramFactory() {
  if(File) File->Close();
  delete File;
}
