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


void HistogramFactory::DefineTH1D(string Path, string Name, string Title, int xBins, double xMin, double xMax)
{
	H1D[Name]	= new TH1D(Name.c_str(), Title.c_str(), xBins, xMin, xMax);
	Store( Name, Path, "TH1D");
}

void HistogramFactory::DefineTH1D_varwidth(string Path, string Name, string Title, vector<double> BinVect)
{
	double *bins = new double[BinVect.size()];
	for ( int i = 0; i < BinVect.size(); i++)
	{
		bins[i]	= BinVect[i];
	}	
	H1D[Name]	= new TH1D(Name.c_str(), Title.c_str(), BinVect.size()-1, bins);
	H1D[Name]->Sumw2();			// use sqrt(sum of weights) instead of sqrt(entries)
	H1D[Name]->SetOption("e1");	// show horizontal error bars
	Store( Name, Path, "TH1D");
}

void HistogramFactory::DefineTH2D(string Path, string Name, string Title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax)
{
	H2D[Name]	= new TH2D(Name.c_str(), Title.c_str(), xBins, xMin, xMax, yBins, yMin, yMax);
	Store( Name, Path, "TH2D");
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

void HistogramFactory::Store( string Name, string Path, string TmpType)
{
	Dir[Name]	= Path;
	Type[Name]	= TmpType;
	File->cd();
	OrderedHists.push_back(Name);
}


void HistogramFactory::Fill( string Name, double Val)
{
	if ( Type[Name] == "TH1D" )
		H1D[Name]->Fill(Val);
	
	if ( Type[Name] == "TH2D" )
		Log->error(" Ouch ! You try to put one value in %s but it is a TH2D.",Name.c_str());
}

void HistogramFactory::Fill( string Name, double Val1, double Val2)
{
	if ( Type[Name] == "TH1D" )
		H1D[Name]->Fill(Val1, Val2);
	
	if ( Type[Name] == "TH2D" )
		H2D[Name]->Fill(Val1, Val2);
}

void HistogramFactory::Fill( string Name, double Val1, double Val2, double Val3)
{
	if ( Type[Name] == "TH1D" )
		Log->error(" Ouch ! You try to put three values in %s but it is a TH1D.",Name.c_str());
	
	if ( Type[Name] == "TH2D" )
		H2D[Name]->Fill(Val1, Val2, Val3);
}

void	HistogramFactory::CutApplied(string LastCut)
{ 
	map<string,TH1D*>::iterator Before = H1D.find("EventsBeforeCut");
	map<string,TH1D*>::iterator Reject = H1D.find("EventsRejectedByCut");
	if( Before != H1D.end() ) 
		for ( int i = 0; i < Cuts[LastCut]+1; i++)
		{
			Fill("EventsBeforeCut",i);
		}
	if( Reject != H1D.end() ) 
	{
		Fill("EventsRejectedByCut",Cuts[LastCut]);
	}
}

void	HistogramFactory::NbCandidateTracks(string ComingCut, EventClass &E)
{
	map<string,TH2D*>::iterator NbTracks = H2D.find("NbExtraTracksPerCut");
	if( NbTracks != H2D.end() ) 
	{
		Fill("NbExtraTracksPerCut",Cuts[ComingCut], E.seltrack.size()-1);
	}
}

void HistogramFactory::DoDirectory( string FilePath)
{
	File->cd();
	string Current	= "";
	if (FilePath == "")
		return;
	vector<const char*> Dirs;

	string tmp = FilePath;
	if (tmp.find("/", 0) == string::npos)
		Dirs.push_back(tmp.c_str());
	else
		while (tmp.find("/", 0) != string::npos)
		{
			Dirs.push_back(tmp.substr( 0, tmp.find("/", 0)).c_str());
			tmp	= tmp.substr( tmp.find("/", 0)+1);
		}

	TKey *K;
	for(vector<const char*>::iterator D=Dirs.begin(); D != Dirs.end(); D++)
	{
		K	= (TKey*)File->GetKey(*D);
		if (K == NULL)
		// The directory doesn't exist. Create it.
		{
			File->mkdir(*D);
			File->cd(*D);
		}
		else
		{
			if ( not K->IsFolder())
			{
			// That's embarrasing. Probably a Histogram has the name of the directory we want. Oups !
				Log->warn("The directory %s cannot be created because a histogram has already this name. The current directory \"%s\" will be used instead.",*D,Current.c_str());
				return;
			}
			else
				File->cd(*D);

		}
		Current += (*D);
		Current += "/";
	}
}


void HistogramFactory::Save()
{
	int H = H1D.size() + H2D.size();
	if ( H > 0)
	{
		TObject *Dummy = new TObject();
		for(vector<string>::iterator N=OrderedHists.begin(); N != OrderedHists.end(); N++)
		{
			Log->info("Storing %s %s",Type[*N].c_str(), (*N).c_str());
			DoDirectory(Dir[*N]);
			if (Type[*N] == "TH1D")
				H1D[*N]->Write();
			else if (Type[*N] == "TH2D")
				H2D[*N]->Write();
			// else if (Type[*N] == "TObjArray")
			// 	H[N].Write(N, Dummy.kSingleKey)
			// else:
			// 	self.H[N].Write()
			File->cd();
		}
	}
}
