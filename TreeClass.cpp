#include "TreeClass.h"

bool TreeClass::InitAll( string Treename, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log		= TmpLog;

	// Initialization of the analyses
	Evt.Init(&Conf, Log);

	if ( not InitTree(Treename) )
	{
		Log->error("Something is wrong with the first tree of the list");
		return false;
	}

	string HistFile;
	// Initialization of the analyses
	if (Conf.read<string>("CommandLine/OutputFile", "")  == "")
		HistFile	= Conf.read<string>("Name")+".root";
	else
		HistFile	= Conf.read<string>("CommandLine/OutputFile");
	Hist	= HistogramFactory();
	Hist.Init(Log, HistFile.c_str());

	Log->info("Registering analysis classes");

	try
	{
		for(vector<ModuleClass*>::iterator i=ModList.begin(); i != ModList.end(); i++)
		{
			if (not (*i)->Init(Evt, Hist, Conf, Log))
			{
				ModList.erase(i);
				i--;
			}
		}
	}
	catch( const ConfigFile::key_not_found &e )
	{
		Log->error("The configuration key \"%s\" could not be found. Check the default configuration in DefaultConfig.cpp", e.key.c_str());
		return false;
	}

	FirstEntry	= Conf.read<int>("FirstEntry", 0);
	LastEntry	= Conf.read<int>("LastEntry", -1);

	CloseTree();

	return true;
}

bool TreeClass::InitTree( string Treename)
{
	string Filename, RunLog;

 	vector<string> Parsed;
	string Run;
	// Is this a directory or a file ?
	struct stat   fileStat;
	stat (Treename.c_str(), &fileStat);

 	if (S_ISDIR (fileStat.st_mode))
 	{
 		Parsed	= StrToStrVect(Treename,'/');
		Run;
 		if (Parsed[Parsed.size()-1] == "" && Parsed.size() > 1)
 			Run = Parsed[Parsed.size()-2];
 		else
 			Run = Parsed[Parsed.size()-1];
 		Filename	= Treename + "/tree0" + Run.substr(3) + ".root";
 		RunLog		= Treename + "/" + Run.substr(3) + "log.txt";
 	}
 	else
 	{
		Filename	= Treename;
		RunLog		= "";
	}

	TreeFile	= new TFile(Filename.c_str());
	if (TreeFile->IsZombie())
	{
		Log->info("The root file is not opened correctly");
		return false;
	}
	// Now take care of the mofia log file.
	// Check only if we decided to.

	Tree		= (TTree*) TreeFile->Get("T");

	// ======> The first entry is loaded to check which variables are in the tree <======
	// Get the next tree in the chain and verify
	int ientry,NbCheck;
	ientry	= Tree->LoadTree(0);
	if (ientry < 0)
	{
		Log->error( "Problem with the first tree ... this is not good.");
		return false;
	}

	NbCheck	= Tree->GetEntry(0);
	if (NbCheck <= 0 )
	{
		Log->error( "Problem with the entry of the first tree ... this is not good.");
		return false;
	}

	Evt.InitVar(Tree);

	return true;
}

void TreeClass::LoopTree()
{
	int ientry,NbCheck;

	Int_t NbEntries	 = Tree->GetEntries();

	// How many entries to analyse ?
	if (LastEntry == -1 or LastEntry > NbEntries)
		LastEntry	= NbEntries;
	// On what entry to start ?
	if (FirstEntry > LastEntry)
		FirstEntry	= LastEntry;

	for (int I = FirstEntry; I < LastEntry; I++)
	{
		ientry	= Tree->LoadTree(I);
		if (ientry < 0)
		{
			break;
		}

		NbCheck	= Tree->GetEntry(I);
		if (NbCheck <= 0 )
			continue;

		// TODO: Check that this is what I should do for this kind of event
		if ( not Evt.Load())
			continue;

		for(int i=0; i < ModList.size(); i++)
		{
			if ( not ModList[i]->Process(Evt, Hist ))
			{
				break;
			}
		}

	}
}

