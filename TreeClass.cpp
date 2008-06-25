#include "TreeClass.h"

bool TreeClass::InitAll( string Treename, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log		= TmpLog;

	// Initialization of the EventClass
	Evt.Init(Conf, Log);

	CheckLog		= true;
	FirstInitTree	= true;
	// Load the first tree to check the available variables in the tree.
	if ( not InitTree(Treename, Conf) )
	{
		Log->error("Something is wrong with the first tree of the list");
		return false;
	}
	FirstInitTree	= false;

	// Read the geometry file only
	if ( not Evt.InitGeometry(Conf) )
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

bool TreeClass::InitTree( string Treename, ConfigFile &Conf)
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
	if (CheckLog && (RunLog != ""))
	{
		if ( ! (ReadMofiaLog(RunLog, Conf)))
		{
			if (FirstInitTree)
			{
				Log->info("The first log file didn't contain the MC information. The rest of the log files will be ignored.");
				CheckLog	= false;
			}
			else
			{
				Log->warn("The log file %s was not read properly. This could be a problem for the nthrown !", RunLog.c_str());
				return false;
			}
		}
	}

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


// --------------------------------------------------------------------------------- \\


bool TreeClass::ReadMofiaLog( string Filename, ConfigFile &Conf)
{
	std::ifstream File(Filename.c_str());
	if (FirstInitTree)
		Log->info("Reading log file %s",Filename.c_str());

	int Nbnthrown = 0;	// Number of nthrown found in this file
 	vector<string> Tmpreqnum;
	if(!File)
	{
		Log->warn("Failed to open log file\n");
		return false;
	}

	string Tmp;
	string Tmpnthrown;
	string Line;
	double BField;
	string GeoFileNum;
	while ( getline( File, Line ) )
	{
 		// Search for the nthrown
		Tmp = ReadParam( " unpMCE3\\(\\): nthrown =\\s+([0-9.]+)", Line);

 		if ( Tmp != "")
		{
 			// There is a serious problem if more than one nthrown is found
 			if (Nbnthrown > 0)
			{
 				Log->warn("Multiple nthrown values found. Not good !");
				// return false;
			}
			Nbnthrown += 1;
			Tmpnthrown	= Tmp;
		}

		// Search for the reqnum
		Tmp = ReadParam( " unpMCB.\\(\\): reqnum =\\s+(\\d+)", Line);
 		if ( Tmp != "")
		{
			// Store all the reqnum found
 			Tmpreqnum.push_back(Tmp);
		}

		// Search for the BField value
		Tmp = ReadParam( "magnet_mod: Field at origin is\\s+([0-9.]+)", Line);
 		if ( Tmp != "")
		{
			// Store the BField found
 			BField = atof(Tmp.c_str());
		}

		// Search for the BField value
		Tmp = ReadParam( "Reading geometry file:.*dt_geo.000(\\d+)", Line);
 		if ( Tmp != "")
		{
			// Store the BField found
 			GeoFileNum = Tmp.c_str();
		}
	}

	File.close();

	// This tree and runlog will be checked again.
	if (FirstInitTree)
	{
		// BField
		if ( Conf.read<double>("Detector/BField") == 2.0)
		{
			Log->warn("The BField was not set in the configuration file. Use the value from the log file = %2.6f",BField);
			Conf.add("Detector/BField", BField);
		}
		// Geometry file
		if ( Conf.read<double>("Detector/GeometryFile") < 0.0 )
		{
			Log->warn("No geometry file number given in the configuration file. Use the CFM number from the log file = %s",GeoFileNum.c_str());
			Conf.add("Detector/GeometryFile", GeoFileNum);
		}

	}
 
	// No nthrown means data run.
	// Return false means at the FirstInitTree that the log files won't be checked again.
	if (Nbnthrown == 0)
	{
		Log->info("No nthrown found.");
		return false;
	}

	// This is the check on the first tree. Don't save everything.
	if ( FirstInitTree)
		return true;

	nthrown.push_back(Tmpnthrown);
	reqnum.insert( reqnum.end(), Tmpreqnum.begin(), Tmpreqnum.end() );

	MofiaLogs.push_back(Filename);

	return true;
}

void TreeClass::StoreExtraValues()
{
	if ( ! CheckLog)
		return;
	Hist.DefineArrayOfStr("", "logdata");

	Hist.ListToTObjArr( "logdata", nthrown);
	Hist.ListToTObjArr( "logdata", reqnum);
	Hist.ListToTObjArr( "logdata", MofiaLogs);
}
