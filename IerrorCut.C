//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class IerrorCut : public ModuleClass{
	public :
		IerrorCut()		{};
		~IerrorCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		vector<int>	ierrorList;

};

bool IerrorCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Ierror Cut");
	//    -------- Name of the cut ---------     //
	Name		= "ierror";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("IerrorCut/Do"))
	{
		Log->info( "Ierror Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "IerrorCut", "ierror_before",	"ierror distribution before the cut",30, -0.5,29.5);
	H.DefineTH1D( "IerrorCut", "ierror_after",	"ierror distribution after the cut", 30, -0.5,29.5);
	H.DefineTH1D( "IerrorCut", "Nb_Good_ierror","Number of good ierror per event", 10, -0.5,9.5);

	//	 --------- Parameters initialization ---------		//
	ierrorList	= StrToIntVect(Conf.read<string>("IerrorCut/SelectIerror"));

	return true;
}

bool IerrorCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int Nb_ierror_OK	= 0;
	bool Track_ierror_OK;

	for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
	{
		Track_ierror_OK	= false;
		H.Fill("ierror_before",E.hefit_ierror[*t]);
		for(vector<int>::iterator i = ierrorList.begin(); i != ierrorList.end(); i++)
		{
			if( E.hefit_ierror[*t] == *i)
			{
				Track_ierror_OK	= true;
				Nb_ierror_OK += 1;
				H.Fill("ierror_after",E.hefit_ierror[*t]);
				break;
			}
		}
		// ===> TRACK CUT HERE
		if( not Track_ierror_OK)
		{
			E.seltrack.erase(t);	// First erase
			t--;					// then decrement to avoid to skip the following track
		}

	}
	// ===> EVENT CUT HERE
	if( not Nb_ierror_OK)
	{
		H.CutApplied(Name);
		return false;
	}

	return true;
}

