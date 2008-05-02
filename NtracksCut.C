//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class NtracksCut : public ModuleClass{
	public :
		NtracksCut()		{};
		~NtracksCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

};

bool NtracksCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register N tracks Cut");
	//    -------- Name of the cut ---------     //
	Name		= "N tracks";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("NtracksCut/Do"))
	{
		Log->info( "N tracks Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "NtracksCut", "ntrks_before",	"Number of tracks before the cut",10, -0.5,9.5);
	H.DefineTH1D( "NtracksCut", "ntrks_after",	"Number of tracks after the cut", 10, -0.5,9.5);

	//	 --------- Parameters initialization ---------		//

	return true;
}

bool NtracksCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ N tracks Cut ____ //
	H.NbCandidateTracks(Name,E);
	H.Fill("ntrks_before",E.seltrack.size());
	if (E.seltrack.size() <= 0)
	{
		H.CutApplied(Name);
		return false;
	}
	H.Fill("ntrks_after",E.seltrack.size());
	return true;
}

