//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class DkWinTimeCut : public ModuleClass{
	public :
		DkWinTimeCut()		{};
		~DkWinTimeCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		double DkWinTimeMin;
		double DkWinTimeMax;
};

bool DkWinTimeCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register DkWinTime Cut");
	//    -------- Name of the cut ---------     //
	Name		= "DkWinTime";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("DkWinTimeCut/Do"))
	{
		Log->info( "DkWinTime Cut turned OFF");
		return false;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "DkWinTimeCut", "DkWinTime_before",	"Decay window time before the cut",4*160, -6000.,10000.);
	H.DefineTH1D( "DkWinTimeCut", "DkWinTime_after",	"Decay window time after the cut", 4*160, -6000.,10000.);

	//	 --------- Parameters initialization ---------		//
	DkWinTimeMin	= Conf.read<double>("DkWinTimeCut/Min");
	DkWinTimeMax	= Conf.read<double>("DkWinTimeCut/Max");


	return true;
}

bool DkWinTimeCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	H.Fill("DkWinTime_before",E.win_time[E.iewin]);

	if( not ( E.win_time[E.iewin] > DkWinTimeMin and E.win_time[E.iewin] < DkWinTimeMax ))
	{
		H.CutApplied(Name);
		return false;
	}

	H.Fill("DkWinTime_after",E.win_time[E.iewin]);

	return true;
}

