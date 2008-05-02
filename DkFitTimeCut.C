//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class DkFitTimeCut : public ModuleClass{
	public :
		DkFitTimeCut()		{};
		~DkFitTimeCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

		double DkFitTimeMin;
		double DkFitTimeMax;

		double c;
};

bool DkFitTimeCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register DkFitTime Cut");
	//    -------- Name of the cut ---------     //
	Name		= "DkFitTimeCut";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("DkFitTimeCut/Do"))
	{
		Log->info( "DkFitTime Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "DkFitTimeCut", "DkFitTime_before",	"Decay Fit Time before the cut",1000, 0.,10000.);
	H.DefineTH1D( "DkFitTimeCut", "DkFitTime_after",	"Decay Fit Time after the cut", 1000, 0.,10000.);

	//	 --------- Parameters initialization ---------		//
	DkFitTimeMin	= Conf.read<double>("DkFitTimeCut/Min");	// in ns
	DkFitTimeMax	= Conf.read<double>("DkFitTimeCut/Max");	// in ns

	c				= Conf.read<double>("Parameters/c");	// cm/ns


	return true;
}

bool DkFitTimeCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int Trk			= E.seltrack[0];
	double DkFitTime	= E.hefit_time[Trk] - ((E.hefit_z[Trk] / E.costh[Trk]) / c);

	H.Fill("DkFitTime_before",DkFitTime);
	if ( not ( DkFitTime > DkFitTimeMin and DkFitTime < DkFitTimeMax ))
	{
		H.CutApplied(Name);
		return false;
	}

	H.Fill("DkFitTime_after",DkFitTime);

	return true;
}

