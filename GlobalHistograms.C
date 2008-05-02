//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class GlobalHistograms : public ModuleClass{
	public :
		GlobalHistograms()		{};
		~GlobalHistograms()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

};

bool GlobalHistograms::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Global Histograms");
	//    -------- Name of the cut ---------     //

	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut("Done");
	
	// if (not Conf.read<bool>("Template/Do", 1))
	// {
	// 	Log->info( " turned OFF");
	// 	return false ;
	// }
	
	//	 --------- Histograms initialization ---------		//
	H.DefineEventsBeforeCut("", "Number of events before each cut");
	H.DefineEventsRejectedByCut("", "Number of events rejected by the cut");
	H.DefineNbExtraTracksPerCut("", "Number of extra decay candidates before cut vs cut", 10, -0.5,9.5);

	//	 --------- Parameters initialization ---------		//

	return true;
}

bool GlobalHistograms::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks("Done",E);
	H.CutApplied("Done");
	return true;
}

