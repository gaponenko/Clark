//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class TemplateCut : public ModuleClass{
	public :
		TemplateCut()		{};
		~TemplateCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

};

bool TemplateCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Name		= "";

	Log	= TmpLog;

	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("Template/Do", 1))
	{
		Log->info( " turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//

	//	 --------- Parameters initialization ---------		//
	// cptime_min		= Conf.read<double>("TemplateCut/TCAPMin", 0.0);
	// cptime_max		= Conf.read<double>("TemplateCut/TCAPMax", 200.0);
	// m12width_min	= Conf.read<double>("TemplateCut/m12widthMin", 0.0);
	// m12width_max	= Conf.read<double>("TemplateCut/m12widthMax", 20000.0);


	Log->info( "Register");
	return true;
}

bool TemplateCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	if (Bad)
	{
		H.CutApplied(Name);
		return false;
	}
	return true;
}

