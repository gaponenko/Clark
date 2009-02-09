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
	Log	= TmpLog;
	Log->info( "Register TemplateCut");
	//    -------- Name of the cut ---------     //
	Name		= "";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("TemplateCut/Do", 1))
	{
		Log->info( "TemplateCut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	// H.DefineTH1D( "Dir", "val_before",	"Value before the cut", 100, -0.5,99.5);
	// H.DefineTH1D( "Dir", "val_after",	"Value after the cut", 100, -0.5,99.5);
	// H.DefineTH2D( "Dir", "val_vs_par",	"Value versus a parameter", 20, 1., 21., 100, -0.5,99.5);

	//	 --------- Parameters initialization ---------		//
	// val_min		= Conf.read<double>("TemplateCut/val_Min", 0.0);
	// val_max		= Conf.read<double>("TemplateCut/val_Max", 200.0);

	return true;
}

bool TemplateCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
//	H.NbCandidateTracks(Name,E);
//	H.Fill("val_before", E.val)
//	H.Fill("val_vs_par", E.par, E.val)
//
//	if (! ( val_min < E.val && E.val < val_max) )
//	{
//		H.CutApplied(Name);
//		return false;
//	}
//	H.Fill("val_after", E.val)
	return true;
}

