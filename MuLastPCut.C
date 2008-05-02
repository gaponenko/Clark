//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class MuLastPCut : public ModuleClass{
	public :
		MuLastPCut()		{};
		~MuLastPCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		int MuonLastPlane;

};

bool MuLastPCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Mu Last Plane Cut");
	//    -------- Name of the cut ---------     //
	Name		= "mu last plane";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("MuLastPCut/Do"))
	{
		Log->info( "Mu Last Plane Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "MuLastPCut", "MuLastP_before",	"Muon last plane before the cut", 58, -0.5,57.5);
	H.DefineTH1D( "MuLastPCut", "MuLastP_after",	"Muon last plane after the cut",  58, -0.5,57.5);

	//	 --------- Parameters initialization ---------		//
	MuonLastPlane	= Conf.read<int>("MuLastPCut/LastPlane");

	return true;
}

bool MuLastPCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	H.Fill("MuLastP_before",E.muon_plast);
	if ((E.muon_plast != MuonLastPlane) || (MuonLastPlane == 5678 && E.muon_plast<27 && E.muon_plast>30))
	{
		H.CutApplied(Name);
		return false;
	}
	H.Fill("MuLastP_after",E.muon_plast);

	return true;
}

