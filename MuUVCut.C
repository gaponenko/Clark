//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class MuUVCut : public ModuleClass{
	public :
		MuUVCut()		{};
		~MuUVCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		double MuonMaxRadius;
};

bool MuUVCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Mu UV Cut");
	//    -------- Name of the cut ---------     //
	Name		= "mu UV";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("MuUVCut/Do"))
	{
		Log->info( "Mu UV Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "MuUVCut", "MuR_before",	"Muon radius before the cut;Radius [cm]", 100, 0.0,10.0);
	H.DefineTH1D( "MuUVCut", "MuR_after",		"Muon radius after the cut;Radius [cm]", 100, 0.0,10.0);

	H.DefineTH2D( "MuUVCut", "MuUV_before",	"Muon last V vs U before the cut;U position [cm];V position [cm]",60,-6.,6.,60,-6.,6.);
	H.DefineTH2D( "MuUVCut", "MuUV_after",	"Muon last V vs U after the cut;U position [cm];V position [cm]",60,-6.,6.,60,-6.,6.);

	//	 --------- Parameters initialization ---------		//
	MuonMaxRadius	= Conf.read<double>("MuUVCut/MaxRadius");


	return true;
}

bool MuUVCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	
	H.Fill("MuUV_before",E.muon_ulast,E.muon_vlast);
	H.Fill("MuR_before",E.muon_radius);
	
	if( E.muon_radius >= MuonMaxRadius)
	{
		H.CutApplied(Name);
		return false;
	}

	H.Fill("MuUV_after",E.muon_ulast,E.muon_vlast);
	H.Fill("MuR_after",E.muon_radius);


	return true;
}

