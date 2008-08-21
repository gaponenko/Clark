//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class TCAPm12widthCut : public ModuleClass{
	public :
		TCAPm12widthCut()		{};
		~TCAPm12widthCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string TCAPName;
		string M12Name;

		double cptime_min;
		double cptime_max;
		double m12width_min;
		double m12width_max;
};

bool TCAPm12widthCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register TCAP and m12width Cut");
	//    -------- Name of the cuts ---------     //
	TCAPName	= "TCAP";
	M12Name		= "M12 width";


	//	 --------- Special list of cuts in this class for the global histograms ---------		//
	H.AddCut(TCAPName);
	H.AddCut(M12Name);
	
	if (not Conf.read<bool>("TCAPm12widthCut/Do"))
	{
		Log->info( "TCAP and m12width Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "TCAPm12widthCut","TCAP_before",		"TCAP before the cut",200, 0,200);
	H.DefineTH1D( "TCAPm12widthCut","TCAP_after",		"TCAP after the cut",200, 0,200);
	H.DefineTH2D( "TCAPm12widthCut","m12VsTCAP_before",	"m12 vs TCAP before the cut",200, 0,200,200, 0, 20000);
	H.DefineTH2D( "TCAPm12widthCut","m12VsTCAP_after",	"m12 vs TCAP after the cut"	,200, 0,200,200, 0, 20000);

	//	 --------- Parameters initialization ---------		//
	cptime_min		= Conf.read<double>("TCAPm12widthCut/TCAPMin");
	cptime_max		= Conf.read<double>("TCAPm12widthCut/TCAPMax");
	m12width_min	= Conf.read<double>("TCAPm12widthCut/m12widthMin");
	m12width_max	= Conf.read<double>("TCAPm12widthCut/m12widthMax");


	return true;
}

bool TCAPm12widthCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ TCAP cut ____ //
	int TCAP_Good	= -1;		// -1 means that none of the three cptime values passed the cut

	H.NbCandidateTracks(TCAPName,E);
	for( int i = 0; i < 3; i++)
	{
		H.Fill("TCAP_before",E.cptime[i]);
		H.Fill("m12VsTCAP_before",E.cptime[i], E.m12width);
		if ( E.cptime[i] >= cptime_min && E.cptime[i] <= cptime_max )
		{
			TCAP_Good = i;
			H.Fill("TCAP_after",E.cptime[i]);
		}
	}

	if (TCAP_Good == -1)
	{
		H.CutApplied(TCAPName);
		return false;
	}

	// ____ m12width cut ____ //
	H.NbCandidateTracks(M12Name,E);
	// cout<<" M12 = "<<E.m12width<<"      "<<m12width_min<<"    "<<m12width_max<<endl;
	if ( E.m12width > m12width_min && E.m12width < m12width_max )
	{
		for( int i = 0; i < 3; i++)
			if ( E.cptime[i] >= cptime_min && E.cptime[i] <= cptime_max )
				H.Fill("m12VsTCAP_after",E.cptime[i], E.m12width);
	}
	else
	{
		H.CutApplied(M12Name);
		return false;
	}

	return true;
}

