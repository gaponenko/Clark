//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class ChargeCut : public ModuleClass{
	public :
		ChargeCut()		{};
		~ChargeCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		int SelCharge;
};

bool ChargeCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Charge Cut");
	//    -------- Name of the cut ---------     //
	Name		= "charge";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("ChargeCut/Do"))
	{
		Log->info( "Charge Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "ChargeCut", "Charge_Up_before",	"Charge distribution, upstream decays, before the cut",3, -1.5,1.5);
	H.DefineTH1D( "ChargeCut", "Charge_Up_after",		"Charge distribution, upstream decays, after the cut", 3, -1.5,1.5);
	H.DefineTH1D( "ChargeCut", "Charge_Down_before",	"Charge distribution, downstream decays, before the cut",3, -1.5,1.5);
	H.DefineTH1D( "ChargeCut", "Charge_Down_after",	"Charge distribution, downstream decays, after the cut", 3, -1.5,1.5);

	//	 --------- Parameters initialization ---------		//
	SelCharge		= Conf.read<int>("ChargeCut/SelectCharge");

	return true;
}

bool ChargeCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int Nb_Charge_OK	= 0;

	for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
	{
		if( E.is_upstreamdk)
			H.Fill("Charge_Up_before",E.hefit_q[*t]);
		else
			H.Fill("Charge_Down_before",E.hefit_q[*t]);
		// ===> TRACK CUT HERE
		if( E.hefit_q[*t]	!= SelCharge)
		{
			E.seltrack.erase(t);	// First erase
			t--;					// then decrement to avoid to skip the following track
			continue;
		}
		else
		{
			Nb_Charge_OK += 1;
			if( E.is_upstreamdk)
				H.Fill("Charge_Up_after",E.hefit_q[*t]);
			else
				H.Fill("Charge_Down_after",E.hefit_q[*t]);
		}
	}

		
	// ===> EVENT CUT HERE
	if( not Nb_Charge_OK)
	{
		H.CutApplied(Name);
		return false;
	}
	return true;
}

