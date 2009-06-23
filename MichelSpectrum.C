//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class MichelSpectrum : public ModuleClass{
	public :
		MichelSpectrum()		{};
		MichelSpectrum(string N, string T)		{Name = N; Title = T;};
		~MichelSpectrum()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		string Title;

		int evt;
		vector<string> EvtTypeStr;
		vector<int> EvtTypeInt;

		int accflag;
};

bool MichelSpectrum::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Log	= TmpLog;

	string N = Name;
	string T = Title;

	int ncbins;
	int nxbins;
	double xmin;
	double xmax;
	
	//	 --------- Parameters initialization ---------		//

	ncbins		= Conf.read<int>("Parameters/NCosThBinsMichel");
	nxbins		= Conf.read<int>("Parameters/NXBinsMichel");
	xmin		= Conf.read<double>("Parameters/XMinMichel");
	xmax		= Conf.read<double>("Parameters/XMaxMichel");

	EvtTypeInt	= StrToIntVect(Conf.read<string>("EventTypeCut/Types"));

	//	 --------- Histograms initialization ---------		//
	H.DefineTH2D( "Michel", "Spectrum_"+Name, "Michel spectrum "+Title+";Momentum [MeV];cos(#theta)",nxbins,xmin,xmax,ncbins, -1.0, 1.0);
	H.DefineTH2D( "Michel", "Spectrum_PzVsPt_"+Name, "Longitudinal vs transverse momentum spectrum "+Title+";P_{t} [MeV];P_{z} [MeV]",nxbins,xmin,xmax,2*nxbins, -1.*xmax, xmax);

	// IMPORTANT for mcfit !
	H.Sumw2("Spectrum_"+Name);
	H.Sumw2("Spectrum_PzVsPt_"+Name);

	string TmpEvt;
	for(vector<int>::iterator evt=EvtTypeInt.begin(); evt != EvtTypeInt.end(); evt++)
	{
		TmpEvt = IntToStr(*evt);
		EvtTypeStr.push_back(TmpEvt);
		H.DefineTH2D( "Michel_PerEvtTypes", "Spectrum_"+Name+"_EvtType"+TmpEvt, "Michel spectrum for event type "+TmpEvt+" "+Title+";Momentum [MeV];cos(#theta)",nxbins,xmin,xmax,ncbins, -1.0, 1.0);
		H.DefineTH2D( "Michel_PerEvtTypes", "Spectrum_PzVsPt_"+Name+"_EvtType"+TmpEvt, "Longitudinal vs transverse momentum spectrum for event type "+TmpEvt+" "+Title+";P_{t} [MeV];P_{z} [MeV]",nxbins,xmin,xmax,2*nxbins, -1.*xmax, xmax);
	}

	Log->info( "Register Michel spectrum ");
	return true;
}

bool MichelSpectrum::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	string N = Name;
	// The selected track
	int Trk = E.seltrack[0];

	accflag = 1; // Default value for data
	
	if (E.Exists("micheld_accflag") && E.Exists("ndecays"))
	// Extract the accflag for MC
	{
		// To store later in the output file
		E.cumulative_accflag |= E.micheld_accflag[0];
		// MC runs with the accflag for the derivatives mostly
		if (! E.ndecays > 0)
		{
			// If there is no decay from the micheld the weight is 0. (Probably still have an entry though.)
			accflag = 0;
		}
		else
		{
			try
			{
				accflag = MichelWeight(E.micheld_accflag[0]);
			}
			catch( const char* Msg)
			{
				Log->crit("MichelSpectrum: %s = 0x%x",Msg,E.micheld_accflag[0]);
				exit(1);
			}
		}
		
	}


	H.Fill("Spectrum_"+N,		E.ptot[Trk],E.costh[Trk],	accflag);
	H.Fill("Spectrum_PzVsPt_"+N,E.pt[Trk],	E.pz[Trk],		accflag);
	// for(vector<int>::iterator et = EvtTypeStr.begin(); et != EvtTypeStr.end(); et++)
	for( int i = 0; i < EvtTypeInt.size(); i++)
	{
		if ( EvtTypeInt[i] == E.type )
		{
			H.Fill("Spectrum_"+N+"_EvtType"+EvtTypeStr[i],		E.ptot[Trk],E.costh[Trk],	accflag);
			H.Fill("Spectrum_PzVsPt_"+N+"_EvtType"+EvtTypeStr[i],E.pt[Trk],	E.pz[Trk],		accflag);
			break;
		}
	}
	return true;
}

