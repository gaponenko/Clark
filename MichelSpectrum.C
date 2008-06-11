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

		int ncbins;
		int nxbins;
		double xmin;
		double xmax;
};

bool MichelSpectrum::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Log	= TmpLog;

	string N = Name;
	string T = Title;
	
	//	 --------- Parameters initialization ---------		//

	ncbins		= Conf.read<int>("Parameters/NCosThBinsMichel", 110);
	nxbins		= Conf.read<int>("Parameters/NXBinsMichel", 110);
	xmin		= Conf.read<double>("Parameters/XMinMichel", 0.0);
	xmax		= Conf.read<double>("Parameters/XMaxMichel", 55.0);

	//	 --------- Histograms initialization ---------		//
	string HistName = "Spectrum_"+Name;
	string HistTitl = "Michel spectrum "+Title+";Momentum [MeV];cos(theta)";
	H.DefineTH2D( "Michel", HistName, HistTitl.c_str(),nxbins,xmin,xmax,ncbins, -1.0, 1.0);

	Log->info( "Register Michel spectrum ");
	return true;
}

bool MichelSpectrum::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	string N = Name;
	// The selected track
	int Trk = E.seltrack[0];
	
	if (E.Exists("micheld_accflag") and E.Exists("ndecays"))
	{
		// MC runs with the accflag for the derivatives mostly
		if (not E.ndecays > 0)
			// If there is no decay from the micheld the weight is 0. (Probably still have an entry though.)
			H.Fill("Spectrum_"+N,E.ptot[Trk],E.costh[Trk], 0);
		try
		{
			H.Fill("Spectrum_"+N,E.ptot[Trk],E.costh[Trk], MichelWeight(E.micheld_accflag[0]));
		}
		catch( const char* Msg)
		{
			Log->crit("MichelSpectrum: %s = 0x%x",Msg,E.micheld_accflag[0]);
			exit(1);
		}
		
	}
	else
		H.Fill("Spectrum_"+N,E.ptot[Trk],E.costh[Trk]);

	return true;
}

