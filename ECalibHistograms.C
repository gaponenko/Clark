//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class ECalibHistograms : public ModuleClass{
	public :
		ECalibHistograms()		{};
		~ECalibHistograms()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		void FillHistos(EventClass &E, HistogramFactory &H, int weight);

		log4cpp::Category *Log;

		int T;		// Index of the selected track
};

bool ECalibHistograms::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register ECalib Histograms");
	
	// if (not Conf.read<bool>("Template/Do", 1))
	// {
	// 	Log->info( " turned OFF");
	// 	return false ;
	// }
	//	 --------- Parameters initialization ---------		//
	int ncbins			= Conf.read<int>("Parameters/NCBins_EC");
	int ncbins_invc		= Conf.read<int>("Parameters/NCBins_ECinvc");
	int nxbins			= Conf.read<int>("Parameters/NXBins_EC");
	double xmin			= Conf.read<double>("Parameters/XMin_EC");
	double xmax			= Conf.read<double>("Parameters/XMax_EC");
	double cosmin		= Conf.read<double>("Parameters/CosMin_EC");
	double cosmax		= Conf.read<double>("Parameters/CosMax_EC");

	// c_bins			= [0.0] * (ncbins + 1)
	// invc_bins		= [0.0] * (ncbins_invc + 1)
	vector<double> c_bins;
	vector<double> invc_bins;
	double d_theta			= M_PI/ncbins;
	double invcosdiff		= 1./cosmin - 1./cosmax;
	double d_sec			= 2.0 * invcosdiff / ncbins_invc;

	for ( int i = 0; i < ncbins+1; i++)
		c_bins.push_back(-cos(i * d_theta));

	for ( int j = 0; j < ncbins_invc+1; j++)
	{
		if (j < ncbins_invc/2)
			invc_bins.push_back(1. / (-1./cosmax - j * d_sec ));
		else
			invc_bins.push_back(1. / ( 1./cosmin - ( j - ncbins_invc/2 ) * d_sec ));
	}

	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH2D_Yvarwidth( "ECalib", "constth",			"Michel spectrum for ECal, binned in #theta",		nxbins, xmin,xmax, c_bins);
	H.DefineTH2D( "ECalib", "constcosth",		"Michel spectrum for ECal, binned in cos#theta",	nxbins, xmin,xmax,100,-1.,1.);
	H.DefineTH2D_Yvarwidth( "ECalib", "constinvcosth",	"Michel spectrum for ECal, binned in 1/cos#theta",	nxbins, xmin,xmax, invc_bins);


	return true;
}

bool ECalibHistograms::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	T =	E.seltrack[0];

	if (E.Exists("micheld_accflag") && E.Exists("ndecays"))
	{
		// MC runs with the accflag.
		if ( ! E.ndecays > 0)
		{
			// If there is no decay from the micheld the weight is 0.
			// (Probably still have an entry though.)
			FillHistos( E, H, 0);
			return true;
		}
		try
		{
			FillHistos( E, H, MichelWeight(E.micheld_accflag[0]));	// weight from accflag
		}
		catch( const char* Msg)
		{
			Log->crit("ECalibHistograms: %s = 0x%x",Msg,E.micheld_accflag[0]);
			exit(1);
		}
	}
	else
		// Data runs or MC runs without accflag.
		FillHistos( E, H, 1);		// Simply weight = 1

	return true;
}

void ECalibHistograms::FillHistos(EventClass &E, HistogramFactory &H, int weight)
{
	H.Fill("constth",E.ptot[T], E.costh[T], weight);
	H.Fill("constcosth",E.ptot[T], E.costh[T], weight);
	H.Fill("constinvcosth",E.ptot[T], E.costh[T], weight);
}
