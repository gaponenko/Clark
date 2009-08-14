//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class ChiSquare : public ModuleClass{
	public :
		ChiSquare()		{};
		~ChiSquare()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog);
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		double RedChi2Min;
		double RedChi2Max;
		double RedChi2;

		int Trk;
};

bool ChiSquare::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Log	= TmpLog;

	if ( ! Conf.read<bool>("ChiSquare/Do") ) //Default OFF
	{
		Log->info( "ChiSquare Module turned OFF");
		return false ;
	}

	//	 --------- Parameters initialization ---------		//
	int ncbins		= Conf.read<int>("ChiSquare/NCosThBinsMichel");
	int nxbins		= Conf.read<int>("ChiSquare/NXBinsMichel");
	double xmin		= Conf.read<double>("ChiSquare/XMinMichel");
	double xmax		= Conf.read<double>("ChiSquare/XMaxMichel");

	RedChi2Min		= Conf.read<double>("ChiSquare/RedChi2CutMin");
	RedChi2Max		= Conf.read<double>("ChiSquare/RedChi2CutMax");

	//	 --------- Histograms initialization ---------		//
	H.DefineTH3D( "ChiSquare", "Chi2vsCosTvsP",		"Reduced #chi^{2} distribution versus cos#theta versus Momentum;P [MeV];cos#theta;#chi^{2}",nxbins,xmin,xmax,ncbins, -1.0, 1.0, 1000, 0,10);
	H.DefineTH3D( "ChiSquare", "Chi2vsCosTvsPFull",	"Reduced #chi^{2} distribution versus cos#theta versus Momentum;P [MeV];cos#theta;#chi^{2}",nxbins,xmin,xmax,ncbins, -1.0, 1.0, 1000, 0,1000);

	Log->info( "Register ChiSquare Module");
	return true;
}

bool ChiSquare::Process(EventClass &E, HistogramFactory &H)
{
	Trk		= E.seltrack[0];
	RedChi2	= E.hefit_chi2[Trk] / E.hefit_ndof[Trk];

	if (E.hefit_ndof[Trk] > 0)
	{
			H.Fill("Chi2vsCosTvsP",		E.ptot[Trk],E.costh[Trk], RedChi2);
			H.Fill("Chi2vsCosTvsPFull",	E.ptot[Trk],E.costh[Trk], RedChi2);
	}

	if ( RedChi2 < RedChi2Min )
		return false;

	// Cut max turned off if CutMax <= 0.0
	if ( RedChi2 > RedChi2Max && RedChi2Max > 0.0 )
		return false;

	return true;
}

