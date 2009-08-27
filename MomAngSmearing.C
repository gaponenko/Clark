//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class MomAngSmearing : public ModuleClass{
	public :
		MomAngSmearing()		{};
		~MomAngSmearing()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		int Trk;
		double tmp_costh;
		double tmp_ptot;

};

bool MomAngSmearing::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register MomAngSmearing");

	
	if (not Conf.read<bool>("MomAngSmearing/Do"))
	{
		Log->info( "MomAngSmearing turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH2D( "MomAngSmearing", "MomVsCosthDifference",	"Difference between after and before the smearing;Momentum difference[MeV];costheta difference", 201, -0.1, 0.1, 101, -0.1, 0.1);
	H.DefineTH2D( "MomAngSmearing", "MomVsThetaDifference",	"Difference between after and before the smearing;Momentum difference[MeV];Theta difference [rad]", 201, -0.1, 0.1, 101, -0.1, 0.1);

	//	 --------- Parameters initialization ---------		//

	return true;
}

bool MomAngSmearing::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	Trk = E.seltrack[0];

	tmp_ptot	= sqrt((E.hefit_pu[Trk]*E.hefit_pu[Trk])+(E.hefit_pv[Trk]*E.hefit_pv[Trk])+(E.hefit_pz[Trk]*E.hefit_pz[Trk]));
	if (not tmp_ptot == 0)
		tmp_costh	= E.hefit_pz[Trk] / tmp_ptot;
	else
		tmp_costh	= 1.0;

	H.Fill("MomVsCosthDifference", E.ptot[Trk] - tmp_ptot, E.costh[Trk] - tmp_costh );
	H.Fill("MomVsThetaDifference", E.ptot[Trk] - tmp_ptot, acos(E.costh[Trk]) - acos(tmp_costh) );


	return true;
}

