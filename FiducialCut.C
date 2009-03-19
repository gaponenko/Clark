//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class FiducialCut : public ModuleClass{
	public :
		FiducialCut()		{};
		~FiducialCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

		double MinAbscosth;
		double MaxAbscosth;
		double MaxMomentum;
		double MinTransMom;
		double MaxTransMom;
		double MinLongiMom;
};

bool FiducialCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Fiducial Cut");
	//    -------- Name of the cut ---------     //
	Name		= "Fiducial";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	//	 --------- Histograms initialization ---------		//

	//	 --------- Parameters initialization ---------		//
	MinAbscosth		= Conf.read<double>("FiducialCut/MinAbsCosTheta");
	MaxAbscosth		= Conf.read<double>("FiducialCut/MaxAbsCosTheta");
	MaxMomentum		= Conf.read<double>("FiducialCut/MaxMomentum");
	MinTransMom		= Conf.read<double>("FiducialCut/MinTransMom");
	MaxTransMom		= Conf.read<double>("FiducialCut/MaxTransMom");
	MinLongiMom		= Conf.read<double>("FiducialCut/MinLongiMom");

	return true;
}

bool FiducialCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int Trk	= E.seltrack[0];	// The selected track
	// cos theta cut
	if ( not (fabs(E.costh[Trk]) > MinAbscosth and fabs(E.costh[Trk]) < MaxAbscosth ))
	{
		H.CutApplied(Name);
		return false;
	}
	// Momentum cut
	if (not (E.ptot[Trk] < MaxMomentum))
	{
		H.CutApplied(Name);
		return false;
	}
	// Transverse momentum cut
	if (not ( E.pt[Trk] > MinTransMom and E.pt[Trk] < MaxTransMom ))
	{
		H.CutApplied(Name);
		return false;
	}
	// Longitudinal momentum cut
	if (not (fabs(E.pz[Trk]) > MinLongiMom))
	{
		H.CutApplied(Name);
		return false;
	}

	return true;
}

