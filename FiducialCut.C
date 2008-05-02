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

};

bool FiducialCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Fiducial Cut");
	//    -------- Name of the cut ---------     //
	Name		= "Fiducial";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	// if (not Conf.read<bool>("FiducialCut/Do", 1))
	// {
	// 	Log->info( "Fiducial Cut turned OFF");
	// 	return false ;
	// }
	
	//	 --------- Histograms initialization ---------		//

	//	 --------- Parameters initialization ---------		//
	// TODO: Put the fiducial parameters in the conf

	return true;
}

bool FiducialCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int Trk	= E.seltrack[0];	// The selected track
	// cos theta cut
	if ( not (fabs(E.costh[Trk]) > 0.5 and fabs(E.costh[Trk]) < 0.92 ))
	{
		H.CutApplied(Name);
		return false;
	}
	// Momentum cut
	if (not (E.ptot[Trk] < 51.5))
	{
		H.CutApplied(Name);
		return false;
	}
	// Transverse momentum cut
	if (not ( E.pt[Trk] > 10.0 and E.pt[Trk] < 39.7 ))
	{
		H.CutApplied(Name);
		return false;
	}
	// Longitudinal momentum cut
	if (not (fabs(E.pz[Trk]) > 13.7))
	{
		H.CutApplied(Name);
		return false;
	}

	return true;
}

