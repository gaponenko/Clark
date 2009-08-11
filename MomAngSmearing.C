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
	//    -------- Name of the cut ---------     //
	Name		= "";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("MomAngSmearing/Do"))
	{
		Log->info( "MomAngSmearing turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	// H.DefineTH1D( "Dir", "val_before",	"Value before the cut", 100, -0.5,99.5);
	// H.DefineTH1D( "Dir", "val_after",	"Value after the cut", 100, -0.5,99.5);
	H.DefineTH2D( "MomAngSmearing", "MomVsCosthDifference",	"Difference between after and before the smearing;Momentum difference[MeV];costheta difference", 201, -0.1, 0.1, 101, -0.1, 0.1);
	H.DefineTH2D( "MomAngSmearing", "MomVsThetaDifference",	"Difference between after and before the smearing;Momentum difference[MeV];Theta difference [rad]", 201, -0.1, 0.1, 101, -0.1, 0.1);

	//	 --------- Parameters initialization ---------		//
	// val_min		= Conf.read<double>("MomAngSmearing/val_Min");
	// val_max		= Conf.read<double>("MomAngSmearing/val_Max");

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


//	H.NbCandidateTracks(Name,E);
//	H.Fill("val_before", E.val)
//	H.Fill("val_vs_par", E.par, E.val)
//
//	if (! ( val_min < E.val && E.val < val_max) )
//	{
//		H.CutApplied(Name);
//		return false;
//	}
//	H.Fill("val_after", E.val)
//
//	// More complicated cut on tracks (for example) in the event.
//	// Example:
//	// The vector of integer track contains the track indices in the event.
//	// The following will loop over the track indices, check a value in the
//	// array valarray for each integer *t and remove the bad tracks.
//	// 
//	// for(vector<int>::iterator t = E.track.begin(); t != E.track.end(); t++)
//	// {
//	// 	// ===> TRACK CUT HERE
//	// 	if( E.valarray[*t]   != GoodVal)
//	// 	{
//	// 		E.seltrack.erase(t);    // First erase
//	// 		t--;                    // then decrement to avoid to skip the following track
//	// 		continue;
//	// 	}
//	// }

//
//
	return true;
}

