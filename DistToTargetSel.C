//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class DistToTargetSel : public ModuleClass{
	public :
		DistToTargetSel()		{};
		~DistToTargetSel()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

};

bool DistToTargetSel::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Distance to target selection");
	//    -------- Name of the cut ---------     //
	Name		= "Select: dplane to tgt";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("DistToTargetSel/Do"))
	{
		Log->info( "Distance to target selection turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
		H.DefineTH1D( "DistToTargetSel", "pstart-trgt_before",	"First plane of the helices before the closest track to target cut",44,0.5,44.5);
		H.DefineTH1D( "DistToTargetSel", "pstart-trgt_after",	"First plane of the helices after the closest track to target cut",44,0.5,44.5);


	//	 --------- Parameters initialization ---------		//

	return true;
}

bool DistToTargetSel::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int d;
	// Skip the cut if there is only one track left
	if ( E.seltrack.size() > 1)
	{
		// KEEP only the closest tracks to the target
		int MinDistToTarget	= 99;
		for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
		{
			if ( E.is_upstreamdk)
				H.Fill("pstart-trgt_before",E.dcmax[*t]);
			else
				H.Fill("pstart-trgt_before",E.dcmin[*t]);
			// Distance of the track to the target
			d = DistanceToTarget( E, *t);
			if ( d < MinDistToTarget)
				MinDistToTarget = d;
		}
	
		for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
		{
			if (DistanceToTarget( E, *t) > MinDistToTarget)
			{
				E.seltrack.erase(t);	// First erase
				t--;					// then decrement to avoid to skip the following track
			}
			else
			{
				if (E.is_upstreamdk)
					H.Fill("pstart-trgt_after",E.dcmax[*t]);
				else
					H.Fill("pstart-trgt_after",E.dcmin[*t]);
			}
		}
	}

	// This should never happen but it is safer to check.
	if ( E.seltrack.size() < 1)
	{
		Log->warn("DistToTargetSel: In event %i, the number of tracks after the selection was %i",E.nevt,E.seltrack.size());
		H.CutApplied(Name);
		return false;
	}
	return true;
}

