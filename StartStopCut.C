//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class StartStopCut : public ModuleClass{
	public :
		StartStopCut()		{};
		~StartStopCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

};

bool StartStopCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register StartStop Cut");
	//    -------- Name of the cut ---------     //
	Name		= "startstop";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("StartStopCut/Do"))
	{
		Log->info( "StartStop Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH2D( "StartStopCut", "StartStop_Up_before",	"DC max vs min, upstream decays, before the cut;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);
	H.DefineTH2D( "StartStopCut", "StartStop_Up_after",		"DC max vs min, upstream decays, after the cut;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);
	H.DefineTH2D( "StartStopCut", "StartStop_Down_before",	"DC max vs min, downstream decays, before the cut;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);
	H.DefineTH2D( "StartStopCut", "StartStop_Down_after",	"DC max vs min, downstream decays, after the cut;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);

	H.DefineTH2D( "StartStopCut", "StartStop_Up_ePlus",		"DC max vs min, upstream decays, positive tracks;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);
	H.DefineTH2D( "StartStopCut", "StartStop_Up_eMinus",	"DC max vs min, upstream decays, negative tracks;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);
	H.DefineTH2D( "StartStopCut", "StartStop_Down_ePlus",	"DC max vs min, downstream decays, positive tracks;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);
	H.DefineTH2D( "StartStopCut", "StartStop_Down_eMinus",	"DC max vs min, downstream decays, negative tracks;DC min;DC max",45,-0.5,44.5,45,-0.5,44.5);


	return true;
}

bool StartStopCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	int Nb_StartStop_OK	= 0;

	for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
	{
		if( E.is_upstreamdk)
		{
			H.Fill("StartStop_Up_before",E.dcmin[*t],E.dcmax[*t]);
			if ( E.hefit_q[*t] < 0 )
				H.Fill("StartStop_Up_eMinus",E.dcmin[*t],E.dcmax[*t]);
			else
				H.Fill("StartStop_Up_ePlus",E.dcmin[*t],E.dcmax[*t]);
		}
		else
		{
			H.Fill("StartStop_Down_before",E.dcmin[*t],E.dcmax[*t]);
			if ( E.hefit_q[*t] < 0 )
				H.Fill("StartStop_Down_eMinus",E.dcmin[*t],E.dcmax[*t]);
			else
				H.Fill("StartStop_Down_ePlus",E.dcmin[*t],E.dcmax[*t]);
		}
		// ===> TRACK CUT HERE
		// Looks complicated but works. Both ends of the track must be on the side
		// we expect the decay to be according to the event type.
		if( (E.is_upstreamdk != (E.dcmin[*t] <= 22)) || (E.is_upstreamdk != (E.dcmax[*t] <= 22)) )
		{
			E.seltrack.erase(t);	// First erase
			t--;					// then decrement to avoid to skip the following track
			continue;
		}
		else
		{
			// Sanity check
			// cout <<"Event = "<<E.nevt<<"   upstream = "<<E.is_upstreamdk<<"   costh = "<<E.costh[*t]<<endl;
			if (E.is_upstreamdk != (E.costh[*t] < 0))
			{
				Log->warn("StartStopCut: cos(theta) sign does not match track position: costh=%f, dcmin=%d, dcmax=%d, wintype=%d, event=%d", E.costh[*t], E.dcmin[*t], E.dcmax[*t], E.win_type[E.iewin], E.nevt);
				E.seltrack.erase(t);	// First erase
				t--;					// then decrement to avoid to skip the following track
				continue;
			}

			Nb_StartStop_OK += 1;
			if( E.is_upstreamdk)
				H.Fill("StartStop_Up_after",E.dcmin[*t],E.dcmax[*t]);
			else
				H.Fill("StartStop_Down_after",E.dcmin[*t],E.dcmax[*t]);
		}
	}

	// ===> EVENT CUT HERE
	if( not Nb_StartStop_OK)
	{
		H.CutApplied(Name);
		return false;
	}

	return true;
}

