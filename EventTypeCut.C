//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class EventTypeCut : public ModuleClass{
	public :
		EventTypeCut()		{};
		~EventTypeCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		vector<int>	TypeList;

};

bool EventTypeCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Event Type Cut");
	//    -------- Name of the cut ---------     //
	Name		= "Evt Type";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("EventTypeCut/Do"))
	{
		Log->info( "Event Type Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "EventTypeCut","EvtType_before",	"Event Type before the cut;Event Type",31, -0.5,31.5);
	H.DefineTH1D( "EventTypeCut","EvtType_after",	"Event Type after the cut;Event Type",31, -0.5,31.5);

	//	 --------- Parameters initialization ---------		//
	TypeList	= StrToIntVect(Conf.read<string>("EventTypeCut/Types"));

	return true;
}

bool EventTypeCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ Event Type Cut ____ //
	H.NbCandidateTracks(Name,E);
	H.Fill("EvtType_before",E.type);
	bool Type_OK	= false;
	for(vector<int>::iterator TR=TypeList.begin(); TR != TypeList.end(); TR++)
	{
		if( E.type == *TR)
		{
			Type_OK = true;
			break;
		}
		if( E.type == (-1*(*TR)))
		{
			Type_OK = false;
			break;
		}
	}

	if( not Type_OK)
	{
		H.CutApplied(Name);
		return false;
	}
		
	H.Fill("EvtType_after",E.type);

	return true;
}

