//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 
// Special class for the PACT analysis
class HitConfig{
	public:
		HitConfig(string name, int fg);
		~HitConfig() {};
		int GetQuad( double x, double y );
		string n56;
		int	n5;
		int	n6;
		string s5;
		string s6;
		int flag;
		int quadrant;
		double slopea;
		double slopeb;
		double intercepta;
		double interceptb;
};

class PACTCut : public ModuleClass{
	public :
		PACTCut()		{};
		~PACTCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

		map<string, int> quadrant;
		map<string, int> region;

		double pc5sum;
		double pc6sum;
		vector<HitConfig> v56;
};

HitConfig::HitConfig(string name, int fg)
{
	n56 = name;
	s5	= n56[0];
	s6	= n56[1];
	n5	= atoi(s5.c_str());
	n6	= atoi(s6.c_str());
	quadrant= 1;
	slopea	= 0.0;
	slopeb	= 0.0;
	intercepta	= 0.0;
	interceptb	= 0.0;
}
int HitConfig::GetQuad( double x, double y )
{
	//returns the quadrant for (x,y)
	double linea = intercepta + x*slopea;
	double lineb = interceptb + x*slopeb;
	if (y > linea && y > lineb) return 2;
	if (y > linea && y < lineb) return 1;
	if (y < linea && y > lineb) return 3;
	if (y < linea && y < lineb) return 4;
	return 0;
}

bool PACTCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register PACT Cut");
	//    -------- Name of the cut ---------     //
	Name		= "PACT";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("PACTCut/Do"))
	{
		Log->info( "PACT Cut turned OFF");
		return false ;
	}
	
	v56.push_back(HitConfig("11",2));
	v56.push_back(HitConfig("12",3));
	v56.push_back(HitConfig("21",4));
	v56.push_back(HitConfig("22",5));

	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D("PACT","PACT_flag","Flag describing the PC hit multiplicity",7,-0.5, 6.5);

	int nb		= Conf.read<int>("PACTCut/NbBins");
	double maxw	= Conf.read<double>("PACTCut/MaximumWidth");	// ns
	string Title;

	for (int i = 0; i < v56.size(); i++)
	{
		// Title = 
		H.DefineTH2D("PACT","pc6vs5_"+v56[i].n56+"_before", "Pulse widths before cut: "+v56[i].s5+" in PC5, "+v56[i].s6+" in PC6;"+ "Total PC5 width (ns); Total PC6 width(ns)",nb,0.,maxw,nb,0.,maxw);
		H.DefineTH2D("PACT","pc6vs5_"+v56[i].n56+"_after", "Pulse widths after cut: "+v56[i].s5+" in PC5, "+v56[i].s6+" in PC6;"+ "Total PC5 width (ns); Total PC6 width(ns)",nb,0.,maxw,nb,0.,maxw);
	}
	////histos to examine wire separation
	H.DefineTH1D("PACT","pc5_hit_sep","Separation of double hits in PC5; number wires apart",20,-0.5,19.5);
	H.DefineTH1D("PACT","pc6_hit_sep","Separation of double hits in PC6; number wires apart",20,-0.5,19.5);

	////special diagnostic histograms
	if (E.Exists("npc7"))
		H.DefineTH2D("PACT","pc7vs6","Pulse widths for PC6,7; PC6 width (ns); PC7 width(ns)",nb,0.,maxw,nb,0.,maxw);
	if (E.Exists("npc7") && E.Exists("npc8"))
		H.DefineTH2D("PACT","pc8vs7","Pulse widths for PC7,8; PC7 width (ns); PC8 width(ns)",nb,0.,maxw,nb,0.,maxw);

	vector<int> s;
	//// (1 quadrant + (2*slope + 2*intercept))*4 histos = 20 parameters!
	for (int i = 0; i < v56.size(); i++)
	{
		v56[i].quadrant  = Conf.read<int>("PACTCut/quadrant_"+v56[i].n56);
		if( v56[i].quadrant < 1 or v56[i].quadrant > 4)
		{
			Log->info("Bad quadrant. PACT Cut de-registered.");
			return false;
		}
		s = StrToIntVect(Conf.read<string>("PACTCut/region_"+v56[i].n56));
		v56[i].slopea		= s[0];
		v56[i].slopeb		= s[1];
		v56[i].intercepta	= s[2];
		v56[i].interceptb	= s[3];
	}

	return true;
}

bool PACTCut::Process(EventClass &E, HistogramFactory &H)
{
	// Useful plot for understanding the PACT quadrants
	//
	// y (pc6)  b       a
	// ^         \  2  /
	// |          \   /
	// |           \ /
	// |        1   X  3
	// |           / \                 .
	// |          /   \                .
	// |         /  4  \               .
	// |
	//  ------------------------------->x (pc5)  
	//
	H.NbCandidateTracks(Name,E);


	//// _______________ fill PC7 histos, regardless of cut _________________
	if( E.Exists("npc7"))
		if( E.npc6 == 1 and E.npc7 == 1) H.Fill("pc7vs6",E.pc6width[0],E.pc7width[0]);
	if( E.Exists("npc7") and E.Exists("npc8"))
		if( E.npc7 == 1 and E.npc8 == 1) H.Fill("pc8vs7",E.pc7width[0],E.pc8width[0]);

	//// _______________ reject evts with no PC5/6 hits, or too many _________________
	if( E.npc5==0 || E.npc6==0)
	{
		H.Fill("PACT_flag",0.0);
		H.CutApplied(Name);
		return false;
	}
	if( E.npc5>2 || E.npc6>2)
	{
		H.Fill("PACT_flag",1.0);
		H.CutApplied(Name);
		return false;
	}

	//// _______________ sum multiple hits on adjacent wires _________________
	////first PC5
	pc5sum = E.pc5width[0];
	if( E.npc5 == 2)
	{
		int abs_sep = abs(E.pc5wire[1]-E.pc5wire[0]);
		H.Fill("pc5_hit_sep",abs_sep);
		if( abs_sep < 1.5)
			pc5sum += E.pc5width[1];
		else
		{
			H.CutApplied(Name);
			return false; ////double hits not adjacent -- can"t analyse
		}
	}
	////then PC6
	pc6sum = E.pc6width[0];
	if( E.npc6 == 2)
	{
		int abs_sep = abs(E.pc6wire[1]-E.pc6wire[0]);
		H.Fill("pc6_hit_sep",abs_sep);
		if( abs_sep < 1.5)
		{
			pc6sum += E.pc6width[1];
		}
		else
		{
			H.CutApplied(Name);
			return false; ////double hits not adjacent -- can"t analyse
		}
	}

	//// _______________ make the cuts _________________
	for (int i = 0; i < v56.size(); i++)
	{
		if( E.npc5==v56[i].n5 and E.npc6 ==v56[i].n6)
		{
			H.Fill("PACT_flag", v56[i].flag );
			H.Fill("pc6vs5_"+v56[i].n56+"_before",pc5sum,pc6sum);
			if( v56[i].GetQuad(pc5sum, pc6sum) == v56[i].quadrant)
			{
				H.Fill("pc6vs5_"+v56[i].n56+"_after",pc5sum,pc6sum);
				return true;
			}
			else
			{
				H.CutApplied(Name);
				return false;
			}
		}
	}
	
	//// _______________ nothing should get here _________________
	H.Fill("PACT_flag",6.0);

	H.CutApplied(Name);
	return false; ////=event fails cut
}

