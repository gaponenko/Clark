//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class Mu_eVertexCut : public ModuleClass{
	public :
		Mu_eVertexCut()		{};
		~Mu_eVertexCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

		bool CalcVsMofia;
		int CutFunction;
		vector<float> CutParameters;

		bool MustErase;
		vector<double> u0;
		vector<double> v0;
		vector<double> u0_Calc;
		vector<double> v0_Calc;
		double du, dv, r;
};

bool Mu_eVertexCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Mu-e vertex cut");
	//    -------- Name of the cut ---------     //
	Name		= "mu-e vertex";

	Log	= TmpLog;

	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("Mu_eVertexCut/Do", 1))
	{
		Log->info( "Mu-e vertex cut turned OFF");
		return false ;
	}
	
	//	 --------- Parameters initialization ---------		//
	CalcVsMofia		= Conf.read<bool>("Mu_eVertexCut/Do_mu_e_Calc-MOFIA");
	CutFunction		= Conf.read<int>("Mu_eVertexCut/CutFunction");
	CutParameters	= StrToFloatVect(Conf.read<string>("Mu_eVertexCut/CutParameters"));

	switch(CutFunction)
	{
		case 0:
			if( CutParameters.size() != 2 )
			{
				Log->error("Mu_eVertexCut: Wrong number of parameters for the function 0. A minimum and a maximum only are required.");
				exit(1);
			}
	}
	//	 --------- Histograms initialization ---------		//
	H.DefineTH2D( "Mu_eVertexCut", "Cut_mu_e_dv_vs_du_up_before",	"Mu-e dV vs dU before mu-e vertex cut, upstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH2D( "Mu_eVertexCut", "Cut_mu_e_dv_vs_du_down_before",	"Mu-e dV vs dU before mu-e vertex cut, downstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH1D( "Mu_eVertexCut", "Cut_mu_e_r_up_before",			"Mu-e ellipse distance, before mu-e vertex cut, upstream tracks", 200, 0.0, 2.0);
	H.DefineTH1D( "Mu_eVertexCut", "Cut_mu_e_r_down_before",			"Mu-e ellipse distance, before mu-e vertex cut, downstream tracks", 200, 0.0, 2.0);
	H.DefineTH2D( "Mu_eVertexCut", "Cut_mu_e_dv_vs_du_up_after",		"Mu-e dV vs dU after mu-e vertex cut, upstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH2D( "Mu_eVertexCut", "Cut_mu_e_dv_vs_du_down_after",	"Mu-e dV vs dU after mu-e vertex cut, downstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH1D( "Mu_eVertexCut", "Cut_mu_e_r_up_after",			"Mu-e ellipse distance, after mu-e vertex cut, upstream tracks", 200, 0.0, 2.0);
	H.DefineTH1D( "Mu_eVertexCut", "Cut_mu_e_r_down_after",			"Mu-e ellipse distance, after mu-e vertex cut, downstream tracks", 200, 0.0, 2.0);
	H.DefineTH2D( "Mu_eVertexCut", "Cut_Track_v_vs_u_up_before",		"Track V vs U  at the target before mu-e vertex cut, upstream tracks", 200, -50.0, 50.0, 200, -50.0, 50.0);
	H.DefineTH2D( "Mu_eVertexCut", "Cut_Track_v_vs_u_down_before",	"Track V vs U  at the target before mu-e vertex cut, downstream tracks", 200, -50.0, 50.0, 200, -50.0, 50.0);


	return true;
}

bool Mu_eVertexCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	
	u0.clear();
	v0.clear();
	u0_Calc.clear();
	v0_Calc.clear();
	double tmpu, tmpv;

	// Little convenient output in the log file
	if (E.seltrack.size() > 1)
		Log->info("Event %d",E.nevt);

	if ( ( not (E.Exists("hefit_u0") and E.Exists("hefit_v0")) ) or CalcVsMofia)
	{
		for ( int t = 0; t < E.ntr; t++)
		{
			// Hardcoded ierror == 0.
			// Otherwise the function Get_uv_at will have divisions by 0
			if (E.hefit_ierror[t] == 0)
			{
				Get_uv_at( &E, t, 0.0, tmpu, tmpv);		// Position target hardcoded at 0.0
				u0_Calc.push_back(tmpu);
				v0_Calc.push_back(tmpv);
			}
			// If bad track, use crazy big numbers
			else
			{
				u0_Calc.push_back(1000);
				v0_Calc.push_back(1000);
			}
		}
	}

	if (E.Exists("hefit_u0") and E.Exists("hefit_v0"))
	// ===========> Position at the target from MOFIA <=========== //
		for ( int t = 0; t < E.ntr; t++)
		{
			u0.push_back(E.hefit_u0[t]);
			v0.push_back(E.hefit_v0[t]);
		}
	else
	// ===========> Position at the target calculated <=========== //
	{
		u0	= u0_Calc;
		v0	= v0_Calc;
	}

	// ===========> Loop over the tracks <============= //

	for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
	{
		MustErase = false;
		du	= E.muon_ulast - u0[*t];
		dv	= E.muon_vlast - v0[*t];
		r	= sqrt(du * du + dv * dv);
		if (E.is_upstreamdk)
		{
			H.Fill("Cut_Track_v_vs_u_up_before",u0[*t],v0[*t]);
			H.Fill("Cut_mu_e_dv_vs_du_up_before",du,dv);
			H.Fill("Cut_mu_e_r_up_before",r);
		}
		else
		{
			H.Fill("Cut_Track_v_vs_u_down_before",u0[*t],v0[*t]);
			H.Fill("Cut_mu_e_dv_vs_du_down_before",du,dv);
			H.Fill("Cut_mu_e_r_down_before",r);
		}

		switch(CutFunction)
		{
			case 0:
				if( ! (CutParameters[0] < r && r < CutParameters[1]))
					MustErase = true;
			default:
				return true;
		}
		if (MustErase)
		{
			E.seltrack.erase(t);	// First erase
			t--;					// then decrement to avoid to skip the following track
		}
		else
		{
			if (E.is_upstreamdk)
			{
				H.Fill("Cut_mu_e_dv_vs_du_up_after",du,dv);
				H.Fill("Cut_mu_e_r_up_after",r);
			}
			else
			{
				H.Fill("Cut_mu_e_dv_vs_du_down_after",du,dv);
				H.Fill("Cut_mu_e_r_down_after",r);
			}
		}
	}


	if (E.seltrack.size() < 1)
	{
		H.CutApplied(Name);
		return false;
	}
	return true;
}
