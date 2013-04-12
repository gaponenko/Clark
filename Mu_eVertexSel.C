//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class Mu_eVertexSel : public ModuleClass{
	public :
		Mu_eVertexSel()		{};
		~Mu_eVertexSel()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

		bool CalcVsMofia;
		vector<double> u0;
		vector<double> v0;
		vector<double> u0_Calc;
		vector<double> v0_Calc;
};

bool Mu_eVertexSel::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Mu-e vertex selection");
	//    -------- Name of the cut ---------     //
	Name		= "Selec: mu-e vertex";

	Log	= TmpLog;

	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("Mu_eVertexSel/Do"))
	{
		Log->info( "Mu-e vertex selection turned OFF");
		return false ;
	}

	if ( ( not (E.Exists("hefit_u0") and E.Exists("hefit_v0")) ) or CalcVsMofia)
		Log->warn("Mu_eVertexSel: Calculating the Mu-e vertex in the treesum. The code is not maintained and could be wrong !");
	
	//	 --------- Parameters initialization ---------		//
	CalcVsMofia	= Conf.read<bool>("Mu_eVertexSel/Do_mu_e_Calc-MOFIA");

	//	 --------- Histograms initialization ---------		//
	if (CalcVsMofia)
	{
		H.DefineTH2D( "Mu_eVertexSel", "Calc-MOFIA_dv_vs_du_up_before",		"Difference between calculation and MOFIA, dV vs dU before mu-e vertex select, upstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
		H.DefineTH2D( "Mu_eVertexSel", "Calc-MOFIA_dv_vs_du_down_before",	"Difference between calculation and MOFIA, dV vs dU before mu-e vertex select, downstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
		H.DefineTH1D( "Mu_eVertexSel", "Calc-MOFIA_r_up_before",				"Difference between calculation and MOFIA, ellipse distance, before mu-e vertex select, upstream tracks", 200, 0.0, 2.0);
		H.DefineTH1D( "Mu_eVertexSel", "Calc-MOFIA_r_down_before",			"Difference between calculation and MOFIA, ellipse distance, before mu-e vertex select, downstream tracks", 200, 0.0, 2.0);
		H.DefineTH1D( "Mu_eVertexSel", "Calc-MOFIA_up_mu_e_agree",			"Agreement between calculation and MOFIA, mu-e vertex select, upstream tracks;Agreement (agree=1)", 2, -0.5, 1.5);
		H.DefineTH1D( "Mu_eVertexSel", "Calc-MOFIA_down_mu_e_agree",			"Agreement between calculation and MOFIA, mu-e vertex select, downstream tracks;Agreement (agree=1)", 2, -0.5, 1.5);
	}

	H.DefineTH2D( "Mu_eVertexSel", "mu_e_dv_vs_du_up_before",	"Mu-e dV vs dU before mu-e vertex select, upstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH2D( "Mu_eVertexSel", "mu_e_dv_vs_du_down_before",	"Mu-e dV vs dU before mu-e vertex select, downstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH1D( "Mu_eVertexSel", "mu_e_r_up_before",			"Mu-e ellipse distance, before mu-e vertex select, upstream tracks", 200, 0.0, 2.0);
	H.DefineTH1D( "Mu_eVertexSel", "mu_e_r_down_before",			"Mu-e ellipse distance, before mu-e vertex select, downstream tracks", 200, 0.0, 2.0);
	H.DefineTH2D( "Mu_eVertexSel", "mu_e_dv_vs_du_up_after",		"Mu-e dV vs dU after mu-e vertex select, upstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH2D( "Mu_eVertexSel", "mu_e_dv_vs_du_down_after",	"Mu-e dV vs dU after mu-e vertex select, downstream tracks", 200, -1.0, 1.0, 200, -1.0, 1.0);
	H.DefineTH1D( "Mu_eVertexSel", "mu_e_r_up_after",			"Mu-e ellipse distance, after mu-e vertex select, upstream tracks", 200, 0.0, 2.0);
	H.DefineTH1D( "Mu_eVertexSel", "mu_e_r_down_after",			"Mu-e ellipse distance, after mu-e vertex select, downstream tracks", 200, 0.0, 2.0);
	H.DefineTH2D( "Mu_eVertexSel", "Track_v_vs_u_up_before",		"Track V vs U  at the target before mu-e vertex select, upstream tracks", 200, -50.0, 50.0, 200, -50.0, 50.0);
	H.DefineTH2D( "Mu_eVertexSel", "Track_v_vs_u_down_before",	"Track V vs U  at the target before mu-e vertex select, downstream tracks", 200, -50.0, 50.0, 200, -50.0, 50.0);



	return true;
}

bool Mu_eVertexSel::Process(EventClass &E, HistogramFactory &H)
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


	// Find the track closest to the muon at the target
	int best			= -1;
	double mindist 		= 9999999.9;
	int best_Calc		= -1;
	double mindist_Calc	= 9999999.9;
	double du, dv, r;
	double du_Calc, dv_Calc, r_Calc;
	double bestdu(0.), bestdv(0.), bestr(0.);
	for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
	{
		du	= E.muon_ulast - u0[*t];
		dv	= E.muon_vlast - v0[*t];
		r	= sqrt(du * du + dv * dv);
		if (E.is_upstreamdk)
		{
			H.Fill("Track_v_vs_u_up_before",u0[*t],v0[*t]);
			H.Fill("mu_e_dv_vs_du_up_before",du,dv);
			H.Fill("mu_e_r_up_before",r);
		}
		else
		{
			H.Fill("Track_v_vs_u_down_before",u0[*t],v0[*t]);
			H.Fill("mu_e_dv_vs_du_down_before",du,dv);
			H.Fill("mu_e_r_down_before",r);
		}
		if (r < mindist)
		{
			mindist = r;
			best	= *t;
			bestdu	= du;
			bestdv	= dv;
			bestr	= r;
		}
		// Option to check the difference between calculated track position at the target
		// and the mofia calculation stored in the tree
		if (CalcVsMofia)
		{
			du_Calc	= E.muon_ulast - u0_Calc[*t];
			dv_Calc	= E.muon_vlast - v0_Calc[*t];
			r_Calc	= sqrt(du_Calc * du_Calc + dv_Calc * dv_Calc);
			if (E.is_upstreamdk)
			{
				H.Fill("Calc-MOFIA_dv_vs_du_up_before",du_Calc-du,dv_Calc-dv);
				H.Fill("Calc-MOFIA_r_up_before",r_Calc-r);
			}
			else
			{
				H.Fill("Calc-MOFIA_dv_vs_du_down_before",du_Calc-du,dv_Calc-dv);
				H.Fill("Calc-MOFIA_r_down_before",r_Calc-r);
			}
			if (r_Calc < mindist_Calc)
			{
				mindist_Calc = r_Calc;
				best_Calc	= *t;
			}
		}
	}

	if (CalcVsMofia)
	{
		if (E.is_upstreamdk)
			H.Fill("Calc-MOFIA_up_mu_e_agree",int(best == best_Calc));
		else
			H.Fill("Calc-MOFIA_down_mu_e_agree",int(best == best_Calc));
	}

	// Best way to make sure that there is only one track after this.
	// First erase the list and add only one element
	E.seltrack.clear();
	E.seltrack.push_back(best);
	if (E.is_upstreamdk)
	{
		H.Fill("mu_e_dv_vs_du_up_after",bestdu,bestdv);
		H.Fill("mu_e_r_up_after",bestr);
	}
	else
	{
		H.Fill("mu_e_dv_vs_du_down_after",bestdu,bestdv);
		H.Fill("mu_e_r_down_after",bestr);
	}

	// These two conditions should never happen but it is safer to check.
	if (E.seltrack.size() < 1)
	{
		Log->warn("Mu_eVertexSel: In event %i, the number of tracks after the selection was %i",E.nevt,E.seltrack.size());
		H.CutApplied(Name);
		return false;
	}
	if (E.seltrack.size() > 1)
	{
		Log->error("Mu_eVertexSel: There is more than one track after the mu-e vertex selection. Abort event");
		return false;
	}
	return true;
}

