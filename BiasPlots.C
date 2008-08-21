//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class BiasPlots : public ModuleClass{
	public :
		BiasPlots()		{};
		BiasPlots(string N, string T)		{Name = N; Title = T;};
		~BiasPlots()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		string Title;

		int Nptotbins;
		double minptot, maxptot, dptot;

		int Trk;

		float p, costh, t, invcosth;
		float pdiff, cosdiff, wtdiff, dtdiff, thetadiff;

		double binlow, binhigh;
		string binlbl;
};

bool BiasPlots::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Log	= TmpLog;

	if ( (not Conf.read<bool>("BiasPlots/Do")) || not ( E.Exists("nmctr") && E.Exists("nmcvtx") ))
	{
		Log->info( "BiasPlots %s turned OFF",Title.c_str());
		return false ;
	}
	Log->info( "Register Bias Plots %s", Title.c_str());

	string N = Name;
	string T = Title;


	//	 --------- Parameters initialization ---------		//
	Nptotbins		= Conf.read<int>("BiasPlots/Nptotbins");
	minptot			= Conf.read<double>("BiasPlots/minptot");
	maxptot			= Conf.read<double>("BiasPlots/maxptot");

	dptot = (maxptot - minptot)/Nptotbins;

	int		Nbinscosth	= Conf.read<int>("BiasPlots/Nbinscosth");
	double	mincosth	= Conf.read<double>("BiasPlots/mincosth");
	double	maxcosth	= Conf.read<double>("BiasPlots/maxcosth");
	int		Nbinsdp		= Conf.read<int>("BiasPlots/Nbinsdp");
	double	mindp		= Conf.read<int>("BiasPlots/mindp");
	double	maxdp		= Conf.read<double>("BiasPlots/maxdp");
	int		Nbinsdcosth	= Conf.read<int>("BiasPlots/Nbinsdcosth");
	double	mindcosth	= Conf.read<double>("BiasPlots/mindcosth");
	double	maxdcosth	= Conf.read<double>("BiasPlots/maxdcosth");

	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH2D( "BiasPlots","diff_dt_cos_"+N,	"Fit - MC t vs. MC cos(theta) "+T,		100,-1.,1.,100,-20.,20.);
    H.DefineTH2D( "BiasPlots","diff_wt_cos_"+N,	"Win - MC t vs. MC cos(theta) "+T,		100,-1.,1.,100,-20.,20.);
    H.DefineTH2D( "BiasPlots","diff_wt_mct_"+N,	"Win - MC t vs. MC t "+T,				100,0.,10000.,100,-20.,20.);
    H.DefineTH2D( "BiasPlots","diff_dt_mct_"+N,	"Fit - MC t vs. MC t "+T,				100,0.,10000.,100,-20.,20.);
    H.DefineTH2D( "BiasPlots","diff_p_cos_"+N,	"Fit - MC p vs. Fit - MC cos(theta) "+T,200,-2.,2.,220,-55.,55.); 
    H.DefineTH1D( "BiasPlots","diff_cos_"+N,	"Fit - MC cos(theta) "+T,				2000,-2.,2.);
    H.DefineTH2D( "BiasPlots","diff_p_2d_"+N,	"Fit - MC p vs. p "+T,					110,0.,55.,220,-55.,55.);
    H.DefineTH2D( "BiasPlots","diff_p_2d2_"+N,	"Fit - MC p vs. cos(theta) "+T,			100,-1.,1.,220,-55.,55.);
    H.DefineTH2D( "BiasPlots","diff_cos_2d_"+N,	"Fit - MC cos(theta) vs. cos(theta) "+T,100,-1.,1.,200,-2.,2.);
    H.DefineTH2D( "BiasPlots","diff_cos_2d2_"+N,"Fit - MC cos(theta) vs. p "+T,			110,0.,55.,220,-2.,2.);
    H.DefineTH1D( "BiasPlots","ndof_all_"+N,	"Degrees of freedom "+T,				101,-0.5,100.5);
    H.DefineTH1D( "BiasPlots","chi2_all_"+N,	"chi square per degrees of freedom "+T,	100,0.,100.);

	char Temp[200];
	for(int istep=-20; istep<20; istep++)
	{
		binlow = istep/20.;
		binhigh = binlow + 0.05;
		if((-0.975 < binlow && binlow < -0.425) || (0.375 < binlow && binlow < 0.925))
		{
			// bin is in region to plot, 0.40 < abs(cos(theta)) < 0.95
			binlbl = IntToStr(int(fabs(binlow)*100.+0.5));
			if(binlow < 0.0)
				binlbl = "neg" + binlbl;

			sprintf(Temp,"ndof vs p for %f < cos(#theta) < %f %s;p [MeV];ndof",						binlow,binhigh, T.c_str());
			H.DefineTH2D( "BiasPlots", "ndof_p_"+binlbl+"_"+N,		Temp,	22,0.,55.,101,-0.5,100.5);

			sprintf(Temp,"chi2 per dof vs p for %f < cos(#theta) < %f %s;p [MeV];reduced #chi^{2}",	binlow,binhigh, T.c_str());
			H.DefineTH2D( "BiasPlots", "chi2_dof_p_"+binlbl+"_"+N,	Temp,	22,0.,55.,200,0.,20.);

			sprintf(Temp,"Fit - MC cos(#theta) vs. p for %f < cos(#theta) < %f %s;p [MeV];#Delta cos(#theta)",	binlow,binhigh, T.c_str());
			H.DefineTH2D( "BiasPlots", "dcos_p_"+binlbl+"_"+N,		Temp,	22,0.,55.,2000,-0.2,0.2);

			sprintf(Temp,"Fit - MC p vs. p for %f < cos(#theta) < %f %s;p [MeV];#Delta p [MeV]",	binlow,binhigh, T.c_str());
			H.DefineTH2D( "BiasPlots", "dp_p_"+binlbl+"_"+N,		Temp,	22,0.,55.,2200,-55.,55.);
		}
	}



	for(int i = 0; i < Nptotbins; i++)
	{
		binlow = minptot + i*dptot;
		binhigh= minptot + (i + 1)*dptot;
		binlbl = IntToStr(int(binlow));

		sprintf( Temp,	"%f < p_{tot} < %f; 1/cos#theta; #Delta p_{tot} (p^{US}_{rec} - p^{US}_{MCTB})", binlow, binhigh);
		H.DefineTH2D( "BiasPlots", "dptot_invc_us_"		+binlbl+"_"+N,		Temp,	Nbinscosth, 1./maxcosth, 1./mincosth, Nbinsdp, mindp, maxdp); 

		sprintf( Temp,	"%f < p_{tot} < %f; 1/cos#theta; #Delta cos#theta (cos#theta^{US}_{rec} - cos#theta^{US}_{MCTB})", binlow, binhigh);
		H.DefineTH2D( "BiasPlots", "dcosth_invc_us_"+binlbl+"_"+N,		Temp,	Nbinscosth, 1./maxcosth, 1./mincosth, Nbinsdcosth, mindcosth, maxdcosth); 

		sprintf( Temp,	"%f < p_{tot} < %f; 1/cos#theta; #Delta #theta (#theta^{US}_{rec} - #theta^{US}_{MCTB})", binlow, binhigh);
		H.DefineTH2D( "BiasPlots", "dtheta_invc_us_"+binlbl+"_"+N,		Temp,	Nbinscosth, 1./maxcosth, 1./mincosth, Nbinsdcosth, mindcosth, maxdcosth); 

		sprintf( Temp,	"%f < p_{tot} < %f; 1/cos#theta; #Delta p_{tot} (p^{DS}_{rec} - p^{DS}_{MCTB})", binlow, binhigh);
		H.DefineTH2D( "BiasPlots", "dptot_invc_ds_"		+binlbl+"_"+N,		Temp,	Nbinscosth, 1./maxcosth, 1./mincosth, Nbinsdp, mindp, maxdp); 

		sprintf( Temp,	"%f < p_{tot} < %f; 1/cos#theta; #Delta cos#theta (cos#theta^{DS}_{rec} - cos#theta^{DS}_{MCTB})", binlow, binhigh);
		H.DefineTH2D( "BiasPlots", "dcosth_invc_ds_"+binlbl+"_"+N,		Temp,	Nbinscosth, 1./maxcosth, 1./mincosth, Nbinsdcosth, mindcosth, maxdcosth); 

		sprintf( Temp,	"%f < p_{tot} < %f; 1/cos#theta; #Delta #theta (#theta^{DS}_{rec} - #theta^{DS}_{MCTB})", binlow, binhigh);
		H.DefineTH2D( "BiasPlots", "dtheta_invc_ds_"+binlbl+"_"+N,		Temp,	Nbinscosth, 1./maxcosth, 1./mincosth, Nbinsdcosth, mindcosth, maxdcosth); 

	}

	return true;
}

bool BiasPlots::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	// There is no vertex bank !
	if ( ! (E.Exists("nmctr") && E.Exists("nmcvtx") ) )
		return true;
	
	// If no decay positron, nothing to do.
	if ( E.tb_e_trki < 0)
		return true;

	// The selected track
	Trk		= E.seltrack[0];

	p		= E.mcvertex_ptot[E.tb_e_firstdcvtx];
	costh	= E.mcvertex_costh[E.tb_e_firstdcvtx];
	invcosth= 1./costh;
	t		= E.mcvertex_time[E.tb_e_firstdcvtx];

	pdiff		= E.ptot[Trk] - p;
	cosdiff		= E.costh[Trk] - costh;
	thetadiff	= acos(E.costh[Trk]) - acos(costh);
	wtdiff		= E.win_time[E.iewin] - t;
	dtdiff		= E.hefit_time[Trk] - t;
	
    if(20.0 < p && p < 50.0)
	{
		H.Fill("diff_cos_"		+Name,	cosdiff			);
		H.Fill("diff_cos_2d_"	+Name,	costh, cosdiff	);

		H.Fill("diff_dt_cos_"	+Name,	costh, dtdiff	);
		H.Fill("diff_wt_cos_"	+Name,	costh, wtdiff	);
		H.Fill("diff_wt_mct_"	+Name,	t, wtdiff		);
		H.Fill("diff_dt_mct_"	+Name,	t, dtdiff		);
	}

    if(0.54 < fabs(costh) && fabs(costh) < 0.80)
		H.Fill("diff_p_2d_"	+Name,	p, pdiff		);

    H.Fill("diff_p_cos_"+Name,		cosdiff, pdiff	);
    H.Fill("diff_p_2d2_"+Name,		costh,pdiff		);
    H.Fill("diff_cos_2d2_"+Name,	p, cosdiff		);

    H.Fill("ndof_all_"+Name,		E.hefit_ndof[Trk]	);
    H.Fill("chi2_all_"+Name,		E.hefit_chi2[Trk]	);

	for(int istep=-20; istep<20; istep++)
	{
		binlow = istep/20.;
		if((-0.975 < binlow && binlow < -0.425) || (0.375 < binlow && binlow < 0.925) &&
			(binlow < costh && costh < binlow + 0.05) )
		{
			// bin is in region to plot, 0.40 < abs(cos(theta)) < 0.95
			binlbl = IntToStr(int(fabs(binlow)*100.+0.5));
			if(binlow < 0.0)
				binlbl = "neg" + binlbl;

			H.Fill( "ndof_p_"		+binlbl+"_"+Name, E.ptot[Trk], E.hefit_ndof[Trk]);
			H.Fill( "chi2_dof_p_"	+binlbl+"_"+Name, E.ptot[Trk], E.hefit_chi2[Trk]/E.hefit_ndof[Trk] );
			H.Fill( "dcos_p_"		+binlbl+"_"+Name, p, cosdiff);
			H.Fill( "dp_p_"			+binlbl+"_"+Name, p, pdiff);
		}
	}


	for(int i = 0; i < Nptotbins; i++)
	{
		binlow = minptot + i*dptot;
		binhigh= minptot + (i + 1)*dptot;
		binlbl = IntToStr(int(binlow));

		if( p > binlow && p <= binhigh)
		{
			if ( E.is_upstreamdk )
			{
				H.Fill( "dptot_invc_us_"	+binlbl+"_"+Name,	invcosth,	pdiff);
				H.Fill( "dcosth_invc_us_"	+binlbl+"_"+Name,	invcosth,	cosdiff);
				H.Fill( "dtheta_invc_us_"	+binlbl+"_"+Name,	invcosth,	thetadiff);
			}
			else
			{
				H.Fill( "dptot_invc_ds_"	+binlbl+"_"+Name,	invcosth,	pdiff);
				H.Fill( "dcosth_invc_ds_"	+binlbl+"_"+Name,	invcosth,	cosdiff);
				H.Fill( "dtheta_invc_ds_"	+binlbl+"_"+Name,	invcosth,	thetadiff);
			}
		}
	}

	return true;
}

