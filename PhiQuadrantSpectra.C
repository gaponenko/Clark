//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class PhiQuadrantSpectra : public ModuleClass{
	public :
		PhiQuadrantSpectra()		{};
		~PhiQuadrantSpectra()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		int evt;
		vector<string> EvtTypeStr;
		vector<int> EvtTypeInt;

		int accflag;
		double u0;
		double v0;
		double x0;
		double y0;
		vector<string> Quad;

		int Trk;
		int q;
};

bool PhiQuadrantSpectra::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Log	= TmpLog;
	if (not Conf.read<bool>("PhiQuadrantSpectra/Do"))
	{
		Log->info( "PhiQuadrantSpectra turned OFF");
		return false ;
	}
	

	int ncbins_m;
	int nxbins_m;
	double xmin_m;
	double xmax_m;
	
	//	 --------- Parameters initialization ---------		//

	ncbins_m	= Conf.read<int>("Parameters/NCosThBinsMichel");
	nxbins_m	= Conf.read<int>("Parameters/NXBinsMichel");
	xmin_m		= Conf.read<double>("Parameters/XMinMichel");
	xmax_m		= Conf.read<double>("Parameters/XMaxMichel");

	//	 --------- Parameters initialization ---------		//
	int ncbins			= Conf.read<int>("Parameters/NCBins_EC");
	int ncbins_invc		= Conf.read<int>("Parameters/NCBins_ECinvc");
	int nxbins			= Conf.read<int>("Parameters/NXBins_EC");
	double xmin			= Conf.read<double>("Parameters/XMin_EC");
	double xmax			= Conf.read<double>("Parameters/XMax_EC");
	double cosmin		= Conf.read<double>("Parameters/CosMin_EC");
	double cosmax		= Conf.read<double>("Parameters/CosMax_EC");

	// c_bins			= [0.0] * (ncbins + 1)
	// invc_bins		= [0.0] * (ncbins_invc + 1)
	vector<double> c_bins;
	vector<double> invc_bins;
	double d_theta			= M_PI/ncbins;
	double invcosdiff		= 1./cosmin - 1./cosmax;
	double d_sec			= 2.0 * invcosdiff / ncbins_invc;

	for ( int i = 0; i < ncbins+1; i++)
		c_bins.push_back(-cos(i * d_theta));

	for ( int j = 0; j < ncbins_invc+1; j++)
	{
		if (j < ncbins_invc/2)
			invc_bins.push_back(1. / (-1./cosmax - j * d_sec ));
		else
			invc_bins.push_back(1. / ( 1./cosmin - ( j - ncbins_invc/2 ) * d_sec ));
	}

	
	//	 --------- Histograms initialization ---------		//
	Quad.push_back("1");
	Quad.push_back("2");
	Quad.push_back("3");
	Quad.push_back("4");


	for ( int i = 0; i < 4; i++)
	{

		H.DefineTH2D( "PhiQuadrantSpectra", "Spectrum_UVquad"+Quad[i], 		"Michel spectrum UV quadrant "+Quad[i]+";Momentum [MeV];cos(#theta)",				nxbins_m,xmin_m,xmax_m,ncbins_m, -1.0, 1.0);
		H.DefineTH2D( "PhiQuadrantSpectra", "Spectrum_PzVsPt_UVquad"+Quad[i],	"Longitudinal vs transverse momentum spectrum UV quadrant "+Quad[i]+";P_{t} [MeV];P_{z} [MeV]",	nxbins_m,xmin_m,xmax_m,2*nxbins_m, -1.*xmax_m, xmax_m);
		H.Sumw2("Spectrum_UVquad"+Quad[i]);
		H.Sumw2("Spectrum_PzVsPt_UVquad"+Quad[i]);
		H.DefineTH2D_Yvarwidth( "PhiQuadrantSpectra", "constth_UVquad"+Quad[i],		"Michel spectrum UV quadrant "+Quad[i]+" for ECal, binned in #theta",		nxbins, xmin,xmax, c_bins);
		H.DefineTH2D_Yvarwidth( "PhiQuadrantSpectra", "constinvcosth_UVquad"+Quad[i],	"Michel spectrum UV quadrant "+Quad[i]+" for ECal, binned in 1/cos#theta",	nxbins, xmin,xmax, invc_bins);
		H.DefineTH2D( "PhiQuadrantSpectra", "constcosth_UVquad"+Quad[i],				"Michel spectrum UV quadrant "+Quad[i]+" for ECal, binned in cos#theta",	nxbins, xmin,xmax,100,-1.,1.);

		H.DefineTH2D( "PhiQuadrantSpectra", "Spectrum_XYquad"+Quad[i], 		"Michel spectrum XY quadrant "+Quad[i]+";Momentum [MeV];cos(#theta)",				nxbins_m,xmin_m,xmax_m,ncbins_m, -1.0, 1.0);
		H.DefineTH2D( "PhiQuadrantSpectra", "Spectrum_PzVsPt_XYquad"+Quad[i],	"Longitudinal vs transverse momentum spectrum XY quadrant "+Quad[i]+";P_{t} [MeV];P_{z} [MeV]",	nxbins_m,xmin_m,xmax_m,2*nxbins_m, -1.*xmax_m, xmax_m);
		H.Sumw2("Spectrum_XYquad"+Quad[i]);
		H.Sumw2("Spectrum_PzVsPt_XYquad"+Quad[i]);
		H.DefineTH2D_Yvarwidth( "PhiQuadrantSpectra", "constth_XYquad"+Quad[i],		"Michel spectrum XY quadrant "+Quad[i]+" for ECal, binned in #theta",		nxbins, xmin,xmax, c_bins);
		H.DefineTH2D_Yvarwidth( "PhiQuadrantSpectra", "constinvcosth_XYquad"+Quad[i],	"Michel spectrum XY quadrant "+Quad[i]+" for ECal, binned in 1/cos#theta",	nxbins, xmin,xmax, invc_bins);
		H.DefineTH2D( "PhiQuadrantSpectra", "constcosth_XYquad"+Quad[i],				"Michel spectrum XY quadrant "+Quad[i]+" for ECal, binned in cos#theta",	nxbins, xmin,xmax,100,-1.,1.);

		H.DefineTH2D("PhiQuadrantSpectra", "u0vsv0_quad"+Quad[i], "UV position of the center of the helix in quad "+Quad[i], 320, -16.,16., 320, -16.,16.);
		// H.DefineTH2D("PhiQuadrantSpectra", "x0vsy0_quad"+Quad[i], "XY position of the center of the helix in quad "+Quad[i], 320, -16.,16., 320, -16.,16.);
		H.DefineTH1D("PhiQuadrantSpectra","Pt_UVquad"+Quad[i],"Transverse momentum for UV quadrant "+Quad[i],1200,10,50);
		H.DefineTH1D("PhiQuadrantSpectra","Chi2_UVquad"+Quad[i],"chisquare for UV quadrant "+Quad[i],1000,0,10);
	}
	

	// H.DefineTH2D("PhiQuadrantSpectra", "x0vsy0", "XY position of the center of the helix", 320, -16.,16., 320, -16.,16.);


	Log->info( "Register PhiQuadrantSpectra spectrum ");
	return true;
}

bool PhiQuadrantSpectra::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
	accflag = 1; // Default value for data
	
	if (E.Exists("micheld_accflag") && E.Exists("ndecays"))
	// Extract the accflag for MC
	{
		// To store later in the output file
		E.cumulative_accflag |= E.micheld_accflag[0];
		// MC runs with the accflag for the derivatives mostly
		if (! E.ndecays > 0)
		{
			// If there is no decay from the micheld the weight is 0. (Probably still have an entry though.)
			accflag = 0;
		}
		else
		{
			try
			{
				accflag = MichelWeight(E.micheld_accflag[0]);
			}
			catch( const char* Msg)
			{
				Log->crit("MichelSpectrum: %s = 0x%x",Msg,E.micheld_accflag[0]);
				exit(1);
			}
		}
		
	}
	// The selected track
	Trk = E.seltrack[0];
	
	H.Fill("Spectrum_UVquad"+Quad[E.hefit_uvquad[Trk]],		E.ptot[Trk],E.costh[Trk], accflag);
	H.Fill("Spectrum_PzVsPt_UVquad"+Quad[E.hefit_uvquad[Trk]],E.pt[Trk],E.hefit_pz[Trk], accflag);
	H.Fill("Pt_UVquad"+Quad[E.hefit_uvquad[Trk]],		E.pt[Trk], accflag);
	H.Fill("Chi2_UVquad"+Quad[E.hefit_uvquad[Trk]],		E.hefit_chi2[Trk]/E.hefit_ndof[Trk], accflag);

	H.Fill("constth_UVquad"+Quad[E.hefit_uvquad[Trk]],		E.ptot[Trk], E.costh[Trk], accflag);
	H.Fill("constcosth_UVquad"+Quad[E.hefit_uvquad[Trk]],		E.ptot[Trk], E.costh[Trk], accflag);
	H.Fill("constinvcosth_UVquad"+Quad[E.hefit_uvquad[Trk]],	E.ptot[Trk], E.costh[Trk], accflag);
	H.Fill("u0vsv0_quad"+Quad[E.hefit_uvquad[Trk]]			,E.hefit_ucenter[Trk],E.hefit_vcenter[Trk]);

	
	// H.Fill("x0vsy0",x0,y0, accflag);
	H.Fill("Spectrum_XYquad"+Quad[E.hefit_xyquad[Trk]],		E.ptot[Trk],E.costh[Trk], accflag);
	H.Fill("Spectrum_PzVsPt_XYquad"+Quad[E.hefit_xyquad[Trk]],E.pt[Trk],E.hefit_pz[Trk], accflag);

	H.Fill("constth_XYquad"+Quad[E.hefit_xyquad[Trk]],		E.ptot[Trk], E.costh[Trk], accflag);
	H.Fill("constcosth_XYquad"+Quad[E.hefit_xyquad[Trk]],		E.ptot[Trk], E.costh[Trk], accflag);
	H.Fill("constinvcosth_XYquad"+Quad[E.hefit_xyquad[Trk]],	E.ptot[Trk], E.costh[Trk], accflag);

	return true;
}

