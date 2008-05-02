//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class AsymmetryPlots : public ModuleClass{
	public :
		AsymmetryPlots()		{};
		~AsymmetryPlots()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;

		double Weighting;
		double costh_min;
		double costh_max;
		double ptot_min;
		double ptot_max;
		double plong_min;
		double ptrans_max;
		double kpmax;
		bool WeightedPlots;
		bool UnweightedPlots;
};

bool AsymmetryPlots::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	//    -------- Name of the cut ---------     //
	Log	= TmpLog;

	
	WeightedPlots	= Conf.read<bool>("AsymmetryPlots/WeightedPlots"); // default ON
	UnweightedPlots= Conf.read<bool>("AsymmetryPlots/UnweightedPlots"); // default OFF
	if ( (! WeightedPlots) && (! UnweightedPlots) )
	{
		Log->info( "Asymmetry Plots turned OFF");
		return false ;
	}

	//	 --------- Parameters initialization ---------		//
	double tmu			= Conf.read<double>("AsymmetryPlots/tmu");			// ns
	// lower edge for first constant-stats bin
	double t_min		= Conf.read<double>("AsymmetryPlots/t_min");		// ns
	// width of first constant-stats bin
	double t_width_min	= Conf.read<double>("AsymmetryPlots/t_width_min");	// ns

	costh_min	= Conf.read<double>("AsymmetryPlots/Fiducial/costh_min");
	costh_max	= Conf.read<double>("AsymmetryPlots/Fiducial/costh_max");
	ptot_min	= Conf.read<double>("AsymmetryPlots/Fiducial/ptot_min");
	ptot_max	= Conf.read<double>("AsymmetryPlots/Fiducial/ptot_max");
	plong_min	= Conf.read<double>("AsymmetryPlots/Fiducial/long_min");
	ptrans_max	= Conf.read<double>("AsymmetryPlots/Fiducial/trans_max");
	kpmax		= Conf.read<double>("Parameters/KinematicPmax");
    
	Weighting	= Conf.read<double>("AsymmetryPlots/Weighting");

	//	 __________ setup variable width bins  __________

	// if (WeightedPlots)
	// {
		vector<double> Vbins;
		Vbins.push_back(0.0);
		Vbins.push_back(t_min);
		Vbins.push_back(t_min+t_width_min);
	
		double e;
		int i_bin = 2;
		while (1)
		{
			i_bin+=1;
			e = 2.*exp(-1.*Vbins[i_bin-1]/tmu) - exp(-1*Vbins[i_bin-2]/tmu); 
			if (e<0) break;
			Vbins.push_back( -1.*tmu*log(e) );
		}
		Vbins[Vbins.size()-1] = 9000.;  // it's approximately 9000ns anyway
		Vbins.push_back(10000.);
	
		//	 --------- Histograms initialization ---------		//
		// first the bins that have constant stat widths
		H.DefineTH1D_varwidth("AsymmetryPlots", "h_gnt_u",      "U_n vs Decay Time (ns)"     ,Vbins);
		H.DefineTH1D_varwidth("AsymmetryPlots", "h_gnt_d",      "D_n vs Decay Time (ns)"     ,Vbins);
		H.DefineTH1D_varwidth("AsymmetryPlots", "h_gnt_deriv1", "Deriv_1 vs Decay Time (ns)" ,Vbins);
		H.DefineTH1D_varwidth("AsymmetryPlots", "h_gnt_deriv2", "Deriv_2 vs Decay Time (ns)" ,Vbins);
		H.DefineTH1D_varwidth("AsymmetryPlots", "h_gnt",        "Gn vs Decay Time (ns)"      ,Vbins);
	// }

	Log->info( "Register Asymmetry Plots");
	return true;
}

bool AsymmetryPlots::Process(EventClass &E, HistogramFactory &H)
{
	// __________ fill as long as fiducial is satisfied __________
	int trk = E.seltrack[0];
	if (WeightedPlots)
	{
		if ( costh_min < fabs(E.costh[trk]) < costh_max and ptot_min < E.ptot[trk] < ptot_max and fabs(E.pz[trk]) > plong_min and E.pt[trk] < ptrans_max)
		{
			double decay_time = E.hefit_time[trk];
			double ai = E.costh[trk]*((E.ptot[trk]/kpmax)-0.5)/(1.5-(E.ptot[trk]/kpmax));
			H.Fill("h_gnt_u",decay_time, ai * pow(fabs(ai),Weighting));
			H.Fill("h_gnt_d",decay_time, pow(fabs(ai),Weighting+1.));
			H.Fill("h_gnt_deriv1",decay_time, pow(fabs(ai),2.*Weighting+2.));
			H.Fill("h_gnt_deriv2",decay_time, ai * pow(fabs(ai),2.*Weighting+1.));
		}
	}

	return true;
}

