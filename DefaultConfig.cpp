#include "ConfigFile.h"

template<class T>
void SetKey(ConfigFile &Conf, string key, const T& t)
{
	if (! Conf.keyExists( key))
	{
		Conf.add(key, t);
	}
}

void SetDefault(ConfigFile &Conf)
{
	// ================ Standard Analysis =============== //

	//	TCAPm12widthCut
	SetKey(Conf, "TCAPm12widthCut/Do", true);
	SetKey(Conf, "TCAPm12widthCut/TCAPMin", 0.0);
	SetKey(Conf, "TCAPm12widthCut/TCAPMax", 200.0);
	SetKey(Conf, "TCAPm12widthCut/m12widthMin", 0.0);
	SetKey(Conf, "TCAPm12widthCut/m12widthMax", 20000.0);

	//	EventTypeCut
	SetKey(Conf, "EventTypeCut/Do", true);
	SetKey(Conf, "EventTypeCut/Types", "1,2,6,7,10,11,21,22");
	
	//	MuLastPCut
	SetKey(Conf, "MuLastPCut/Do", true);
	SetKey(Conf, "MuLastPCut/LastPlane", 28);
	
	//	MuUVCut
	SetKey(Conf, "MuUVCut/Do", true);
	SetKey(Conf, "MuUVCut/MaxRadius", 2.5);
	
	//	DkWinTimeCut
	SetKey(Conf, "DkWinTimeCut/Do", true);
	SetKey(Conf, "DkWinTimeCut/Min", -100000.0);
	SetKey(Conf, "DkWinTimeCut/Max", 100000.0);
	
	//	NtracksCut
	SetKey(Conf, "NtracksCut/Do", true);
	
	//	IerrorCut
	SetKey(Conf, "IerrorCut/Do", true);
	SetKey(Conf, "IerrorCut/SelectIerror", "0");
	
	//	StartStopCut
	SetKey(Conf, "StartStopCut/Do", true);

	//	ChargeCut
	SetKey(Conf, "ChargeCut/Do", true);
	SetKey(Conf, "ChargeCut/SelectCharge", 1);

	//	PairMatchingCut
	SetKey(Conf, "PairMatchingCut/Do", true);
	SetKey(Conf, "PairMatchingCut/MinCDA", 0.5);
	SetKey(Conf, "PairMatchingCut/MinDT", 60.0);
	SetKey(Conf, "PairMatchingCut/tolerance_abs", 1.e-4); /*cm*/
	SetKey(Conf, "PairMatchingCut/tolerance_rel", 0.);
	SetKey(Conf, "PairMatchingCut/max_iter", 100);
	
	//	DistToTargetSel
	SetKey(Conf, "DistToTargetSel/Do", true);
	
	//	Mu_eVertexSel
	SetKey(Conf, "Mu_eVertexSel/Do", true);
	SetKey(Conf, "Mu_eVertexSel/Do_mu_e_Calc-MOFIA", 0);
	
	//	DkFitTimeCut
	SetKey(Conf, "DkFitTimeCut/Do", true);
	SetKey(Conf, "DkFitTimeCut/Min", 1050.0);	// in ns
	SetKey(Conf, "DkFitTimeCut/Max", 9000.0);	// in ns
	
	//	PACTCut
	SetKey(Conf, "PACTCut/Do", 0);
	SetKey(Conf, "PACTCut/NbBins", 500);
	SetKey(Conf, "PACTCut/MaximumWidth", 500.0);	// ns
	SetKey(Conf, "PACTCut/quadrant_11", 1);
	SetKey(Conf, "PACTCut/quadrant_12", 1);
	SetKey(Conf, "PACTCut/quadrant_21", 1);
	SetKey(Conf, "PACTCut/quadrant_22", 1);
	SetKey(Conf, "PACTCut/region_11", "1.634, -5.428, -61.26, 642.3");
	SetKey(Conf, "PACTCut/region_12", "1.634, -5.428, -61.26, 642.3");
	SetKey(Conf, "PACTCut/region_21", "1.634, -5.428, -61.26, 642.3");
	SetKey(Conf, "PACTCut/region_22", "1.634, -5.428, -61.26, 642.3");
	
	// ==================== Special Analyses ======================== //
	
	// AsymmetryPlots
	SetKey(Conf, "AsymmetryPlots/WeightedPlots", 1); // default ON
	SetKey(Conf, "AsymmetryPlots/UnweightedPlots", 0); // default OFF
	SetKey(Conf, "AsymmetryPlots/tmu", 2197.03);			// ns
	SetKey(Conf, "AsymmetryPlots/t_min", 245.28);		// ns
	SetKey(Conf, "AsymmetryPlots/t_width_min", 68.435);	// ns
	SetKey(Conf, "AsymmetryPlots/Fiducial/costh_min", 	0.50);
	SetKey(Conf, "AsymmetryPlots/Fiducial/costh_max", 	0.94);
	SetKey(Conf, "AsymmetryPlots/Fiducial/ptot_min", 	31.0);
	SetKey(Conf, "AsymmetryPlots/Fiducial/ptot_max", 	51.5);
	SetKey(Conf, "AsymmetryPlots/Fiducial/long_min", 	13.7);
	SetKey(Conf, "AsymmetryPlots/Fiducial/trans_max", 	43.0);
	SetKey(Conf, "AsymmetryPlots/Weighting", 2.718);

	
	// ================== General Parameters ======================== //
	// Kinematic end point
	SetKey(Conf, "Parameters/KinematicPmax",		 	52.828);

	// Speed of light
	SetKey(Conf, "Parameters/c",29.979245800);	// cm/ns

	// Michel spectrum histogram
	SetKey(Conf, "Parameters/NCosThBinsMichel", 110);
	SetKey(Conf, "Parameters/NXBinsMichel", 110);
	SetKey(Conf, "Parameters/XMinMichel", 0.0);
	SetKey(Conf, "Parameters/XMaxMichel", 55.0);
}
