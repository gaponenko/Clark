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
	// MuCapture, if enabled, can prevent the standard analysis from running
	SetKey(Conf, "MuCapture/Do", false);
	SetKey(Conf, "MuCapture/MCType", "G3");
        SetKey(Conf, "AnalyzeAllEvents", false); // Michel analysis code crashes on empty events

	// ================ Standard Analysis =============== //

	//	Energy Calibration
	SetKey(Conf, "EnergyCalibration/Mode", 0);
	//// This is turned on or off using the command line
	//// The possible modes:
	////	0 => alpha and beta corrections = shift
	////	1 => 1)beta correction = scale; 2)alpha correction = shift
	////	2 => 1)alpha correction = shift; 2)beta correction = scale
	////	3 => alpha and beta corrections = scale
	
	// Momentum and angle smearing
	SetKey(Conf, "MomAngSmearing/Do", false);
	SetKey(Conf, "MomAngSmearing/MomMean", 0.0);
	SetKey(Conf, "MomAngSmearing/MomSigma", 0.0);
	SetKey(Conf, "MomAngSmearing/AngMean", 0.0);
	SetKey(Conf, "MomAngSmearing/AngSigma", 0.0);

	// Dimension scaling
	SetKey(Conf, "MomentumScale/Do", false);
	//// This scales the momentum values pz or pu and pv
	//// If MomentumScale/Pt = S
	//// Then pu=pu*S and pv=pv*S
	SetKey(Conf, "MomentumScale/Pt", 1.0);
	SetKey(Conf, "MomentumScale/Pz", 1.0);
	SetKey(Conf, "MomentumScale/Pt_Upstream", 1.0);
	SetKey(Conf, "MomentumScale/Pz_Upstream", 1.0);
	SetKey(Conf, "MomentumScale/Pt_Downstream", 1.0);
	SetKey(Conf, "MomentumScale/Pz_Downstream", 1.0);
	

	//	TCAPm12widthCut
	SetKey(Conf, "TCAPm12widthCut/Do", true);
	//// For 2006 data: TCAPMin = 26; TCAPMax = 52
	//// For 2007 data: TCAPMin = 29; TCAPMax = 55
	SetKey(Conf, "TCAPm12widthCut/TCAPMin", -10000.0);
	SetKey(Conf, "TCAPm12widthCut/TCAPMax", 10000.0);
	SetKey(Conf, "TCAPm12widthCut/m12widthMin", -100000.0);
	SetKey(Conf, "TCAPm12widthCut/m12widthMax", 100000.0);

	//	TECCut
	SetKey(Conf, "TECCut/Do", false);
	SetKey(Conf, "TECCut/TECxMin", -1000.0);
	SetKey(Conf, "TECCut/TECxMax", 1000.0);
	SetKey(Conf, "TECCut/TECyMin", -1000.0);
	SetKey(Conf, "TECCut/TECyMax", 1000.0);
	SetKey(Conf, "TECCut/TECtanthxMin", -1000.0);
	SetKey(Conf, "TECCut/TECtanthxMax", 1000.0);
	SetKey(Conf, "TECCut/TECtanthyMin", -1000.0);
	SetKey(Conf, "TECCut/TECtanthyMax", 1000.0);
	SetKey(Conf, "TECCut/TECsigmaxMin", -1.0);
	SetKey(Conf, "TECCut/TECsigmaxMax", 1000.0);
	SetKey(Conf, "TECCut/TECsigmayMin", -1.0);
	SetKey(Conf, "TECCut/TECsigmayMax", 1000.0);
	SetKey(Conf, "TECCut/TECifailxMin", -1.0);
	SetKey(Conf, "TECCut/TECifailxMax", 1000.0);
	SetKey(Conf, "TECCut/TECifailyMin", -1.0);
	SetKey(Conf, "TECCut/TECifailyMax", 1000.0);

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
	//// The module selects the tracks on the right side
	//// of the target according to the classification.
	//// The StartPlane and StopPlane can be defined
	//// precisely. For example for the tracks starting
	//// right at the target and using 22 planes:
	//// StartStopCut/StartPlane = 22,23
	//// StartStopCut/StopPlane = 1,44
	SetKey(Conf, "StartStopCut/Do", true);
	SetKey(Conf, "StartStopCut/StartPlane", "-1,-1");
	SetKey(Conf, "StartStopCut/StopPlane", "-1,-1");

	//	ChargeCut
	SetKey(Conf, "ChargeCut/Do", true);
	SetKey(Conf, "ChargeCut/SelectCharge", 1);

	//	PairMatchingCut
	SetKey(Conf, "PairMatchingCut/Do", true);
	SetKey(Conf, "PairMatchingCut/MinCDA", 0.5);				// For all the other tracks
	SetKey(Conf, "PairMatchingCut/MinCDA_BrokTrk", 2.0);		// For broken tracks
	SetKey(Conf, "PairMatchingCut/MinCDA_Beame", 0.5);			// For beam positron
	SetKey(Conf, "PairMatchingCut/BrokTrk_Z_Range", "22.0,48.0");	// For beam positron
	SetKey(Conf, "PairMatchingCut/Beame_Z_Range", "0.0,6.0");	// For beam positron
	SetKey(Conf, "PairMatchingCut/MinDT", 60.0);
	SetKey(Conf, "PairMatchingCut/tolerance_abs", 1.e-4); /*cm*/
	SetKey(Conf, "PairMatchingCut/tolerance_rel", 0.);
	SetKey(Conf, "PairMatchingCut/max_iter", 100);
	SetKey(Conf, "PairMatchingCut/CalculateCDA", false);
	
	//	Mu_eVertexCut
	SetKey(Conf, "Mu_eVertexCut/Do", true);
	SetKey(Conf, "Mu_eVertexCut/Do_mu_e_Calc-MOFIA", 0);
	//// Function 0:
	////     min < r < max
	////     CutParameters = min,max
	//// Function 1:
	////     Bmin + Amin / |cos(theta)| < r < Bmax + Amax / |cos(theta)|
	////     CutParameters = Bmin,Amin,Bmax,Amax
	SetKey(Conf, "Mu_eVertexCut/CutFunction", 1);
	SetKey(Conf, "Mu_eVertexCut/CutParameters", "0.0,0.0,1.0,1.0");
	
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
	
	// Fiducial Cut
	SetKey(Conf, "FiducialCut/MinAbsCosTheta", 0.54);
	SetKey(Conf, "FiducialCut/MaxAbsCosTheta", 0.96);
	SetKey(Conf, "FiducialCut/MaxMomentum", 52.0);			// Mev/c
	SetKey(Conf, "FiducialCut/MinTransMom", 10.0);			// Mev/c
	SetKey(Conf, "FiducialCut/MaxTransMom", 38.0);			// Mev/c
	SetKey(Conf, "FiducialCut/MinLongiMom", 14.0);			// Mev/c


	//	StatusHistograms
	SetKey(Conf, "StatusHistograms/Do", true);
	SetKey(Conf, "StatusHistograms/PerWindowType", false);

	// ==================== Special Analyses ======================== //
	
	// Truth Bank
	SetKey(Conf, "TruthBank/Do",			false);
	//// This flag turns the truth bank analysis off. No vertex bank analysed
	//// but the accflag is still used for the Michel spectrum.

	// Bias Plots
	SetKey(Conf, "BiasPlots/Do",			true);
	SetKey(Conf, "BiasPlots/Nptotbins",		8);
	SetKey(Conf, "BiasPlots/minptot",		15.0);			// MeV
	SetKey(Conf, "BiasPlots/maxptot",		55.0);			// MeV

	SetKey(Conf, "BiasPlots/Nbinscosth",	100);			//
	SetKey(Conf, "BiasPlots/mincosth",		0.3);			//
	SetKey(Conf, "BiasPlots/maxcosth",		-0.3);			//

	SetKey(Conf, "BiasPlots/Nbinsdp",		401);			//
	SetKey(Conf, "BiasPlots/mindp",			-2.0);			// MeV
	SetKey(Conf, "BiasPlots/maxdp",			2.0);			// MeV
	
	SetKey(Conf, "BiasPlots/Nbinsdcosth",	200);			//
	SetKey(Conf, "BiasPlots/mindcosth",		-0.05);			//
	SetKey(Conf, "BiasPlots/maxdcosth",		0.05);			//
	

	// Asymmetry Plots
	SetKey(Conf, "AsymmetryPlots/WeightedPlots",	1); // default ON
	SetKey(Conf, "AsymmetryPlots/UnweightedPlots",	0); // default OFF
	SetKey(Conf, "AsymmetryPlots/tmu",				2197.03);			// ns
	SetKey(Conf, "AsymmetryPlots/t_min",			245.28);		// ns
	SetKey(Conf, "AsymmetryPlots/t_width_min",		68.435);	// ns
	SetKey(Conf, "AsymmetryPlots/Fiducial/costh_min", 	0.54);
	SetKey(Conf, "AsymmetryPlots/Fiducial/costh_max", 	0.96);
	SetKey(Conf, "AsymmetryPlots/Fiducial/ptot_min", 	31.0);
	SetKey(Conf, "AsymmetryPlots/Fiducial/ptot_max", 	52.0);
	SetKey(Conf, "AsymmetryPlots/Fiducial/long_min", 	14.0);
	SetKey(Conf, "AsymmetryPlots/Fiducial/trans_max", 	38.0);
	SetKey(Conf, "AsymmetryPlots/Weighting", 1.0);

	// ChiSquare
	SetKey(Conf, "ChiSquare/Do", 0);				// default OFF
	SetKey(Conf, "ChiSquare/RedChi2CutMin", -1.0);			// default OFF
	SetKey(Conf, "ChiSquare/RedChi2CutMax", -1.0);			// default OFF

	// Quadrant Spectra
	SetKey(Conf, "PhiQuadrantSpectra/Do",			false);

	// Michel spectrum histogram
	SetKey(Conf, "ChiSquare/NCosThBinsMichel", 100);
	SetKey(Conf, "ChiSquare/NXBinsMichel", 110);
	SetKey(Conf, "ChiSquare/XMinMichel", 0.0);
	SetKey(Conf, "ChiSquare/XMaxMichel", 55.0);
	
	// ================== Detector Parameters ======================== //
	SetKey(Conf, "Detector/BField",	2.0);
	SetKey(Conf, "Detector/GeometryFile",	-57);	// CFM number of the file

	// ================== General Parameters ======================== //
	// Truth Bank
	SetKey(Conf, "TruthBank/MuUVMaxRadius", 	4.0);	// 
	SetKey(Conf, "TruthBank/TriggerTimeJitter", 100.0);	// 
	SetKey(Conf, "TruthBank/MuZAroundTarget",	1.0);	// 
	SetKey(Conf, "TruthBank/MueVertexEpsilon",	0.01);	// 
	SetKey(Conf, "TruthBank/MinDkTimeAfterMu",	300.0);	// 
	SetKey(Conf, "TruthBank/MaxDkTimeAfterMu",	9100.0);// 
	
	// Kinematic end point
	SetKey(Conf, "Parameters/KinematicPmax",		 	52.828);

	// Speed of light
	SetKey(Conf, "Parameters/c",29.979245800);	// cm/ns

	// Michel spectrum histogram
	SetKey(Conf, "Parameters/NCosThBinsMichel", 100);
	SetKey(Conf, "Parameters/NXBinsMichel", 110);
	SetKey(Conf, "Parameters/XMinMichel", 0.0);
	SetKey(Conf, "Parameters/XMaxMichel", 55.0);

	// Energy Calibration histograms
	SetKey(Conf, "Parameters/NCBins_EC",	360);
	SetKey(Conf, "Parameters/NCBins_ECinvc",72);
	SetKey(Conf, "Parameters/NXBins_EC",	4000);
	SetKey(Conf, "Parameters/XMin_EC",		30.0);
	SetKey(Conf, "Parameters/XMax_EC",		70.0);
	SetKey(Conf, "Parameters/CosMin_EC",	0.3);
	SetKey(Conf, "Parameters/CosMax_EC",	0.95);

	SetKey(Conf, "Parameters/MuonWinType",	1);

	SetKey(Conf, "Parameters/UpDkWinType",		"2,7,9,14");
	SetKey(Conf, "Parameters/DownDkWinType",	"3,8,10,15");
	SetKey(Conf, "Parameters/DCplanes",		"5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52");
	SetKey(Conf, "Parameters/PCplanes",		"1, 2, 3, 4, 27, 28, 29, 30, 53, 54, 55, 56");
	SetKey(Conf, "Parameters/PlaneType",	"PC,PC,PC,PC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,PC,PC,PC,PC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,DC,PC,PC,PC,PC");
	SetKey(Conf, "Parameters/DCzposition",	"-49.7924,-49.3923,-48.9923,-48.5922,-48.1921,-47.7921,-47.3920,-46.9919,-42.1920,-41.7921,-36.9923,-36.5924,-29.7941,-29.3941,-22.5960,-22.1961,-15.3980,-14.9981,-10.1987, -9.7988, -4.9995, -4.5995,  4.5998,  4.9998,  9.7994, 10.1993, 14.9976, 15.3976, 22.1960, 22.5959, 29.3943, 29.7942, 36.5927, 36.9926, 41.7924, 42.1924, 46.9924, 47.3923, 47.7922, 48.1922, 48.5922, 48.9922, 49.3922, 49.7922");


        // ================== MuCapture analysis ======================== //
        SetKey(Conf, "MuCapture/loadDefaultTWISTVars", true);
        SetKey(Conf, "MuCapture/fillXtalkPC", true);
        SetKey(Conf, "MuCapture/fillXtalkDC", false);
        SetKey(Conf, "MuCapture/debugEventList", "");
        SetKey(Conf, "Detector/Geometry/dc_ppc", 25);
        SetKey(Conf, "Detector/Geometry/pc_ppc", 5);

        // The value in dt_geo.00061: -0.0904, however this is the
        // 25um foil position, and 71um target is glued on the "-z"
        // side (in G3's PC6 volume).
        // -0.0904 - 0.0250/2 - 0.071/2 = -0.1384 mm = -0.01384 cm
        SetKey(Conf, "Detector/Geometry/zTargetCenter", -0.01384);

        SetKey(Conf, "MuCapture/killPC6DeadWire", true);

        SetKey(Conf, "MuCapture/HitPreproc/PC/applyXTalk", true);
        SetKey(Conf, "MuCapture/HitPreproc/DC/applyXTalk", true);

        SetKey(Conf, "MuCapture/HitPreproc/PC/processor", "SameWireHitDiscarder");
        SetKey(Conf, "MuCapture/HitPreproc/PC/NarrowHitDiscarder/cutMinTDCWidth", 40.);
        SetKey(Conf, "MuCapture/HitPreproc/PC/SameWireHitDiscarder/cutSameWireDt", 200.);

        SetKey(Conf, "MuCapture/HitPreproc/DC/processor", "NarrowHitDiscarder");
        SetKey(Conf, "MuCapture/HitPreproc/DC/NarrowHitDiscarder/cutMinTDCWidth", 50.);

        // If defined, ignore all other input events
        SetKey(Conf, "MuCapture/inputEventNumberFile", "");
        // Print out details about events on this list
        SetKey(Conf, "MuCapture/debugEventList", "");

        SetKey(Conf, "MuCapture/winPCLength", 100.);
        // Trigger window has to be within max dt from 0
        SetKey(Conf, "MuCapture/winTrigMaxdt", 30.);
        SetKey(Conf, "MuCapture/winPCPreTrigSeparation", 1050.);

        SetKey(Conf, "MuCapture/cutMuonFirstPlane", 1);

        SetKey(Conf, "MuCapture/Accidentals/tmax", -1100.);
        // Cyclotron RF=23.058 MHz [Glen Marshall in the 2013-06-12 muminus meeting]
        // See e.g.
        // R.~H.~M.~Gummer,
        // %``Accelerating Voltage Control and Stabilization in the TRIUMF Cyclotron,''
        // IEEE Trans.\ Nucl.\ Sci.\  {\bf 22}, 1257 (1975).
        // http://inspirehep.net/record/102708
        // http://dx.doi.org/10.1109/TNS.1975.4327860
        SetKey(Conf, "MuCapture/Accidentals/cycleLength", 1000./23.058);
        SetKey(Conf, "MuCapture/Accidentals/numCycles", 104);
        SetKey(Conf, "MuCapture/Accidentals/maxSubdivisions", 3);

        SetKey(Conf, "MuCapture/winDCStart", -100.); // W.r.t PC win start
        SetKey(Conf, "MuCapture/winDCEnd", 1050.); // W.r.t PC win start
        SetKey(Conf, "MuCapture/winDCDoHistos", true);
        SetKey(Conf, "MuCapture/muStopRMax", 2.5); // in cm

        SetKey(Conf, "MuCapture/cutWinTimeMin", 400.); // ns
        SetKey(Conf, "MuCapture/cutWinTimeMax", 10000. - 1050.); // ns, make sure DC window ends.

        SetKey(Conf, "MuCapture/PACT/slopea_11", 1.634);
        SetKey(Conf, "MuCapture/PACT/slopeb_11", -1.59);
        SetKey(Conf, "MuCapture/PACT/intercepta_11", -90);
        SetKey(Conf, "MuCapture/PACT/interceptb_11", 350);

        SetKey(Conf, "MuCapture/PACT/slopea_12", 2.941);
        SetKey(Conf, "MuCapture/PACT/slopeb_12", -0.857);
        SetKey(Conf, "MuCapture/PACT/intercepta_12", -235.3);
        SetKey(Conf, "MuCapture/PACT/interceptb_12", 300);

        SetKey(Conf, "MuCapture/PACT/slopea_21", 1.25);
        SetKey(Conf, "MuCapture/PACT/slopeb_21", -1.40);
        SetKey(Conf, "MuCapture/PACT/intercepta_21", -62.5);
        SetKey(Conf, "MuCapture/PACT/interceptb_21", 350);

        SetKey(Conf, "MuCapture/PACT/slopea_22", 1.852);
        SetKey(Conf, "MuCapture/PACT/slopeb_22", -1.0);
        SetKey(Conf, "MuCapture/PACT/intercepta_22", -148.1);
        SetKey(Conf, "MuCapture/PACT/interceptb_22", 350);

        SetKey(Conf, "MuCapture/cutBeamVetoMaxPCplanes", 1);

        SetKey(Conf, "MuCapture/DIOUp/cutMinTime", 300.);
        SetKey(Conf, "MuCapture/DIODn/cutMinTime", 300.);

        SetKey(Conf, "MuCapture/cutMultiwinNextdt", 1050.);

        SetKey(Conf, "MuCapture/HistDriftTime/cutEffTrackHitDtPC", 100.);
        SetKey(Conf, "MuCapture/HistDriftTime/cutEffTrackHitDtDC", 1500.);

        // Empty file name disables the output
        SetKey(Conf, "MuCapture/uvOutFileName", "");
        SetKey(Conf, "MuCapture/commonSkimOutFileName", "");

        // make sure the planes are not adjacent at dz=0.4cm
        SetKey(Conf, "MuCapture/ProtonWindow/Containment1D/minPlanedz", 0.5);

        //----------------
        SetKey(Conf, "MuCapture/UVAnalysis/cutTrackRmax", 6.3);//cm
        // Use kinematic cuts from the DIO analysis (2009), execpt PtMax superseded by Rmax above
        SetKey(Conf, "MuCapture/UVAnalysis/cutCosThetaMin", 0.54);
        SetKey(Conf, "MuCapture/UVAnalysis/cutCosThetaMax", 0.92);
        SetKey(Conf, "MuCapture/UVAnalysis/cutPtMin", 11.); // MeV/c
        SetKey(Conf, "MuCapture/UVAnalysis/cutPzMin", 14.); // MeV/c
        SetKey(Conf, "MuCapture/UVAnalysis/cutPtotMin", 17.5); // MeV/c
        SetKey(Conf, "MuCapture/UVAnalysis/cutPtotMax", 73.5); // MeV/c
        //
        SetKey(Conf, "MuCapture/UVAnalysis/cutTrackMuonOffset", 0.4);//cm
        // Empty file name disables the output
        //        SetKey(Conf, "MuCapture/UVAnalysis/uvOutFileName", "");

        //----------------
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutCharge", +1);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutTrackWinTimedt", 100.);//ns
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutTrackRmax", 99999.);//cm
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutCosThetaMin", 0.5);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutCosThetaMax", 0.98);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutPtMin", 11.9); // MeV/c = 2 cm R
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutPzMin", 28.4); // MeV/c = 30cm L
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutPtotMin", 0.); // MeV/c
        // NB: remember about "name helixfit helixfitmommax=..." in the KCM file.
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutPtotMax", 250.); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/pos/cutTrackMuonOffset", 1.5);//cm

        //----------------
        // The kinematic cuts for the DIO normalization sample
        // follow the TWIST DIO Phys.Rev.D.
        // Here we use ptot rather than E.  But the biggest difference
        // is our use of the wire center fit.
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutCharge", -1);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutTrackWinTimedt", 100.);//ns
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutTrackRmax",6.34);//cm = 38.0 MeV/c pt
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutCosThetaMin", 0.54);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutCosThetaMax", 0.92);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutPtMin", 11.0); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutPzMin", 14.0); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutPtotMin", 17.5); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutPtotMax", 73.5); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioNorm/cutTrackMuonOffset", 1.5);//cm

        //----------------
        // A loose set of cuts to get more DIOs for vetoing capture events
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutCharge", -1);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutTrackWinTimedt", 100.);//ns
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutTrackRmax",81.7);//cm, approx 49 MeV/c pt
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutCosThetaMin", 0.0);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutCosThetaMax", 0.98);
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutPtMin", 0.0); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutPzMin", 4.0); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutPtotMin", 0.); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutPtotMax", 73.5); // MeV/c
        SetKey(Conf, "MuCapture/TrkAnalysisHF/dioVeto/cutTrackMuonOffset", 1.5);//cm

        //----------------
        // Require the track to NOT hit the last plane of the dense stack.
        // So that it is not sensitive to hits in the outer PCs.
        SetKey(Conf, "MuCapture/dnPosTrkContainment/cutMaxPlane", 51);
        SetKey(Conf, "MuCapture/dnPosTrkContainment/cutMaxRout", 15.);

        //----------------
        SetKey(Conf, "MuCapture/channels/versions", "legacy");

        SetKey(Conf, "MuCapture/channels/legacy/contained/ccut/cutMaxPlane", 55);
        SetKey(Conf, "MuCapture/channels/legacy/contained/ccut/cutMaxRout", 15.);

        SetKey(Conf, "MuCapture/channels/legacy/numGeneratorBins", 500);
        SetKey(Conf, "MuCapture/channels/legacy/genpmin", 0.);
        SetKey(Conf, "MuCapture/channels/legacy/genpmax", 500.);

        SetKey(Conf, "MuCapture/channels/legacy/contained/xvarnbins", 88);
        SetKey(Conf, "MuCapture/channels/legacy/contained/xvarmin", 30.);
        SetKey(Conf, "MuCapture/channels/legacy/contained/xvarmax", 250.);
        SetKey(Conf, "MuCapture/channels/legacy/contained/yvarnbins", 6);
        SetKey(Conf, "MuCapture/channels/legacy/contained/yvarmin", 8.);
        SetKey(Conf, "MuCapture/channels/legacy/contained/yvarmax", 32.);

        SetKey(Conf, "MuCapture/channels/legacy/uncontained/recopnbins", 88);
        SetKey(Conf, "MuCapture/channels/legacy/uncontained/recopmin", 30.);
        SetKey(Conf, "MuCapture/channels/legacy/uncontained/recopmax", 250.);

        SetKey(Conf, "MuCapture/channels/legacy/hitbased/cwiresnbins", 200);
        SetKey(Conf, "MuCapture/channels/legacy/hitbased/cwiresmin", -0.5);
        SetKey(Conf, "MuCapture/channels/legacy/hitbased/cwiresmax", 199.5);
        SetKey(Conf, "MuCapture/channels/legacy/hitbased/cplanesnbins", 28);
        SetKey(Conf, "MuCapture/channels/legacy/hitbased/cplanesmin", 0.5);
        SetKey(Conf, "MuCapture/channels/legacy/hitbased/cplanesmax", 28.5);

        //----------------
        //# Use downstream PC planes as a veto
        //NB: need implement ganged wires to include PC9,10 in the range.
        SetKey(Conf, "MuCapture/ProtonWindow/maxPlane", 52);
        // See Andrei's slides for 2013-05-01.
        // From purity vs eff Want proton eff~=0.2,
        // from slide 3 this is <~13 cm
        SetKey(Conf, "MuCapture/ProtonWindow/RextMax", 13.);

        // This cut is inactive by default
        SetKey(Conf, "HitBasedAnalysis/maxClusterWiresFilterCutPC", 999);

        //================================================================
}
