//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Include de C++
#include <math.h>
// using namespace std;

// Include Program
#include "TreeClass.h"

// Modules

#include "TCAPm12widthCut.C"
#include "TECCut.C"
#include "EventTypeCut.C"
#include "MuLastPCut.C"
#include "MuUVCut.C"
#include "DkWinTimeCut.C"
#include "NtracksCut.C"
#include "IerrorCut.C"
#include "StartStopCut.C"
#include "ChargeCut.C"
#include "PairMatchingCut.C"
#include "DistToTargetSel.C"
#include "Mu_eVertexSel.C"
#include "Mu_eVertexCut.C"
#include "DkFitTimeCut.C"
#include "PACTCut.C"

#include "ChiSquare.C"
#include "MomAngSmearing.C"

#include "AsymmetryPlots.C"

#include "MichelSpectrum.C"
#include "ECalibHistograms.C"
#include "PhiQuadrantSpectra.C"

#include "FiducialCut.C"

#include "StatusHistograms.C"
#include "BiasPlots.C"
#include "GlobalHistograms.C"

#include "ConfigFile.h"
#include "MuCapture.h"

void LoadAnalysisClasses(const ConfigFile& conf, TreeClass *AnaObj)
{
  const bool doMuCapture = conf.read<bool>("MuCapture/Do");
  const bool doDefaultTWIST = !doMuCapture || conf.read<bool>("MuCapture/doDefaultTWIST");

  if(doMuCapture) {
    AnaObj->Register(new MuCapture());
  }
  if(!doDefaultTWIST) {
    return;
  }

	AnaObj->Register( new StatusHistograms("Beginning",	"at the beginning of the treesumming"));

	AnaObj->Register( new TCAPm12widthCut());
	AnaObj->Register( new TECCut());
	AnaObj->Register( new EventTypeCut());
	AnaObj->Register( new MuLastPCut());
	AnaObj->Register( new MuUVCut());
	AnaObj->Register( new PACTCut());
	AnaObj->Register( new DkWinTimeCut());
	AnaObj->Register( new NtracksCut());
	AnaObj->Register( new IerrorCut());
	AnaObj->Register( new StartStopCut());
	AnaObj->Register( new ChargeCut());
	AnaObj->Register( new PairMatchingCut());
	AnaObj->Register( new Mu_eVertexCut());
	AnaObj->Register( new DistToTargetSel());
	AnaObj->Register( new Mu_eVertexSel());
	AnaObj->Register( new DkFitTimeCut());

	AnaObj->Register( new ChiSquare());
	AnaObj->Register( new MomAngSmearing());
	AnaObj->Register( new PhiQuadrantSpectra());

	AnaObj->Register( new AsymmetryPlots());

	AnaObj->Register( new MichelSpectrum("Selected", "of the selected events"));
	AnaObj->Register( new ECalibHistograms());
	AnaObj->Register( new StatusHistograms("InSpectrum",	"of the events in the spectrum"));
	AnaObj->Register( new BiasPlots("InSpectrum",	"of the events in the spectrum"));

	AnaObj->Register( new FiducialCut());
	AnaObj->Register( new MichelSpectrum("Fiducial", "after the fiducial cut"));
	AnaObj->Register( new StatusHistograms("InFiducial",	"of the events in the fiducial"));
	AnaObj->Register( new BiasPlots("InFiducial",	"of the events in the fiducial"));

	AnaObj->Register( new GlobalHistograms());
}
