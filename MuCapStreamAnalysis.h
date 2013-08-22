// Andrei Gaponenko, 2013

#ifndef MuCapStreamAnalysis_h
#define MuCapStreamAnalysis_h

#include <string>

#include "EventClass.h"
#include "WireCluster.h"
#include "TimeWindow.h"
#include "TDCHitWP.h"

#include "HistTDCWidth.h"
//#include "HistStreamAnalysis.h"
#include "MuCapUVAnalysis.h"
#include "MuCapContainmentCheck.h"
#include "HistMuCapTruth.h"
#include "HistTDCParticleClassifier.h"
#include "HistHitStream.h"
#include "HistOccupancy.h"
#include "HistDriftTime.h"
#include "HistMuCapFinal.h"

#include "Math/Point2D.h"
#include "TAxis.h"
class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;

//================================================================
class MuCapStreamAnalysis {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_NOHITS, "No hits");
    ax->SetBinLabel(1+CUT_BEAM_VETO, "Beam veto");
    ax->SetBinLabel(1+CUT_WINTIME, "win time");
    ax->SetBinLabel(1+CUT_MULTIWIN_NEXTDT, "multiwin time");

    ax->SetBinLabel(1+CUT_Z_CONTAINED, "Z contained");

    ax->SetBinLabel(1+CUT_PC7_HIT, "PC7 hit");
    ax->SetBinLabel(1+CUT_PC7_COORDINATE, "PC7 V"); // uncomment once implemented

    ax->SetBinLabel(1+CUTS_LOOSE_PROTONS, "Loose protons");
  }

public:
  enum EventCutNumber {
    CUT_NOHITS,

    CUT_BEAM_VETO, // veto PC1-4 here to make wintime smooth (get rid of cyclotron structure)

    CUT_MULTIWIN_NEXTDT, // plot t2 and t2 - t1 before

    CUT_WINTIME,

    // === UV analysis goes here ===

    CUT_Z_CONTAINED,

    // Just do the downstream analysis, but at earlier times
    // Can use the already-studies TDC width cut to select protons
    // Can allow some upstream hits...  Note that upstream DC hits
    // would be eaten by the muon window for early times anyway.
    CUT_PC7_HIT,
    CUT_PC7_COORDINATE,

    CUTS_LOOSE_PROTONS,

    // plot lastPlane distribution here
    // plot number of ranges
    // plot begin-end of ranges
    // plot begin-end of holes
    // plot missing planes
    // Plot TDC widths/particle ID distributions

    // // Per 20130501-muminus/slides.pdf
    // // the Rext cut can only marginally improve proton purity
    // // at the cost of a large drop in efficiency.
    // // Do not use the "tight" protons for now.
    // // Glen argues that we need events with protons
    // // ranging out in the stack and not hitting the glass
    // // to estimate their energy spectrum.  However even
    // // hitting the glass protons give information about
    // // their MIN possible energy.
    // CUT_MIN_RANGE,
    // CUT_RANGE_GAPS, // in the "stream" direction
    // CUT_REXT,
    // CUTS_TIGHT_PROTONS,

    CUTS_END
  };

  void init(HistogramFactory &hf, const std::string& hdir,
            const DetectorGeo& geom, const ConfigFile &conf,
            TimeWindow::StreamType cutWinStream, double cutWinTimeMin);

  void process(const EventClass& evt,
               const TimeWindowingResults& wres,
               const ROOT::Math::XYPoint& muStopUV,
               const std::vector<ClustersByPlane>& afterTrigGlobalClusters);

  MuCapStreamAnalysis()
    : doMCTruth_(false)
    , cutStream_(TimeWindow::DOWNSTREAM)
    , cutBeamVetoMaxPCplanes_()
    , cutWinTimeMin_()
    , cutWinTimeMax_()
    , cutMultiwinNextdt_()
    , cutZContainedNumToCheck_()
    , cutZContainedMaxHitPlanes_()
    , cutRextMax_()
    , geom_()
    , h_cuts_r()
    , h_cuts_p()
    , hBeamVetoNumHitPlanes_()
    , hHitPCsAterBeamVeto_()
    , hWindowTimeBefore_()
    , hWindowTimeAfter_()
    , hNumAfterTrigWindows_()
    , hWindow2Time_()
    , hWindow2dt_()
    , hZContaintedNumHitPlanesUp_()
    , hZContaintedNumHitPlanesDn_()
    , hNumPC7Clusters_()
    , hNumPC7WiresVsClusters_()
    , hPC7DistanceToMuStop_()
    , hLastPlaneLoose_()
  {}

private :
  bool doMCTruth_;
  TimeWindow::StreamType cutStream_;

  int cutBeamVetoMaxPCplanes_;

  double cutWinTimeMin_;
  double cutWinTimeMax_;

  double cutMultiwinNextdt_; // min t2-t1

  // The minimal required number of planes
  // at each end to check for "z containment".
  int cutZContainedNumToCheck_;
  //
  int cutZContainedMaxHitPlanes_; // in the range on each side

  double cutPC7MaxDistanceToMuStop_;

  double cutRextMax_;

  const DetectorGeo *geom_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1 *hBeamVetoNumHitPlanes_;
  TH1 *hHitPCsAterBeamVeto_; // occupancy after the cut

  TH1 *hWindowTimeBefore_;
  TH1 *hWindowTimeAfter_;

  TH1 *hNumAfterTrigWindows_;
  TH1 *hWindow2Time_;
  TH1 *hWindow2dt_;

  TH1 *hZContaintedNumHitPlanesUp_;
  TH1 *hZContaintedNumHitPlanesDn_;

  TH1 *hNumPC7Clusters_;
  TH1 *hNumPC7WiresVsClusters_;

  TH1 *hPC7DistanceToMuStop_;

  TH1 *hLastPlaneLoose_;

  MuCapUVAnalysis uvan_;
  HistPlaneRanges hRangeDIO_;
  HistDriftTime hdriftPCFiltered_;

  HistHitStream hhsZContained_;
  HistHitStream hhsLooseProtons_;

  HistPlaneRanges hRangeAfterPC7Cuts_;

  HistTDCWidth hwidthPCDIO_;
  HistTDCWidth hwidthDCDIO_;
  HistTDCWidth hwidthPCLooseProtons_;
  HistTDCWidth hwidthDCLooseProtons_;

  HistTDCParticleClassifier hcdio_;
  HistTDCParticleClassifier hcLooseProtons_;

  HistMuCapFinal hfLoose_;

  HistMuCapTruth htruthLoose_;

  EventCutNumber analyze(const EventClass& evt,
                         const TimeWindowingResults& wres,
                         const ROOT::Math::XYPoint& muStopUV,
                         const std::vector<ClustersByPlane>& afterTrigGlobalClusters);
};

#endif/*MuCapStreamAnalysis_h*/
