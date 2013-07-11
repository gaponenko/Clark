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
#include "HistMuCapRTruth.h"
#include "HistMuCapTruth.h"
#include "HistTDCParticleClassifier.h"
#include "HistHitStream.h"
#include "HistOccupancy.h"
#include "HistDriftTime.h"

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
    // ax->SetBinLabel(1+CUT_PC7_COORDINATE, "PC7 V"); // uncomment once implemented

    ax->SetBinLabel(1+CUTS_LOOSE_PROTONS, "Loose protons");

    //ax->SetBinLabel(1+CUT_MIN_RANGE, "Min range");
    //ax->SetBinLabel(1+CUT_RANGE_GAPS, "RANGE_GAPS");
    //ax->SetBinLabel(1+CUT_REXT, "Rext");

    ax->SetBinLabel(1+CUTS_TIGHT_PROTONS, "Tight protons");
  }

public:
  enum EventCutNumber {
    CUT_NOHITS,

    CUT_BEAM_VETO, // veto PC1-4 here to make wintime smooth (get rid of cyclotron structure)

    CUT_WINTIME,

    CUT_MULTIWIN_NEXTDT, // plot t2 and t2 - t1 before

    // === UV analysis goes here ===

    CUT_Z_CONTAINED,

    CUT_PC7_HIT,
    CUT_PC7_COORDINATE,

    CUTS_LOOSE_PROTONS,

    // plot lastPlane distribution here
    // plot number of ranges
    // plot begin-end of ranges
    // plot begin-end of holes
    // plot missing planes
    // Plot TDC widths/particle ID distributions

    CUT_MIN_RANGE,
    CUT_RANGE_GAPS, // in the "stream" direction
    CUT_REXT,
    CUTS_TIGHT_PROTONS,
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

  double cutRextMax_;

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

  TH1 *hLastPlaneLoose_;

  MuCapUVAnalysis uvan_;
  HistPlaneRanges hRangeDIO_;
  HistDriftTime hdriftPCFiltered_;

  HistHitStream hhsZContained_;

  HistPlaneRanges hRangeAfterPC7Cuts_;


  HistTDCWidth hwidthPCTightDIO_;
  HistTDCWidth hwidthDCTightDIO_;
  HistTDCWidth hwidthPCLooseProtons_;
  HistTDCWidth hwidthDCLooseProtons_;
  HistTDCWidth hwidthPCTightProtons_;
  HistTDCWidth hwidthDCTightProtons_;

  HistTDCParticleClassifier hcdio_;
  HistTDCParticleClassifier hcLooseProtons_;
  HistTDCParticleClassifier hcTightProtons_;

  EventCutNumber analyze(const EventClass& evt,
                         const TimeWindowingResults& wres,
                         const ROOT::Math::XYPoint& muStopUV,
                         const std::vector<ClustersByPlane>& afterTrigGlobalClusters);
};

#endif/*MuCapStreamAnalysis_h*/
