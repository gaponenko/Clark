// Andrei Gaponenko, 2013

#ifndef MuCapProtonWindow_h
#define MuCapProtonWindow_h

#include "EventClass.h"
#include "WireCluster.h"
#include "TimeWindow.h"
#include "TDCHitWP.h"

#include "HistTDCWidth.h"
#include "HistProtonWindow.h"
#include "MuCapUVAnalysis.h"
#include "MuCapContainmentCheck.h"
#include "HistMuCapRTruth.h"
#include "HistMuCapTruth.h"

#include "Math/Point2D.h"
#include "TAxis.h"
class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;

//================================================================
class MuCapProtonWindow {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_UPSTREAM, "Upstream");
    ax->SetBinLabel(1+CUT_NOPC7, "NOPC7");
    ax->SetBinLabel(1+CUT_RANGE_GAPS, "RANGE_GAPS");
    ax->SetBinLabel(1+CUT_MAX_RANGE, "MAX_RANGE");

    ax->SetBinLabel(1+CUTS_LOOSE_PROTONS, "Loose protons");

    ax->SetBinLabel(1+CUT_MIN_RANGE, "Min range");
    ax->SetBinLabel(1+CUT_REXT, "Rext");

    ax->SetBinLabel(1+CUTS_TIGHT_PROTONS, "Tight protons");
  }

public:
  enum EventCutNumber {
    CUT_UPSTREAM,
    CUT_NOPC7,
    CUT_RANGE_GAPS,
    CUT_MAX_RANGE,

    CUTS_LOOSE_PROTONS,

    CUT_MIN_RANGE,
    CUT_REXT,

    CUTS_TIGHT_PROTONS,
    CUTS_END
  };

  void init(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile &conf);

  void process(const ROOT::Math::XYPoint& muStopUV,
               const TimeWindow& protonWindowPC,
               const TimeWindow& protonWindowDC,
               const TDCHitWPPtrCollection& unassignedDCHits,
               const EventClass& evt
               );

  MuCapProtonWindow()
    : doMCTruth_(false)
    , cutMaxPlane_()
    , cutRextMax_()
    , h_cuts_r()
    , h_cuts_p()
    , hStartPos_()
    , hStartOffset_()
    , hStartClusterSize_()
    , hLastPlane_()
    , hProtonTime_()
    , hCCRvsPlaneDIO_()
    , hCCRvsPlaneProtons_()
    , hLastPlaneVsMCPstart_()
  {}

private :
  bool doMCTruth_;
  int cutMaxPlane_;
  double cutRextMax_;

  std::string tightProtonsOutFileName_;
  std::ofstream tightProtonsOutFile_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH2 *hNumClusters_;

  TH2 *hStartPos_;
  TH2 *hStartOffset_; // from the muon stop position
  TH2 *hStartClusterSize_;

  TH1 *hLastPlane_;
  TH1 *hProtonTime_;

  TH2 *hCCRvsPlaneDIO_;
  TH2 *hCCRvsPlaneProtons_;

  TH2 *hLastPlaneVsMCPstart_;

  HistTDCWidth hwidthPCProtonWin_;
  HistTDCWidth hwidthDCProtonWin_;
  MuCapUVAnalysis uvan_;
  HistProtonWindow hpw_;

  MuCapContainmentCheck rcheckDIO_;
  MuCapContainmentCheck rcheckProtonCandidates_;
  HistMuCapRTruth hrtruth_;
  HistMuCapTruth htruthLoose_;
  HistMuCapTruth htruthMinRange_;
  HistMuCapTruth htruthTight_;

  EventCutNumber analyze(const ROOT::Math::XYPoint& muStopUV,
                         const TimeWindow& protonWindowPC,
                         const TimeWindow& protonWindowDC,
                         const TDCHitWPPtrCollection& unassignedDCHits,
                         const EventClass& evt
                         );
};

#endif/*MuCapProtonWindow_h*/
