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
#include "HistTDCParticleClassifier.h"

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
    ax->SetBinLabel(1+CUT_TRIGSEP, "dt trig");
    ax->SetBinLabel(1+CUT_STREAM, "stream");
    ax->SetBinLabel(1+CUT_NOPC7, "PC7");
    ax->SetBinLabel(1+CUT_RANGE_GAPS, "RANGE_GAPS");
    ax->SetBinLabel(1+CUT_MAX_RANGE, "MAX_RANGE");

    ax->SetBinLabel(1+CUTS_LOOSE_PROTONS, "Loose protons");

    ax->SetBinLabel(1+CUT_NOPC8, "PC8");

    ax->SetBinLabel(1+CUT_MIN_RANGE, "Min range");
    ax->SetBinLabel(1+CUT_REXT, "Rext");

    ax->SetBinLabel(1+CUTS_TIGHT_PROTONS, "Tight protons");
  }

public:
  enum EventCutNumber {
    CUT_STREAM,
    CUT_TRIGSEP,
    CUT_NOPC7,
    CUT_RANGE_GAPS,
    CUT_MAX_RANGE,

    CUTS_LOOSE_PROTONS,

    CUT_NOPC8,
    CUT_MIN_RANGE,
    CUT_REXT,

    CUTS_TIGHT_PROTONS,
    CUTS_END
  };

  void init(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile &conf,
            TimeWindow::StreamType cutWinStream, double cutAfterTrigTimeSep);

  void process(const ROOT::Math::XYPoint& muStopUV,
               const TimeWindowingResults& wres,
               const EventClass& evt
               );

  MuCapProtonWindow()
    : doMCTruth_(false)
    , cutStream_(TimeWindow::DOWNSTREAM)
    , cutAfterTrigTimeSep_()
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
  TimeWindow::StreamType cutStream_;
  double cutAfterTrigTimeSep_;
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

  HistTDCWidth hwidthPCTightProtons_;
  HistTDCWidth hwidthDCTightProtons_;
  HistTDCWidth hwidthPCTightDIO_;
  HistTDCWidth hwidthDCTightDIO_;

  HistTDCParticleClassifier hcLooseProtons_;
  HistTDCParticleClassifier hcTightProtons_;
  HistTDCParticleClassifier hcdio_;

  MuCapUVAnalysis uvan_;
  HistProtonWindow hpw_;

  MuCapContainmentCheck rcheckDIO_;
  MuCapContainmentCheck rcheckProtonCandidates_;
  HistMuCapRTruth hrtruth_;
  HistMuCapTruth htruthLoose_;
  HistMuCapTruth htruthPC8_;
  HistMuCapTruth htruthMinRange_;
  HistMuCapTruth htruthTight_;

  EventCutNumber analyze(const ROOT::Math::XYPoint& muStopUV,
                         const TimeWindowingResults& wres,
                         const EventClass& evt
                         );
};

#endif/*MuCapProtonWindow_h*/
