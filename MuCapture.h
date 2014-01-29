// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include <string>
#include <fstream>

#include "ModuleClass.h"

#include "TDCHitWP.h"
#include "WireCluster.h"

#include "MuCapPACTCut.h"
#include "MuCapProtonWindow.h"
#include "MuCapStreamAnalysis.h"
#include "MuCapUVAnalysis.h"
#include "HistDriftTime.h"
#include "HistTDCWidth.h"
#include "HistOccupancy.h"
#include "HistAfterPulsing.h"
#include "HistXtalk.h"
#include "HistMuCapTruth.h"
#include "HistMuStopTruth.h"
#include "HistAccidentals.h"
#include "HistWinTime.h"
#include "TimeWindowingPC.h"
#include "TimeWindowingDC.h"
#include "EventList.h"
#include "TDCHitPreprocessing.h"
#include "RecoResMuCapTrk.h"

#include "RooUnfold/RooUnfoldResponse.h"

#include "TAxis.h"
class TH1;
class TH2;

//================================================================
class MuCapture : public ModuleClass {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_EVENT_NUMBER, "Event number");
    ax->SetBinLabel(1+CUT_NOPCHITS, "NoPCHits");
    ax->SetBinLabel(1+CUT_NOTRIGWIN, "NoTrigWin");
    ax->SetBinLabel(1+CUT_PCWIN_TRIGSEPPAST, "Pre-trig hits");
    ax->SetBinLabel(1+CUT_UNASSIGNEDDCHITS, "Unassigned DC hits");

    ax->SetBinLabel(1+CUT_MUON_FIRST_PLANE, "mu first plane");
    ax->SetBinLabel(1+CUT_MUON_RANGE_GAPS, "mu range gaps");
    ax->SetBinLabel(1+CUT_MUON_LAST_PLANE, "mu last plane");

    ax->SetBinLabel(1+CUT_MUSTOP_SINGLECLUSTER, "Mu single cluster");
    ax->SetBinLabel(1+CUT_MUSTOP_UV, "Mu stop UV");

    ax->SetBinLabel(1+CUT_MUSTOP_PACT, "Mu stop PACT");

    ax->SetBinLabel(1+CUTS_MUSTOP_ACCEPTED, "Accepted stop");
  }

public :

  enum EventCutNumber {
    CUT_EVENT_NUMBER,
    CUT_NOPCHITS,
    CUT_NOTRIGWIN,
    CUT_PCWIN_TRIGSEPPAST,
    CUT_UNASSIGNEDDCHITS,

    CUT_MUON_FIRST_PLANE,
    CUT_MUON_RANGE_GAPS,
    CUT_MUON_LAST_PLANE,

    CUT_MUSTOP_SINGLECLUSTER,
    CUT_MUSTOP_UV,

    CUT_MUSTOP_PACT,

    CUTS_MUSTOP_ACCEPTED,

    CUTS_END
  };

  virtual bool Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log);
  virtual bool Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : fillXtalkPC_(true)
    , fillXtalkDC_(false)
    , doMCTruth_(false)
    , doDefaultTWIST_(false)
    , winPCPreTrigSeparation_()
    , maxUnassignedDCHits_()
    , cutMuonFirstPlane_()
    , cutMuonRangeGapsEnabled_()
    , muStopRMax_()

    , pcHitProcessor_()
    , dcHitProcessor_()

    , h_cuts_r()
    , h_cuts_p()
    , hPCPreTrigSeparation_()
    , hWinDCUnassignedCount_()

    , hMuonFirstPlane_()
    , hMuonLastPlaneBeforeGaps_()
    , hMuonLastPlaneAfterGaps_()
    , hMuonRangeGaps_()
    , hMuonMissingPlanes_()

    , hMuStopUVCell_()
    , hMuStopUVPos_()
    , hMuStopRadius_()
  {}

private :
  bool fillXtalkPC_;
  bool fillXtalkDC_;
  bool doMCTruth_;
  bool doDefaultTWIST_;
  EventList inputNumberList_;

  std::string uvOutFileName_;
  std::ofstream uvOutFile_;

  double winPCPreTrigSeparation_;

  int maxUnassignedDCHits_;

  int cutMuonFirstPlane_;
  bool cutMuonRangeGapsEnabled_;

  double muStopRMax_;

  TDCHitPreprocessing::IProcessor *pcHitProcessor_;
  TDCHitPreprocessing::IProcessor *dcHitProcessor_;

  RecoResMuCapTrk anDnLateRes_;
  RooUnfoldResponse anDnLateResponse_;
  TH1D* hTruthMomentum_;
  TH1D* hTruthMomentumReco_;
  TH2D* hMeasVsTruthMomentum_;
  TH1D* hTruthMomentumNotReco_;
  TH1D* hMeasuredMomentum_;

  TimeWindowingPC pcWindowing_;
  TimeWindowingDC dcWindowing_;
  MuCapPACTCut pactCut_;
  MuCapProtonWindow protonWindow_;
  MuCapStreamAnalysis anUpLate_;
  MuCapStreamAnalysis anDnLate_;
  MuCapStreamAnalysis anDnEarly_;

  TH1D *h_cuts_r;
  TH1D *h_cuts_p;

  TH1 *hPCPreTrigSeparation_;
  TH1 *hWinDCUnassignedCount_;

  TH1 *hMuonFirstPlane_;
  TH1 *hMuonLastPlaneBeforeGaps_;
  TH1 *hMuonLastPlaneAfterGaps_;
  TH2 *hMuonRangeGaps_;
  TH1 *hMuonMissingPlanes_;

  TH2 *hMuStopUVCell_;
  TH2 *hMuStopUVPos_;
  TH1 *hMuStopRadius_;

  HistTDCWidth hwidthPCall_;
  HistTDCWidth hwidthPCfiltered_;

  HistTDCWidth hwidthDCall_;
  HistTDCWidth hwidthDCfiltered_;

  HistAfterPulsing hAfterPulsingPCAll_;
  HistAfterPulsing hAfterPulsingPCFiltered_;
  HistAfterPulsing hAfterPulsingDCAll_;
  HistAfterPulsing hAfterPulsingDCFiltered_;

  HistXtalk hXtalkSameWirePC_;
  HistXtalk hXtalk1PC_;
  HistXtalk hXtalkPlanePC_;
  HistXtalk hXtalkSameWireDC_;
  HistXtalk hXtalk1DC_;
  HistXtalk hXtalkPlaneDC_;

  HistOccupancy hOccupancyPCAll_;
  HistOccupancy hOccupancyDCAll_;
  HistAccidentals haccidentalsTrig_;
  HistAccidentals haccidentalsStop_;

  HistWinTime winTimeBeforeNoTrigWin_;
  HistWinTime winTimeMuStop_;

  MuCapUVAnalysis dioUp_;
  MuCapUVAnalysis dioDn_;

  HistDriftTime hdriftPCAll_;
  HistDriftTime hdriftPCFiltered_;
  HistDriftTime hdriftDCAll_;
  HistDriftTime hdriftDCFiltered_;

  HistMuStopTruth hmuStopTruthAll_;
  HistMuStopTruth hmuStopTruthAfterGaps_;

  HistMuCapTruth hTruthAll_;
  HistMuCapTruth hTruthAfterPreTrigHits_;
  HistMuCapTruth hTruthAfterUnassignedDCHits_;
  HistMuCapTruth hTruthAfterMuRangeGaps_;
  HistMuCapTruth hTruthAfterMuLastPlane_;
  HistMuCapTruth hTruthAfterMuSingleCluster_;
  HistMuCapTruth hTruthAfterMuStopUV_;
  HistMuCapTruth hTruthMuStop_;

  EventCutNumber analyze(EventClass &E, HistogramFactory &H);

  static TDCHitPreprocessing::IProcessor*
  makeTDCHitPreprocessor(WirePlane::DetType d, HistogramFactory& hf, const DetectorGeo& geom, ConfigFile& conf);
};

#endif/*MuCapture_h*/
