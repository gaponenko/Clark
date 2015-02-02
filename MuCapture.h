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
#include "MuCapUVAnalysis.h"
#include "MuCapTrkAnalysisHF.h"
#include "HistMuCapAnalysisChannels.h"
#include "HistDriftTime.h"
#include "HistTDCWidth.h"
#include "HistTDCBCSWidth.h"
#include "HistOccupancy.h"
#include "HistAfterPulsing.h"
#include "HistXtalk.h"
#include "HistXT2.h"
#include "HistMuStopTruth.h"
#include "HistAccidentals.h"
#include "HistWinDCUnassigned.h"
#include "HistProtonPID.h"
#include "TimeWindowingPC.h"
#include "TimeWindowingDC.h"
#include "EventList.h"
#include "TDCHitPreprocessing.h"
#include "RecoResMuCapTrk.h"
#include "MuCapTrkContainmentCut.h"

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

    ax->SetBinLabel(1+CUT_PCWIN_TRIGSEPPAST, "Pre-trig win");

    ax->SetBinLabel(1+CUT_MUON_FIRST_PLANE, "mu first plane");
    ax->SetBinLabel(1+CUT_MUON_LAST_PLANE, "mu last plane");

    ax->SetBinLabel(1+CUT_MUSTOP_SINGLECLUSTER, "Mu single cluster");
    ax->SetBinLabel(1+CUT_MUSTOP_UV, "Mu stop UV");

    ax->SetBinLabel(1+CUT_MUSTOP_PACT, "Mu stop PACT");

    //ax->SetBinLabel(1+CUTS_MUSTOP_ACCEPTED, "Accepted stop");

    ax->SetBinLabel(1+CUT_DOWNSTREAM_PCHITS, "Dn PC");

    ax->SetBinLabel(1+CUT_BEAM_VETO, "Beam veto");

    ax->SetBinLabel(1+CUT_WIN_TIME, "Win time");

    ax->SetBinLabel(1+CUT_MULTIWIN_NEXTDT, "Multiwin time");

    ax->SetBinLabel(1+CUTS_DOWNSTREAM_ACCEPTED, "Dn candidate");
  }

public :

  enum EventCutNumber {
    CUT_EVENT_NUMBER,
    CUT_NOPCHITS,
    CUT_NOTRIGWIN,
    CUT_PCWIN_TRIGSEPPAST,
    CUT_MUON_FIRST_PLANE,
    CUT_MUON_LAST_PLANE,

    CUT_MUSTOP_SINGLECLUSTER,
    CUT_MUSTOP_UV,

    CUT_MUSTOP_PACT,
    // CUTS_MUSTOP_ACCEPTED,

    CUT_DOWNSTREAM_PCHITS,

    CUT_BEAM_VETO,

    // cut on the same var for DIO and protons
    // avoid norm systematic due to time resolution
    CUT_WIN_TIME,

    CUT_MULTIWIN_NEXTDT,

    CUTS_DOWNSTREAM_ACCEPTED,

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
    , cutMuonFirstPlane_()
    , muStopRMax_()
    , cutBeamVetoMaxPCplanes_()
    , cutMultiwinNextdt_()

    , pcXTalkProcessor_()
    , dcXTalkProcessor_()

    , pcHitProcessor_()
    , dcHitProcessor_()

    , h_cuts_r()
    , h_cuts_p()
    , hPCPreTrigSeparation_()

    , hMuonFirstPlane_()
    , hMuonLastPlaneAfterGaps_()
    , hStoppedMuonRangeGaps_()
    , hStoppedMuonMissingPlanes_()

    , hMuStopUVCell_()
    , hMuStopUVPos_()
    , hMuStopRadius_()

    , hBeamVetoNumHitPlanes_()
    , hHitPCsAterBeamVeto_()

    , hWindowTimeBefore_()
    , hWindowTimeAfter_()
    , hNumAfterTrigWindows_()
    , hWindow2Time_()
    , hWindow2dt_()

    , hPosNegMom_()
    , hPosNegCosth_()
    , hVetoedPCosth_()
    , hVetoingPCosth_()

  {}

private :
  bool fillXtalkPC_;
  bool fillXtalkDC_;
  bool doMCTruth_;
  bool doDefaultTWIST_;
  EventList inputNumberList_;

  std::string uvOutFileName_;
  std::ofstream uvOutFile_;
  std::string commonSkimOutFileName_;
  std::ofstream commonSkimOutFile_;

  double winPCPreTrigSeparation_;

  int cutMuonFirstPlane_;

  double muStopRMax_;

  int cutBeamVetoMaxPCplanes_;

  double cutWinTimeMin_;
  double cutWinTimeMax_;
  double cutMultiwinNextdt_; // min t2-t1

  TDCHitPreprocessing::IProcessor *pcXTalkProcessor_;
  TDCHitPreprocessing::IProcessor *dcXTalkProcessor_;

  TDCHitPreprocessing::IProcessor *pcHitProcessor_;
  TDCHitPreprocessing::IProcessor *dcHitProcessor_;

  HistMuCapAnalysisChannels channels_;
  RecoResMuCapTrk anDnLateRes_;
  RooUnfoldResponse anDnLateResponse_;
  TH1D* hTruthMomentum_;
  TH1D* hTruthMomentumReco_;
  TH2D* hMeasVsTruthMomentum_;
  TH1D* hTruthMomentumNotReco_;
  TH1D* hMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponseWithPID_;
  TH2D* hWithPIDTruthMomentum_;
  TH2D* hWithPIDTruthMomentumReco_;
  TH2D* hWithPIDMeasVsTruthMomentumTruProtons_;
  TH2D* hWithPIDMeasVsTruthMomentumTruDeuterons_;
  TH2D* hWithPIDTruthMomentumNotReco_;
  TH2D* hWithPIDMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponseWithPIDAllTrks_;
  TH2D* hWithPIDAllTrksTruthMomentum_;
  TH2D* hWithPIDAllTrksTruthMomentumReco_;
  TH2D* hWithPIDAllTrksMeasVsTruthMomentumTruProtons_;
  TH2D* hWithPIDAllTrksMeasVsTruthMomentumTruDeuterons_;
  TH2D* hWithPIDAllTrksTruthMomentumNotReco_;
  TH2D* hWithPIDAllTrksMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponseContained_;
  TH1D* hContainedTruthMomentum_;
  TH1D* hContainedTruthMomentumReco_;
  TH2D* hContainedMeasVsTruthMomentum_;
  TH1D* hContainedTruthMomentumNotReco_;
  TH1D* hContainedMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnRngCutPln_;
  TH1D* hPlnRngCutPlnTruthMomentum_;
  TH1D* hPlnRngCutPlnTruthMomentumReco_;
  TH2D* hPlnRngCutPlnMeasVsTruthMomentum_;
  TH1D* hPlnRngCutPlnTruthMomentumNotReco_;
  TH1D* hPlnRngCutPlnMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone1_;
  TH1D* hPlnVsPCutZone1TruthMomentum_;
  TH1D* hPlnVsPCutZone1TruthMomentumReco_;
  TH2D* hPlnVsPCutZone1MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone1TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone1MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone2_;
  TH1D* hPlnVsPCutZone2TruthMomentum_;
  TH1D* hPlnVsPCutZone2TruthMomentumReco_;
  TH2D* hPlnVsPCutZone2MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone2TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone2MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone3_;
  TH1D* hPlnVsPCutZone3TruthMomentum_;
  TH1D* hPlnVsPCutZone3TruthMomentumReco_;
  TH2D* hPlnVsPCutZone3MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone3TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone3MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone4_;
  TH1D* hPlnVsPCutZone4TruthMomentum_;
  TH1D* hPlnVsPCutZone4TruthMomentumReco_;
  TH2D* hPlnVsPCutZone4MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone4TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone4MeasuredMomentum_;

  TimeWindowingPC pcWindowing_;
  TimeWindowingDC dcWindowing_;
  MuCapPACTCut pactCut_;
  MuCapProtonWindow protonWindow_;

  TH1D *h_cuts_r;
  TH1D *h_cuts_p;

  TH1 *hPCPreTrigSeparation_;

  TH1 *hMuonFirstPlane_;
  TH1 *hMuonLastPlaneAfterGaps_;

  TH2 *hStoppedMuonRangeGaps_;
  TH1 *hStoppedMuonMissingPlanes_;

  TH2 *hMuStopUVCell_;
  TH2 *hMuStopUVPos_;
  TH1 *hMuStopRadius_;

  TH1 *hBeamVetoNumHitPlanes_;
  TH1 *hHitPCsAterBeamVeto_; // occupancy after the cut

  TH1 *hWindowTimeBefore_;
  TH1 *hWindowTimeAfter_;

  TH1 *hNumAfterTrigWindows_;
  TH1 *hWindow2Time_;
  TH1 *hWindow2dt_;

  TH2 *hPosNegMom_; // events with both positive and negative tracks accepted
  TH2 *hPosNegCosth_; // events with both positive and negative tracks accepted
  TH2 *hVetoedPCosth_;
  TH2 *hVetoingPCosth_;

  HistTDCWidth hwidthPCall_;
  HistTDCWidth hwidthPCfiltered_;

  HistTDCWidth hwidthDCall_;
  HistTDCWidth hwidthDCfiltered_;

  HistTDCBCSWidth hwidthMuHits_;
  HistTDCBCSWidth hwidthDIOHits_;

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
  HistXT2   hXT2PlanePC_;
  HistXT2   hXT2PlaneDC_;

  HistOccupancy hOccupancyPCAll_;
  HistOccupancy hOccupancyDCAll_;
  HistAccidentals haccidentalsTrig_;
  HistAccidentals haccidentalsStop_;

  HistWinDCUnassigned winDCUnassignedAfterWindowing_;
  HistWinDCUnassigned winDCUnassignedMuStop_;
  HistWinDCUnassigned winDCUnassignedDnDecay_;

  MuCapUVAnalysis dioUp_;
  MuCapUVAnalysis dioDn_;

  MuCapTrkAnalysisHF dnPosTracks_;
  MuCapTrkAnalysisHF dnDIOVetoTracks_;
  MuCapTrkAnalysisHF dnDIONormTracks_;

  MuCapTrkContainmentCut dnPosTrkContainment_;

  HistProtonPID hProtonPID_;
  HistProtonPID hContainedProtonPID_;

  HistDriftTime hdriftPCAll_;
  HistDriftTime hdriftPCFiltered_;
  HistDriftTime hdriftDCAll_;
  HistDriftTime hdriftDCFiltered_;

  HistMuStopTruth hmuStopTruthAll_;
  HistMuStopTruth hmuStopTruthAfterGaps_;

  HistMuCapTruth hTruthAll_;
  HistMuCapTruth hTruthMuStop_;
  HistMuCapTruth hTruthDnCandidate_;

  EventCutNumber analyze(EventClass &E, HistogramFactory &H);

  static TDCHitPreprocessing::IProcessor*
  makeXTalkPreprocessor(const std::string& topdir, WirePlane::DetType d, HistogramFactory& hf, const DetectorGeo& geom, ConfigFile& conf);

  static TDCHitPreprocessing::IProcessor*
  makeTDCHitPreprocessor(const std::string& topdir, WirePlane::DetType d, HistogramFactory& hf, const DetectorGeo& geom, ConfigFile& conf);
};

#endif/*MuCapture_h*/
