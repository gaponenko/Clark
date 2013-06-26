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
#include "HistTDCWidth.h"
#include "HistOccupancy.h"
#include "HistAfterPulsing.h"
#include "HistMuCapTruth.h"
#include "HistAccidentals.h"
#include "HistWinTime.h"
#include "TimeWindowingPC.h"
#include "TimeWindowingDC.h"
#include "EventList.h"

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
    ax->SetBinLabel(1+CUT_TRIGPCWIN_TYPE, "TrigPCWinType");
    ax->SetBinLabel(1+CUT_TRIGPCWIN_GAPS, "TrigPCWinGaps");
    ax->SetBinLabel(1+CUT_TRIGPCWIN_END_PLANE, "TrigPCWinEndPlane");
    ax->SetBinLabel(1+CUT_TRIGPCWIN_START_PLANE, "TrigPCWinStartPlane");

    ax->SetBinLabel(1+CUT_TRIGDCWIN_TYPE, "TrigDCWinType");
    ax->SetBinLabel(1+CUT_UNASSIGNEDDCHITS, "Unassigned DC hits");

    ax->SetBinLabel(1+CUT_MU_RANGE_GAPS, "Mu range gaps");

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

    CUT_TRIGPCWIN_TYPE,
    CUT_TRIGPCWIN_END_PLANE,
    CUT_TRIGPCWIN_START_PLANE,
    CUT_TRIGPCWIN_GAPS,

    CUT_TRIGDCWIN_TYPE,
    CUT_UNASSIGNEDDCHITS,

    CUT_MU_RANGE_GAPS,

    CUT_MUSTOP_SINGLECLUSTER,
    CUT_MUSTOP_UV,

    CUT_MUSTOP_PACT,

    CUTS_MUSTOP_ACCEPTED,

    CUTS_END
  };

  virtual bool Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) ;
  virtual bool Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : doMCTruth_(false)
    , doDefaultTWIST_(false)
    , cutMinTDCWidthPC_()
    , cutMinTDCWidthDC_()
    , winPCPreTrigSeparation_()
    , maxUnassignedDCHits_()
    , muStopRMax_()

    , h_cuts_r()
    , h_cuts_p()
    , hPCPreTrigSeparation_()
    , hTrigPCWinStartPlane_()
    , hTrigPCWinGaps_()

    , hWinDCUnassignedCount_()

    , hMuStopUVCell_()
    , hMuStopUVPos_()
    , hMuStopRadius_()
  {}

private :
  bool doMCTruth_;
  bool doDefaultTWIST_;
  EventList inputNumberList_;

  std::string muStopOutFileName_;
  std::ofstream muStopOutFile_;

  double cutMinTDCWidthPC_;
  double cutMinTDCWidthDC_;

  double winPCPreTrigSeparation_;

  bool cutTrigPCWinGapsEnabled_;
  int cutTrigPCWinStartPlane_;

  int maxUnassignedDCHits_;

  double muStopRMax_;

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
  TH1 *hTrigPCWinStartPlane_;
  TH2 *hTrigPCWinGaps_;

  TH1 *hWinDCUnassignedCount_;

  TH2 *hMuStopUVCell_;
  TH2 *hMuStopUVPos_;
  TH1 *hMuStopRadius_;

  HistTDCWidth hwidthPCall_;
  HistTDCWidth hwidthDCall_;

  HistAfterPulsing hAfterPulsingPCAll_;
  HistAfterPulsing hAfterPulsingPCFiltered_;

  HistOccupancy hOccupancyPCAll_;
  HistOccupancy hOccupancyDCAll_;
  HistAccidentals haccidentals_;

  HistWinTime winTimeBeforeNoTrigWin_;
  HistWinTime winTimeBeforeTrigPCWinType_;
  HistWinTime winTimeBeforeTrigPCWinGaps_;

  HistWinTime winTimeBeforeTrigPCWinEndPlane_;
  HistWinTime winTimeBeforeTrigPCWinStartPlane_;

  HistWinTime winTimeBeforeTrigDCWinType_;
  HistWinTime winTimeMuStop_;

  HistMuCapTruth hTruthAll_;

  EventCutNumber analyze(EventClass &E, HistogramFactory &H);
  TDCHitWPPtrCollection selectHits(const TDCHitWPCollection& hits, double minWidthCut);
};

#endif/*MuCapture_h*/
