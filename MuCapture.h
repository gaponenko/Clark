// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include "ModuleClass.h"

#include "TDCHitWP.h"
#include "WireCluster.h"

#include "MuCapPACT.h"
#include "MuCapProtonWindow.h"
#include "HistTDCWidth.h"
#include "HistOccupancy.h"
#include "HistMuCapTruth.h"
#include "TimeWindowingPC.h"
#include "TimeWindowingDC.h"

#include "TAxis.h"
class TH1;
class TH2;

//================================================================
class MuCapture : public ModuleClass {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_NOPCHITS, "NoPCHits");
    ax->SetBinLabel(1+CUT_NOTRIGWIN, "NoTrigWin");
    ax->SetBinLabel(1+CUT_PCWIN_TRIGSEPPAST, "Pre-trig hits");
    ax->SetBinLabel(1+CUT_TRIGPCWIN_TYPE, "TrigPCWinType");
    ax->SetBinLabel(1+CUT_TRIGPCWIN_GAPS, "TrigPCWinGaps");
    ax->SetBinLabel(1+CUT_TRIGPCWIN_RANGE, "TrigPCWinRange");

    ax->SetBinLabel(1+CUT_TRIGDCWIN_TYPE, "TrigDCWinType");
    ax->SetBinLabel(1+CUT_UNASSIGNEDDCHITS, "Unassigned DC hits");

    ax->SetBinLabel(1+CUT_MU_RANGE_GAPS, "Mu range gaps");

    ax->SetBinLabel(1+CUT_MUSTOP_SINGLECLUSTER, "Mu single cluster");
    ax->SetBinLabel(1+CUT_MUSTOP_UV, "Mu stop UV");

    ax->SetBinLabel(1+CUT_MUSTOP_PACT, "Mu stop PACT");

    ax->SetBinLabel(1+CUT_MU_UV_PC, "mu UV PC");
    ax->SetBinLabel(1+CUT_MU_UV_DC, "mu UV DC");

    ax->SetBinLabel(1+CUTS_MUSTOP_ACCEPTED, "mu stop accepted");

    ax->SetBinLabel(1+CUT_WIN_NUMAFTERTRIG, "Num after-trig windows");

    ax->SetBinLabel(1+CUT_PROTONWIN_TYPE, "Proton win type");

    ax->SetBinLabel(1+CUTS_ACCEPTED, "mu int accepted");
  }

public :

  enum EventCutNumber {
    CUT_NOPCHITS,
    CUT_NOTRIGWIN,
    CUT_PCWIN_TRIGSEPPAST,

    CUT_TRIGPCWIN_TYPE,
    CUT_TRIGPCWIN_GAPS,
    CUT_TRIGPCWIN_RANGE,

    CUT_TRIGDCWIN_TYPE,
    CUT_UNASSIGNEDDCHITS,

    CUT_MU_RANGE_GAPS,

    CUT_MUSTOP_SINGLECLUSTER,
    CUT_MUSTOP_UV,

    CUT_MUSTOP_PACT,

    CUT_MU_UV_PC,
    CUT_MU_UV_DC,

    CUTS_MUSTOP_ACCEPTED,

    CUT_WIN_NUMAFTERTRIG,
    CUT_PROTONWIN_TYPE,

    CUTS_ACCEPTED,
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
    , muUVPCCellMin_()
    , muUVPCCellMax_()
    , muUVDCCellMin_()
    , muUVDCCellMax_()
    , muStopRMax_()

    , h_cuts_r()
    , h_cuts_p()
    , hPCPreTrigSeparation_()
    , hNumAfterTrigWindows_()

    , hWinDCUnassignedCount_()

    , hMuUVLimitsPCUp_()
    , hMuUVLimitsDC_()

    , hMuStopUVCell_()
    , hMuStopUVPos_()
    , hMuStopRadius_()
  {}

private :
  bool doMCTruth_;
  bool doDefaultTWIST_;

  double cutMinTDCWidthPC_;
  double cutMinTDCWidthDC_;

  double winPCPreTrigSeparation_;

  int maxUnassignedDCHits_;

  int muUVPCCellMin_;
  int muUVPCCellMax_;
  int muUVDCCellMin_;
  int muUVDCCellMax_;
  double muStopRMax_;

  TimeWindowingPC pcWindowing_;
  TimeWindowingDC dcWindowing_;
  MuCapPACT pactCut_;
  MuCapProtonWindow protonWindow_;

  TH1D *h_cuts_r;
  TH1D *h_cuts_p;

  TH1 *hPCPreTrigSeparation_;
  TH1 *hNumAfterTrigWindows_;

  TH1 *hWinDCUnassignedCount_;

  TH2 *hMuUVLimitsPCUp_;
  TH2 *hMuUVLimitsDC_;

  TH2 *hMuStopUVCell_;
  TH2 *hMuStopUVPos_;
  TH1 *hMuStopRadius_;

  HistTDCWidth hwidthPCall_;
  HistTDCWidth hwidthDCall_;

  HistOccupancy hOccupancyPCAll_;
  HistOccupancy hOccupancyDCAll_;
  HistMuCapTruth hTruthAll_;

  EventCutNumber analyze(EventClass &E, HistogramFactory &H);
  TDCHitWPPtrCollection selectHits(const TDCHitWPCollection& hits, double minWidthCut);
};

#endif/*MuCapture_h*/
