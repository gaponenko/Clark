// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include "ModuleClass.h"

#include "TDCHitWP.h"
#include "TimeWindow.h"
#include "WireCluster.h"

#include "MuCapPACT.h"
#include "MuCapProtonWindow.h"
#include "HistTDCWidth.h"

#include "TAxis.h"
class TH1;
class TH2;


//================================================================
class MuCapture : public ModuleClass {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_NOPCWIN ,"NOPCs");

    ax->SetBinLabel(1+CUT_PCWIN_TRIGSEPPAST ,"PCTRIGSEP PAST");

    ax->SetBinLabel(1+CUT_PCWIN_NUMAFTERTRIG ,"Num after-trig windows");

    ax->SetBinLabel(1+CUT_PCWIN_TRIGSEPFUTURE ,"PCTRIGSEP FUTURE");

    ax->SetBinLabel(1+CUT_UNASSIGNEDDCHITS ,"Unassigned DC hits");

    ax->SetBinLabel(1+CUT_MU_RANGE_GAPS ,"Mu range gaps");

    ax->SetBinLabel(1+CUT_MU_PC_RANGE ,"Mu PC range");
    ax->SetBinLabel(1+CUT_MU_DC_RANGE ,"Mu DC range");

    ax->SetBinLabel(1+CUT_MU_UV_PC,"mu UV PC");
    ax->SetBinLabel(1+CUT_MU_UV_DC,"mu UV DC");

    ax->SetBinLabel(1+CUT_MUSTOP_SINGLECLUSTER,"Mu single cluster");
    ax->SetBinLabel(1+CUT_MUSTOP_UV,"Mu stop UV");

    ax->SetBinLabel(1+CUT_MUSTOP_PACT,"Mu stop PACT");

    ax->SetBinLabel(1+CUTS_ACCEPTED ,"Accepted");
  }

public :

  enum EventCutNumber {
    CUT_NOPCWIN,
    CUT_PCWIN_TRIGSEPPAST,
    CUT_PCWIN_NUMAFTERTRIG,
    CUT_PCWIN_TRIGSEPFUTURE,
    CUT_UNASSIGNEDDCHITS,
    CUT_MU_RANGE_GAPS,
    CUT_MU_PC_RANGE,
    CUT_MU_DC_RANGE,
    CUT_MU_UV_PC,
    CUT_MU_UV_DC,

    CUT_MUSTOP_SINGLECLUSTER,
    CUT_MUSTOP_UV,

    CUT_MUSTOP_PACT,

    CUTS_ACCEPTED,
    CUTS_END
  };

  virtual bool Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) ;
  virtual bool Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : doDefaultTWIST_(false)
    , cutMinTDCWidthPC_()
    , cutMinTDCWidthDC_()
    , winPCLength_()
    , winPCSeparation_()
    , winDCLength_()
    , winDCEarlyMargin_()
    , maxUnassignedDCHits_()
    , muUVPCCellMin_()
    , muUVPCCellMax_()
    , muUVDCCellMin_()
    , muUVDCCellMax_()
    , muStopRMax_()

    , h_cuts_r()
    , h_cuts_p()
    , hNumPCWin_()
    , hWinPCTimeAll_()
    , hWinPCTimeTrig_()
    , hWinPCTStartBeforeTrig_()
    , hWinPCTStartAfterTrig_()
    , hWinDCUnassignedEarly_()
    , hWinDCUnassignedLate_()
    , hWinDCUnassignedCount_()

    , hMuRangePCFirst_()
    , hMuRangePCLast_()
    , hMuRangeDCFirst_()
    , hMuRangeDCLast_()

    , hMuUVLimitsPCUp_()
    , hMuUVLimitsDC_()

    , hMuStopUVPos_()
    , hMuStopRadius_()
  {}

private :
  bool doDefaultTWIST_;

  double cutMinTDCWidthPC_;
  double cutMinTDCWidthDC_;

  double winPCLength_;
  double winPCSeparation_;

  double winDCLength_;
  double winDCEarlyMargin_;
  int maxUnassignedDCHits_;

  int muUVPCCellMin_;
  int muUVPCCellMax_;
  int muUVDCCellMin_;
  int muUVDCCellMax_;
  double muStopRMax_;

  MuCapPACT pactCut_;
  MuCapProtonWindow protonWindow_;

  TH1D *h_cuts_r;
  TH1D *h_cuts_p;

  TH1* hNumPCWin_;
  TH1* hWinPCTimeAll_;
  TH1* hWinPCTimeTrig_;
  TH1* hWinPCTStartBeforeTrig_;
  TH1* hWinPCTStartAfterTrig_;

  TH1 *hWinDCUnassignedEarly_;
  TH1 *hWinDCUnassignedLate_;
  TH1 *hWinDCUnassignedCount_;

  TH1 *hMuRangePCFirst_;
  TH1 *hMuRangePCLast_;
  TH1 *hMuRangeDCFirst_;
  TH1 *hMuRangeDCLast_;

  TH2 *hMuUVLimitsPCUp_;
  TH2 *hMuUVLimitsDC_;

  TH2 *hMuStopUVPos_;
  TH1 *hMuStopRadius_;

  HistTDCWidth hwidthPCall_;
  HistTDCWidth hwidthDCall_;

  HistTDCWidth hwidthPCMuWin_;
  HistTDCWidth hwidthDCMuWin_;

  EventCutNumber analyze(EventClass &E, HistogramFactory &H);

  TDCHitWPPtrCollection selectHits(const TDCHitWPCollection& hits, double minWidthCut);

  TimeWindowCollection constructTimeWindows(const TDCHitWPPtrCollection& timeSortedHits, double winLength);

  // Returns the index of the window with start time closest to t=0
  int findTriggerWindow(const TimeWindowCollection& windows);

  TimeWindowCollection assignDCHits(TDCHitWPPtrCollection* unassignedDCHits, const TDCHitWPPtrCollection& timeSortedDCHits, const TimeWindowCollection& winpcs);
};

#endif/*MuCapture_h*/
