// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include "ModuleClass.h"

#include "TDCHitWP.h"
#include "TimeWindow.h"

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

    ax->SetBinLabel(1+CUT_MUZ ,"mu Z");
    ax->SetBinLabel(1+CUT_MU_UV,"mu UV");
    ax->SetBinLabel(1+CUTS_ACCEPTED ,"Accepted");
  }

public :

  enum EventCutNumber {
    CUT_NOPCWIN,
    CUT_PCWIN_TRIGSEPPAST,
    CUT_PCWIN_NUMAFTERTRIG,
    CUT_PCWIN_TRIGSEPFUTURE,
    CUT_MUZ,
    CUT_MU_UV,

    CUTS_ACCEPTED,
    CUTS_END
  };

  virtual bool Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) ;
  virtual bool Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : doDefaultTWIST_(false)
    , winPCLength_()
    , winPCSeparation_()
    , h_cuts_r()
    , h_cuts_p()
    , hNumPCWin_()
    , hWinPCTimeAll_()
    , hWinPCTimeTrig_()
    , hWinPCTStartBeforeTrig_()
    , hWinPCTStartAfterTrig_()
  {}

private :
  bool doDefaultTWIST_;
  double winPCLength_;
  double winPCSeparation_;

  TH1D *h_cuts_r;
  TH1D *h_cuts_p;

  TH1* hNumPCWin_;
  TH1* hWinPCTimeAll_;
  TH1* hWinPCTimeTrig_;
  TH1* hWinPCTStartBeforeTrig_;
  TH1* hWinPCTStartAfterTrig_;

  EventCutNumber analyze(EventClass &E, HistogramFactory &H);

  // Returns the index of the window with start time closest to t=0
  int findTriggerWindow(const TimeWindowCollection& windows);
};

#endif/*MuCapture_h*/
