// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include "ModuleClass.h"

#include "TDCHitWP.h"
#include "TimeWindow.h"

class TH1;
class TH2;

class MuCapture : public ModuleClass {
public :
  bool    Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) ;
  bool    Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : doDefaultTWIST_(false)
    , winPCLength_()
    , winPCSeparation_()
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

  TH1* hNumPCWin_;
  TH1* hWinPCTimeAll_;
  TH1* hWinPCTimeTrig_;
  TH1* hWinPCTStartBeforeTrig_;
  TH1* hWinPCTStartAfterTrig_;

  // Returns the index of the window with start time closest to t=0

  int findTriggerWindow(const TimeWindowCollection& windows);
};

#endif/*MuCapture_h*/
