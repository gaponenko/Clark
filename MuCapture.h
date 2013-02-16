// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include "ModuleClass.h"

class TH1;
class TH2;

class MuCapture : public ModuleClass {
public :
  bool    Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) ;
  bool    Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : doDefaultTWIST_(false)
    , windtpc_()
    , hNumPCWin_()
    , hWinPCTimeAll_()
    , hWinPCTimeTrig_()
  {}

private :
  bool doDefaultTWIST_;
  double windtpc_;

  TH1* hNumPCWin_;
  TH1* hWinPCTimeAll_;
  TH1* hWinPCTimeTrig_;
};

#endif/*MuCapture_h*/
