// Andrei Gaponenko, 2013

#ifndef MuCapture_h
#define MuCapture_h

#include "ModuleClass.h"

class TH1;
class TH2;

class MuCapture : public ModuleClass {
public :
  bool    Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
  bool    Process(EventClass &E, HistogramFactory &H);

  MuCapture()
    : Log(0)
    , doDefaultTWIST(false)
    , h_dt_pc_(0)
    , h_dt_dc_(0)
    , h_dt_any_(0)
    , h_dt_any2_(0)
    , h_m12_t0_(0)
  {}

private :
  log4cpp::Category *Log;
  string Name;
  bool doDefaultTWIST;

  TH1* h_dt_pc_;
  TH1* h_dt_dc_;
  TH1* h_dt_any_;
  TH2* h_dt_any2_;
  TH1* h_m12_t0_;
};

#endif/*MuCapture_h*/
