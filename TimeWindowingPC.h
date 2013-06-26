// Andrei Gaponenko, 2013

#ifndef TimeWindowingPC_h
#define TimeWindowingPC_h

#include <string>

#include "TDCHitWP.h"
#include "TimeWindow.h"

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class TH1;
class TH2;

class TimeWindowingPC {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void assignPCHits(const TDCHitWPPtrCollection& hits, TimeWindowingResults *out);

  TimeWindowingPC()
    : winPCLength_()
    , winTrigMaxdt_()
    , hNumPCWin_()
    , hNumPCHitsPerWin_()
    , hWinPCTimeAll_()
    , hWinPCTimeTrig_()
    , hWinPCTStartBeforeTrig_()
    , hWinPCTStartAfterTrig_()
    , hWinTimeTypes_()
    , hWinMultTypes_()
  {}

private:
  double winPCLength_;
  double winTrigMaxdt_;
  TH1* hNumPCWin_;
  TH1* hNumPCHitsPerWin_;
  TH1* hWinPCTimeAll_;
  TH1* hWinPCTimeTrig_;
  TH1* hWinPCTStartBeforeTrig_;
  TH1* hWinPCTStartAfterTrig_;
  TH2* hWinTimeTypes_;
  TH2* hWinMultTypes_;

  // Returns the index of the trigger window or (-1u) if none is found
  unsigned findTriggerWindow(const TimeWindowCollection& windows);
};

#endif/*TimeWindowingPC_h*/
