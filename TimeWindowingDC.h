// Andrei Gaponenko, 2013

#ifndef TimeWindowingDC_h
#define TimeWindowingDC_h

#include <string>

#include "TDCHitWP.h"
#include "TimeWindow.h"
#include "HistOccupancy.h"

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class TH1;
class TH2;

class TimeWindowingDC {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void assignDCHits(const TDCHitWPPtrCollection& dcHits, TimeWindowingResults *inout);

  TimeWindowingDC()
    : winDCLength_()
    , winDCEarlyMargin_()
    , winDCDoHistos_()
    , hWinDCUnassignedEarly_()
    , hWinDCUnassignedLate_()
  {}

private:
  double winDCLength_;
  double winDCEarlyMargin_;
  bool   winDCDoHistos_;

  HistOccupancy hWinDCMuMixStreamMap_;
  HistOccupancy hWinDCProtonMixStreamMap_;

  void fillDiagnostics(const TimeWindowingResults& wres);
  HistOccupancy hOccupancyDCUnassigned_;
  TH1 *hWinDCUnassignedEarly_;
  TH1 *hWinDCUnassignedLate_;
};

#endif/*TimeWindowingDC_h*/
