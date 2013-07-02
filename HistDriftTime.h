// Andrei Gaponenko, 2013

#ifndef HistDriftTime_h
#define HistDriftTime_h

#include <string>
#include <vector>

#include "TDCHitWP.h"

class TH1;
class TH2;
class TProfile2D;
class TEfficiency;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistDriftTime {
public:
  void init(const std::string& hdir,
            HistogramFactory& hf,
            unsigned maxPlaneNumber,
            double driftTimeHistLimit,
            double cutEffTrackHitDt,
            const ConfigFile &conf);

  void fill(const EventClass& evt, int idio, const TDCHitWPPtrCollection& hits);

  HistDriftTime() : cutEffTrackHitDt_(), hexpectedPlaneHits_(), hobservedPlaneHits_() {}

private :
  std::vector<TH1*> hdriftTime_;

  double cutEffTrackHitDt_;
  TH1 *hexpectedPlaneHits_;
  TH1 *hobservedPlaneHits_;

  TH2 *hexpectedPlaneHitsVsTime_;
  TH2 *hobservedPlaneHitsVsTime_;
};

#endif/*HistDriftTime_h*/
