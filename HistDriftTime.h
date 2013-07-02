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
            const ConfigFile &conf);

  void fill(const EventClass& evt, int idio, const TDCHitWPPtrCollection& hits);

private :
  std::vector<TH1*> hdriftTime_;
};

#endif/*HistDriftTime_h*/
