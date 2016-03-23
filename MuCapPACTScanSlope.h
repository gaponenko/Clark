// Andrei Gaponenko, 2013

#ifndef MuCapPACTScanSlope_h
#define MuCapPACTScanSlope_h

#include <string>
#include <vector>

#include "HistMuCapMCTgtStops.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class WireCluster;
class EventClass;
class DetectorGeo;

//================================================================
class MuCapPACTScanSlope {
public:
  void fill(const EventClass& evt, const WireCluster& pc5cluster, const WireCluster& pc6cluster);

  MuCapPACTScanSlope(HistogramFactory& hf,
                     const std::string& hdir,
                     const DetectorGeo& geom,
                     const ConfigFile& conf,
                     const std::string& suffix);

private :
  double slopea_;
  std::vector<double> slopeb_;
  std::vector<TH2*> hia_vs_ib_scan_;
};

#endif/*MuCapPACTScanSlope_h*/
