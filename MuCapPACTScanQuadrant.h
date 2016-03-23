// Andrei Gaponenko, 2013

#ifndef MuCapPACTScanQuadrant_h
#define MuCapPACTScanQuadrant_h

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
class MuCapPACTScanQuadrant {
public:
  void fill(const EventClass& evt, const WireCluster& pc5cluster, const WireCluster& pc6cluster);

  MuCapPACTScanQuadrant(HistogramFactory& hf,
                        const std::string& hdir,
                        const DetectorGeo& geom,
                        const ConfigFile& conf,
                        const std::string& suffix);

private :
  bool doMCTruth_;

  double slopea_;
  double intercepta_;

  double slopeb_;
  double interceptb_min_;
  double interceptb_max_;
  double interceptb_npoints_;
  std::vector<double> interceptb_;

  TH2 *hpc6vs5widthAll_;
  std::vector<TH2*> hpc6vs5widthQ1_scanb_;
  std::vector<HistMuCapMCTgtStops*> hMcMuStops_scanb_;
};

#endif/*MuCapPACTScanQuadrant_h*/
