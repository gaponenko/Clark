// Andrei Gaponenko, 2013

#ifndef MuCapPACTScan_h
#define MuCapPACTScan_h

#include <vector>
#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class WireCluster;
class EventClass;
class DetectorGeo;

//================================================================
template<class T> class MuCapPACTScan {
public:

  // Does Clark allow book histograms in the constructor?  Perhaps not...
  void init(HistogramFactory &hf, const std::string& hdir, const DetectorGeo& geom, const ConfigFile &conf);

  void fill(const EventClass& evt, const WireCluster& pc5cluster, const WireCluster& pc6cluster);

  MuCapPACTScan() : hClusterSize_() {}

private :
  TH2 *hClusterSize_;
  typedef std::vector<std::vector<T> > Quadrants;
  Quadrants qq_;
};

#endif/*MuCapPACTScan_h*/
