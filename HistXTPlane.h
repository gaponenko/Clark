// Look at cross talk between different planes.
// The inputs are filtered hits from a single time window.
//
// Andrei Gaponenko, 2014

#ifndef HistXTPlane_h
#define HistXTPlane_h

#include <string>
#include <vector>

#include "WireCluster.h"
#include "HistOccupancy.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistXTPlane {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf,
            double cutBroadHitWidth);

  void fill(const EventClass& evt, const ClustersByPlane& globalClusters);

  HistXTPlane() {}

private :
  const DetectorGeo *geom_;
  double cutNarrowHitWidth_[2];
  double cutBroadHitWidth_[2];

  // hdt_[planeWithInducedXt][inducingPlane]
  std::vector<std::vector<TH1*> > hdt_;

  void bookHisto(HistogramFactory& hf, const std::string& hdir, int induced, int inducing);
};

#endif/*HistXTPlane_h*/
