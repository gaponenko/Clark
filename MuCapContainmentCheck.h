// Andrei Gaponenko, 2013

#ifndef MuCapContainmentCheck_h
#define MuCapContainmentCheck_h

#include <string>
#include <limits>

#include "MuCapContainment1D.h"
#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class MuCapContainmentCheck {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  // Computes a bound on R at the next plane
  // throws an exception if the clusters do not allow to extrapolate
  double rmax(int lastPlane,
              const ClustersByPlane& globalPlaneClusters);

  MuCapContainmentCheck()
    : geom_()
    , huv_()
    , hr_()
  {}

private :
  const DetectorGeo *geom_;
  MuCapContainment1D cu_;
  MuCapContainment1D cv_;
  TH2 *huv_;
  TH1 *hr_;
};

#endif/*MuCapContainmentCheck_h*/
