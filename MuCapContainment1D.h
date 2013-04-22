// Andrei Gaponenko, 2013

#ifndef MuCapContainment1D_h
#define MuCapContainment1D_h

#include <string>

#include "WireCluster.h"
#include "DetectorGeo.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;

//================================================================
class MuCapContainment1D {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  // Returns max(abs(extrapolated coordinate))
  // throws an exception if the clusters do not allow to extrapolate
  WirePlane::Measurement limit(int lastPlane,
                               double extrapolateToZ,
                               const ClustersByPlane& globalClusters);

  MuCapContainment1D()
    : geom_()
    , minPlanedz_()
    , hNumClusters_()
    , hSlopeAll_()
    , hSlopeMulticluster_()
  {}

private :
  const DetectorGeo *geom_;
  double minPlanedz_;
  TH2 *hNumClusters_;
  TH1 *hSlopeAll_;
  TH1 *hSlopeMulticluster_;
};

#endif/*MuCapContainment1D_h*/
