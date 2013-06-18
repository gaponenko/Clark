// Andrei Gaponenko, 2013

#ifndef HistHitStream_h
#define HistHitStream_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;


//================================================================
class HistHitStream {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const ClustersByPlane& globalClusters);

  HistHitStream()
    : geom_()
    , hNumPlanesUpVsDn_() , hNumPlanesDnMinusUp_()
    , hNumPCsUpVsDn_() , hNumPCsDnMinusUp_()
    , hNumDCsUpVsDn_() , hNumDCsDnMinusUp_()
  {}

private :
  const DetectorGeo *geom_;

  TH2 *hNumPlanesUpVsDn_;
  TH1 *hNumPlanesDnMinusUp_;

  TH2 *hNumPCsUpVsDn_;
  TH1 *hNumPCsDnMinusUp_;

  TH2 *hNumDCsUpVsDn_;
  TH1 *hNumDCsDnMinusUp_;
};

#endif/*HistHitStream_h*/
