// Andrei Gaponenko, 2013

#ifndef HistTDCParticleClassifier_h
#define HistTDCParticleClassifier_h

#include <vector>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;

//================================================================
class HistTDCParticleClassifier {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const ClustersByPlane& globalClusters);

  HistTDCParticleClassifier()
    : hpc8vs7maxWidth_()
    , hpc8vs7meanWidth_()
    , hpc8vs7medianWidth_()
    , hpc8vs7maxHits_()
    , hLastVsRestMaxWidth_()
    , hLastVsRestMeanWidth_()
    , hLastVsRestMedianWidth_()
  {}

private :
  TH2 *hpc8vs7maxWidth_;
  TH2 *hpc8vs7meanWidth_;
  TH2 *hpc8vs7medianWidth_;
  TH2 *hpc8vs7maxHits_;

  TH2 *hLastVsRestMaxWidth_;
  TH2 *hLastVsRestMeanWidth_;
  TH2 *hLastVsRestMedianWidth_;
};

#endif/*HistTDCParticleClassifier_h*/
