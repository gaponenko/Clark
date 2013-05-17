// Andrei Gaponenko, 2013

#ifndef HistTDCParticleClassifier_h
#define HistTDCParticleClassifier_h

#include <vector>

#include "WireCluster.h"
#include "TimeWindow.h"

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
            const ConfigFile &conf,
            TimeWindow::StreamType stream);

  void fill(const ClustersByPlane& globalClusters);

  HistTDCParticleClassifier()
    : geom_()
    , stream_(TimeWindow::MIXED)
    , hpc8vs7maxWidth_()
    , hpc8vs7meanWidth_()
    , hpc8vs7medianWidth_()
    , hpc8vs7maxHits_()
    , htgtpcmaxWidth_()
    , htgtpcmeanWidth_()
    , htgtpcmedianWidth_()
    , hMaxWidthDC_()
    , hMeanWidthDC_()
    , hMedianWidthDC_()
    , hMaxWidthPC_()
    , hMeanWidthPC_()
    , hMedianWidthPC_()
  {}

private :
  const DetectorGeo *geom_;
  TimeWindow::StreamType stream_;
  TH2 *hpc8vs7maxWidth_;
  TH2 *hpc8vs7meanWidth_;
  TH2 *hpc8vs7medianWidth_;
  TH2 *hpc8vs7maxHits_;
  TH1 *htgtpcmaxWidth_;
  TH1 *htgtpcmeanWidth_;
  TH1 *htgtpcmedianWidth_;

  TH1 *hMaxWidthDC_;
  TH1 *hMeanWidthDC_;
  TH1 *hMedianWidthDC_;

  TH1 *hMaxWidthPC_;
  TH1 *hMeanWidthPC_;
  TH1 *hMedianWidthPC_;
};

#endif/*HistTDCParticleClassifier_h*/
