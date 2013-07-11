// Andrei Gaponenko, 2013

#ifndef HistPlaneRanges_h
#define HistPlaneRanges_h

#include <string>

#include "PlaneRange.h"

class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;


//================================================================
class HistPlaneRanges {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const PlaneRange& r);

  HistPlaneRanges()
    : geom_()
    , hNumRanges_(), hSingleRange_(), hDoubleRangeGap_()
    , hDoubleRangeMissingPlanes_(), hDoubleRangeSizes_()
  {}

private :
  const DetectorGeo *geom_;

  TH1 *hNumRanges_;
  TH2 *hSingleRange_;
  TH2 *hDoubleRangeGap_;
  TH1 *hDoubleRangeMissingPlanes_;
  TH2 *hDoubleRangeSizes_;
};

#endif/*HistPlaneRanges_h*/
