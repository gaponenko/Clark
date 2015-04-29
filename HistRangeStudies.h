// Various plots to understand the range vs p_rec discrepancy in the
// contained channel.
//
// Andrei Gaponenko, 2015

#ifndef HistRangeStudies_h
#define HistRangeStudies_h

#include <string>

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
class HistRangeStudies {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const EventClass& evt, int itrack, const ClustersByPlane& globalPlaneClusters);

private :
  TH2* extendedVsTrackRange_;
  TH2* rangeDiffVsPrec_;
  TProfile2D *clusterMultiplicity_;

  TH1* tdcWidthTrack_;
  TH1* tdcWidthExtended_;
  HistOccupancy hocc_;
};

#endif/*HistRangeStudies_h*/
