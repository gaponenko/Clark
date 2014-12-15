// Andrei Gaponenko, 2014

#ifndef HistHitBasedAnalysis_h
#define HistHitBasedAnalysis_h

#include <string>

#include "HistMuCapTruth.h"
#include "HistMuCapFinal.h"
#include "HistHotSpot.h"
#include "TimeWindow.h"
#include "HistTDCBCSWidth.h"

class TH1;
class TH2;
class TH3;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistHitBasedAnalysis {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  bool accepted(const EventClass& evt, const ClustersByPlane& globalPlaneClusters);

  HistHitBasedAnalysis() : geom_(0), doMCTruth_(false) {}

private :
  const DetectorGeo *geom_;
  bool doMCTruth_;

  TH2* lastconPlaneVsCWires_; // last plane in a contiguous range vs sum(largest cluster size)

  TH3* migration_;

  //----------------------------------------------------------------
  // Extra histograms to plot efficiencies etc.  Not essential for the
  // unfolding.

  // same content as the reco hists, but which one is filled depends on MC PID
  TH2* lastconPlaneVsCWires_mcproton_;
  TH2* lastconPlaneVsCWires_mcdeuteron_;
  TH2* lastconPlaneVsCWires_mcdio_;

  HistMuCapFinal noncontiguous_;

  HistMuCapTruth hTruth_in_;
  HistMuCapTruth hTruth_accepted_;

  TH1* hOuterVetoNumHitPlanes_;
  TH1* hNumPC7Clusters_;
  TH2* hClusterMultiplicity_;

  //----------------------------------------------------------------
  // Understanding the "hot spot" feature in data

  HistHotSpot hshot_;
  HistHotSpot hscold_;

  HistTDCBCSWidth htdcwidthMaxWires_;
  HistTDCBCSWidth htdcwidthMaxTDCWidth_;
  //----------------------------------------------------------------
};

#endif/*HistHitBasedAnalysis_h*/
