// Why does the max TDC width distribution in DC23
// has a peak at ~200ns even if we fill it only for
// presumably good proton (or heavier) tracks?
//
// Andrei Gaponenko, 2017

#ifndef Hist200nsStudy_h
#define Hist200nsStudy_h

#include <string>
#include <vector>

#include "TDCHitWP.h"
#include "WireCluster.h"

#include "MuCapContainedVars.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class EventClass;

//================================================================
class Hist200nsStudy {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const EventClass& evt, int iTrack, const ClustersByPlane& clarkClusters);

  Hist200nsStudy() {}

private :
  MuCapContainedVars::RangeCosVsP cvp_;

  TH1 *hitTrackdtDC_;
  TH1 *hitTrackdtPC_;

  //----------------------------------------------------------------
  // Proton/Deuteron PID from range vs momentum.

  TH2 *rangemom_contained_;
  TH2 *rangemom_proton_; // events that we call protons
  TH2 *rangemom_deuteron_; // events that we call deuterons

  //----------------
  TH1 *hits23MaxWidth_contained_; // max width(hits in plane)
  TH1 *hits23MaxWidth_proton_; // max width(hits in plane)
  TH1 *hits23MaxWidth_deuteron_; // max width(hits in plane)

  //----------------------------------------------------------------
  TH2 *hitsMaxWidth24vs23_all_;
  TH2 *hitsMaxWidth24vs23_contained_;

  TH2 *hitsMaxWidthPC8vsDC23_all_;
  TH2 *hitsMaxWidthPC8vsDC23_contained_;

  //----------------------------------------------------------------
};

#endif/*Hist200nsStudy_h*/
