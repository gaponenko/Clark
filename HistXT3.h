// Look at events with an accepted downstream track.
// They guaranteed to have a real track hit in DC23.
//
// Andrei Gaponenko, 2017

#ifndef HistXT3_h
#define HistXT3_h

#include <string>
#include <vector>

#include "TDCHitWP.h"
#include "WireCluster.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class ConfigFile;
class EventClass;


//================================================================
// a helper class
class HistXT3TwoCluster {
public:
  void init(const std::string& hdir, const std::string& suffix, HistogramFactory &hf, const ConfigFile &conf);
  void fill(const WireClusterCollection& planeClusters);
private :
  TH1 *twocluster_dt_;
  TH2 *twocluster_ww_;
  TH2 *twocluster_wt_;
  TH2 *twocluster_ss_;

  // loop over hits in the second cluster
  TH1 *hitcluster_dt_all_;
  TH2 *hitcluster_ww_all_;
  TH2 *hitcluster_wt_all_;

  // loop over hits in the second cluster, fill for XT hits
  TH1 *hitcluster_dt_xt_;
  TH2 *hitcluster_ww_xt_;
  TH2 *hitcluster_wt_xt_;

  // loop over hits in the second cluster, fill for non-XT hits
  TH1 *hitcluster_dt_nxt_;
  TH2 *hitcluster_ww_nxt_;
  TH2 *hitcluster_wt_nxt_;
};

//================================================================
class HistXT3 {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const EventClass& evt, int iTrack, const ClustersByPlane& clarkClusters);

  HistXT3() {}

private :
  //----------------------------------------------------------------
  // preliminaries
  TH1 *trackCosth_;
  TH1 *trackAngle_;
  TH1 *trackMomentum_;

  TH1 *inputHitMultiplicity23_;
  TH1 *inputHitMultiplicity23xt_;
  TH1 *inputHitMultiplicity23nxt_;

  TH1 *hitWidthAll_;
  TH1 *hitWidthXT_;
  TH1 *hitWidthNXT_;

  TH1 *hitTrackTimeAll_;

  TH1 *inputWireMultiplicity23_;

  TH1 *inputClusterMultiplicity23_;

  TH1 *clusterMultiplicity23All_;
  TH1 *clusterMultiplicity23NXT_;
  TH1 *clusterMultiplicity23Clark_;

  //----------------------------------------------------------------
  // Real hit params
  TH1 *singleHit23Dt_; // t(hit)-t(track) for unique DC23 hits
  TH1 *singleHit23Width_; // width(hit)

  // likely real hit params
  TH1 *hits23MaxWidth_; // max width(hits in plane)
  TH1 *hits23FirstWidth_; // width(earliest hit in plane)

  //----------------------------------------------------------------
  // multicluster preliminaries

  TH2 *doubleClusterSizeAll_;
  TH2 *doubleClusterSizeNXT_;
  TH2 *doubleClusterSizeClark_;

  //----------------
  // 2-cluster distributions

  HistXT3TwoCluster twoClusterAll_;
  HistXT3TwoCluster twoClusterNXT_;

  //----------------------------------------------------------------

};

#endif/*HistXT3_h*/
