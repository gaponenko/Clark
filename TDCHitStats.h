// Andrei Gaponenko, 2013

#ifndef TDCHitStats_h
#define TDCHitStats_h

#include <vector>
#include <map>

#include "MuCapUtilities.h"
#include "TDCHitWP.h"
#include "WireCluster.h"

class TDCHitStats {
public:
  typedef std::map<WireCellId,TDCHitWPPtrCollection> WireMap;

  void fill(const TDCHitWPPtr& hit);

  const MuCapUtilities::Stats& widthStats() const { return sw_; }

  const WireMap& hitsByWire() const { return wm_; }

  int maxHitsPerWire() const;

private:
  MuCapUtilities::Stats sw_;
  WireMap wm_;
};

//================================================================
class WireClusterStats {
public:
  void fill(const WireClusterCollection& clusters);
  const TDCHitStats& hitStats() const { return hitStats_; }
  const MuCapUtilities::Stats& clusterSizeStats() const { return clusterSizeStats_; }

private:
  TDCHitStats hitStats_;
  MuCapUtilities::Stats clusterSizeStats_;
};

#endif/*TDCHitStats_h*/
