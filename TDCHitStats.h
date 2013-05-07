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
  void fill(const WireClusterCollection& clusters);

  const MuCapUtilities::Stats& widthStats() const { return sw_; }

  const WireMap& hitsByWire() const { return wm_; }

  int maxHitsPerWire() const;

  TDCHitStats();

private:
  MuCapUtilities::Stats sw_;
  WireMap wm_;
};

#endif/*TDCHitStats_h*/
