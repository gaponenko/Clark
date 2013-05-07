#include "TDCHitStats.h"
#include "TDCHitWP.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <limits>

//================================================================
TDCHitStats::TDCHitStats() {}

//================================================================
void TDCHitStats::fill(const TDCHitWPPtr& hit) {
  sw_.fill(hit->width());
  wm_[hit->cid()].push_back(hit);
}

//================================================================
void TDCHitStats::fill(const WireClusterCollection& clusters) {
  for(WireClusterCollection::const_iterator ic = clusters.begin(); ic != clusters.end(); ++ic) {
    for(TDCHitWPPtrCollection::const_iterator ih = ic->hits().begin(); ih != ic->hits().end(); ++ih) {
      fill(*ih);
    }
  }
}

//================================================================
int TDCHitStats::maxHitsPerWire() const {
  int res=0;
  for(WireMap::const_iterator i=wm_.begin(); i!=wm_.end(); ++i) {
    if(res < i->second.size()) res = i->second.size();
  }
  return res;
}

//================================================================
