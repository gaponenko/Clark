// Andrei Gaponenko, 2013

#include "WireCluster.h"

#include <algorithm>
#include <iterator>
#include <cassert>

WireCluster::WireCluster(const TDCHitWPPtrCollection& in)
  : hits_(in)
  , plane_()
  , centralCell_()
  , numCells_()
{
  assert(!hits_.empty());
  plane_ = hits_.front()->plane;

  int cmin = hits_.front()->cell;
  int cmax = hits_.front()->cell;
  for(TDCHitWPPtrCollection::const_iterator i= ++hits_.begin(); i!=hits_.end(); ++i) {
    if( (*i)->cell < cmin) {
      cmin = (*i)->cell;
    }
    if(cmax < (*i)->cell) {
      cmax = (*i)->cell;
    }
  }

  centralCell_ = (cmin + cmax)/2.;
  numCells_ = 1 + cmax - cmin;
}

//================================================================
std::ostream& operator<<(std::ostream& os, const WireCluster& cl) {
  return os<<"Cluster(plane="<<cl.plane()<<", hits = "<<cl.hits()<<" )";
}

std::ostream& operator<<(std::ostream& os, const WireClusterCollection& coll) {
  std::copy(coll.begin(), coll.end(), std::ostream_iterator<WireCluster>(os, "\n"));
  return os;
}

std::ostream& operator<<(std::ostream& os, const ClustersByPlane& cbp) {
  for(ClustersByPlane::const_iterator i=cbp.begin(); i!=cbp.end(); ++i) {
    os<<i->second;
  }
  return os;
}
