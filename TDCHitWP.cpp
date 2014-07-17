// Andrei Gaponenko, 2013

#include "TDCHitWP.h"

#include <algorithm>
#include <iterator>

std::ostream& operator<<(std::ostream& os, const WireCellId& cid) {
  return os<<"cid("<<cid.plane<<", "<<cid.cell<<")";
}

std::ostream& operator<<(std::ostream& os, const TDCHitWP& hit) {
  return os<<"TDCHit(cid="<<hit.cid()
           <<", time="<<hit.time()
           <<", width="<<hit.width()
           <<", xt="<<hit.xtalk()
           <<" )";
}

std::ostream& operator<<(std::ostream& os, const TDCHitWPCollection& coll) {
  TDCHitWPPtrCollection hits;
  hits.reserve(coll.size());
  for(unsigned i=0; i<coll.size(); ++i) {
    hits.push_back(TDCHitWPPtr(coll, i));
  }

  std::sort(hits.begin(), hits.end(), TDCHitWPCmpTime());
  std::stable_sort(hits.begin(), hits.end(), TDCHitWPCmpGeom());
  // here we have hits sorted by cell then time

  for(TDCHitWPPtrCollection::const_iterator i=hits.begin(); i!=hits.end(); ++i) {
    os<<**i<<" ";
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const TDCHitWPPtrCollection& origHits) {
  TDCHitWPPtrCollection hits(origHits);
  std::sort(hits.begin(), hits.end(), TDCHitWPCmpTime());
  std::stable_sort(hits.begin(), hits.end(), TDCHitWPCmpGeom());
  // here we have hits sorted by cell then time

  for(TDCHitWPPtrCollection::const_iterator i=hits.begin(); i!=hits.end(); ++i) {
    os<<**i<<" ";
  }
  return os;
}
