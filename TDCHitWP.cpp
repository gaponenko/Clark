// Andrei Gaponenko, 2013

#include "TDCHitWP.h"

#include <algorithm>
#include <iterator>

std::ostream& operator<<(std::ostream& os, const TDCHitWP& hit) {
  return os<<"TDCHit(plane="<<hit.plane<<", cell="<<hit.cell<<", time="<<hit.time<<", width="<<hit.width<<" )";
}

std::ostream& operator<<(std::ostream& os, const TDCHitWPCollection& coll) {
  std::copy(coll.begin(), coll.end(), std::ostream_iterator<TDCHitWP>(os, " "));
  return os;
}

std::ostream& operator<<(std::ostream& os, const TDCHitWPPtrCollection& hits) {
  for(TDCHitWPPtrCollection::const_iterator i=hits.begin(); i!=hits.end(); ++i) {
    os<<**i<<" ";
  }
  return os;
}
