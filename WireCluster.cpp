// Andrei Gaponenko, 2013

#include "WireCluster.h"

#include <algorithm>
#include <iterator>

std::ostream& operator<<(std::ostream& os, const WireCluster& cl) {
  return os<<"Cluster(plane="<<cl.plane<<", hits = "<<cl.hits<<" )";
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
