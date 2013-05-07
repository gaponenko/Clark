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
  plane_ = hits_.front()->plane();

  int cmin = hits_.front()->cell();
  int cmax = hits_.front()->cell();
  for(TDCHitWPPtrCollection::const_iterator i= ++hits_.begin(); i!=hits_.end(); ++i) {
    if( (*i)->cell() < cmin) {
      cmin = (*i)->cell();
    }
    if(cmax < (*i)->cell()) {
      cmax = (*i)->cell();
    }
  }

  centralCell_ = (cmin + cmax)/2.;
  numCells_ = 1 + cmax - cmin;
}

//================================================================
double WireCluster::totalTDCWidth() const {
  double res(0);
  for(TDCHitWPPtrCollection::const_iterator i= hits_.begin(); i!=hits_.end(); ++i) {
    res += (*i)->width;
  }
  return res;
}

//================================================================
ClustersByPlane constructPlaneClusters(int maxPlaneNumber, const TDCHitWPPtrCollection& hits) {
  ClustersByPlane res(1+maxPlaneNumber);

  for(unsigned i=0; i<hits.size(); ) {

    TDCHitWPPtrCollection tmphits;
    tmphits.push_back(hits[i]);
    const int plane = hits[i]->plane();

    for(++i; i<hits.size();++i) {
      if(plane != hits[i]->plane()) break;
      if(tmphits.back()->cell() + 1 <  hits[i]->cell()) break;
      tmphits.push_back(hits[i]);
    }

    res[plane].push_back(WireCluster(tmphits));
  }

  return res;
}

//================================================================
ClustersByPlane globalPlaneClusters(const ClustersByPlane& pcClusters, const ClustersByPlane& dcClusters) {
  assert(pcClusters.size() == 13);
  assert(dcClusters.size() == 45);
  ClustersByPlane res(57);

  for(int iGlobal=1; iGlobal <=4; ++iGlobal) {
    res[iGlobal] = pcClusters[iGlobal];
  }

  for(int iGlobal=5; iGlobal <=26; ++iGlobal) {
    res[iGlobal] = dcClusters[iGlobal - 4];
  }

  for(int iGlobal=27; iGlobal <=30; ++iGlobal) {
    res[iGlobal] = pcClusters[iGlobal - 22];
  }

  for(int iGlobal=31; iGlobal <=52; ++iGlobal) {
    res[iGlobal] = dcClusters[iGlobal - 8];
  }

  for(int iGlobal=53; iGlobal <=56; ++iGlobal) {
    res[iGlobal] = pcClusters[iGlobal - 44];
  }

  return res;
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
  os<<"ClustersByPlane: num planes = "<<cbp.size()<<",\n";
  for(unsigned i=0; i < cbp.size(); ++i) {
    os<<"plane "<<i<<" clusters:\n";
    os<<cbp[i];
  }
  return os;
}
