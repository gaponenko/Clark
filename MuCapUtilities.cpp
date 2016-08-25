#include "MuCapUtilities.h"
#include "EventClass.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <set>

namespace MuCapUtilities {
  //================================================================
  double mass(int pdgId, const EventClass& evt) {
    static const double protonMass =    938.272 /* MeV/c^2 */;
    static const double neutronMass =   939.565 /* MeV/c^2 */;
    static const double deuteronMass = 1875.613 /* MeV/c^2, from G3 */;
    static const double tritonMass =   2809.25  /* MeV/c^2, from G3 */;
    static const double alphaMass =    3727.417 /* MeV/c^2, from G3 */;

    static const double electronMass = 0.510999 /* MeV/c^2 */;
    static const double muonMass = 105.6584 /* MeV/c^2 */;

    switch(evt.mctype) {
    case EventClass::G4: //----------------------------------------------------------------
      switch(std::abs(pdgId)) {

      case PID_G4_PROTON: return protonMass;
      case PID_G4_NEUTRON: return neutronMass;
      case PID_G4_EMINUS: return electronMass;
      case PID_G4_MUMINUS: return muonMass;

      default:
        std::ostringstream os;
        os<<"MuCapUtilities::mass(): unknown G4 pdgId = "<<pdgId;
        throw std::runtime_error(os.str());
      }
      break;

    case EventClass::G3: //----------------------------------------------------------------
      switch(pdgId) {
      case PID_G3_PROTON: return protonMass;
      case PID_G3_DEUTERON: return deuteronMass;
      case PID_G3_TRITON: return tritonMass;
      case PID_G3_ALPHA: return alphaMass;
      case PID_G3_MUMINUS: return muonMass;
      case PID_G3_EMINUS: return electronMass;
      case PID_G3_PHOTON: return 0.;
      case PID_G3_NEUTRON: return neutronMass;

      default:
        std::ostringstream os;
        os<<"MuCapUtilities::mass(): unknown G3 pdgId = "<<pdgId;
        throw std::runtime_error(os.str());
      }
      break;

    default:
      throw std::runtime_error("MuCapUtilities::mass(): unknown evt.mctype");
    }
  }

  //================================================================
  double kineticEnergy(int pdgId, double ptot, const EventClass& evt) {
    const double m = mass(pdgId, evt);
    const double etot = std::sqrt(std::pow(ptot,2) + std::pow(m, 2));
    return etot - m;
  }

  //================================================================
  unsigned numPlanes(const TDCHitWPPtrCollection& hits) {
    std::set<int> planes;
    for(TDCHitWPPtrCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
      planes.insert((**i).plane());
    }
    return planes.size();
  }

  //================================================================
  unsigned numWires(const WireClusterCollection& clusters) {
    unsigned res(0);
    for(WireClusterCollection::const_iterator i = clusters.begin(); i!=clusters.end(); ++i) {
      res  += i->numCells();
    }
    return res;
  }

  //================================================================
  int findExtendedLastPlane(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
    // Find the last plane contiguous with the track
    int lastPlane = evt.hefit_pstop[itrack];
    while(++lastPlane < protonGlobalClusters.size() && !protonGlobalClusters[lastPlane].empty())
      {}

    --lastPlane;
    return lastPlane;
  }

  //================================================================
  int findTrackLastPlane(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
    return evt.hefit_pstop[itrack];
  }

  //================================================================
  Stats::Stats()
    : sum_(0.)
    , min_(std::numeric_limits<double>::max())
    , max_(std::numeric_limits<double>::min())
    , numEntries_(0)
 {}

  //================================================================
  double Stats::mean() const {
    if(!numEntries_) {
      throw std::runtime_error("Stats::mean() called with numEntries==0");
    }
    return sum_/numEntries_;
  }

  //================================================================
  double Stats::median() const {
    if(!values_.size()) {
      throw std::runtime_error("Stats::median() called with numEntries==0");
    }

    // Want to keep median() const, not impose overhead on filling,
    // and do not want to implement caching.
    std::vector<double> tmp(values_);
    std::sort(tmp.begin(), tmp.end());
    unsigned im1 = (tmp.size()-1)/2;
    double res = tmp[im1];
    if(!(tmp.size()%2)) {
      res += tmp[im1+1];
      res /= 2.;
    }
    return res;
  }

  //================================================================
  double Stats::min() const {
    if(!numEntries_) {
      throw std::runtime_error("Stats::min() called with numEntries==0");
    }
    return min_;
  }

  //================================================================
  double Stats::max() const {
    if(!numEntries_) {
      throw std::runtime_error("Stats::max() called with numEntries==0");
    }
    return max_;
  }

  //================================================================
  void Stats::fill(double x) {
    ++numEntries_;
    sum_ += x;
    if(x < min_) min_ = x;
    if(x > max_) max_ = x;
    values_.push_back(x);
  }

  //================================================================
  double computeHitRMax(const TDCHitWPPtrCollection& hits, WirePlane::DetType pt, const DetectorGeo& geo) {
    double rmax = 0.;
    for(unsigned ihit = 0; ihit < hits.size(); ++ihit) {
      const TDCHitWPPtr& hit = hits[ihit];
      const double rhit = std::abs(geo.planes(pt)[hit->plane()].measurement(hit->cell()).coordinate);
      if(rmax < rhit) {
        rmax = rhit;
      }
    }
    return rmax;
  }

  //================================================================
}
