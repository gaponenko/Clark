#include "MuCapUtilities.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <limits>

namespace MuCapUtilities {
  //================================================================
  double mass(int pdgId) {
    switch(std::abs(pdgId)) {

    case 2212: /*proton*/ return 938.272 /* MeV/c^2 */;

    case 11: /*electron*/ return 0.510999 /* MeV/c^2 */;

    default:
      std::ostringstream os;
      os<<"MuCapUtilities::mass(): unknown pdgId = "<<pdgId;
      throw std::runtime_error(os.str());
    }
  }

  //================================================================
  double kineticEnergy(int pdgId, double ptot) {
    const double m = mass(pdgId);
    const double etot = std::sqrt(std::pow(ptot,2) + std::pow(m, 2));
    return etot - m;
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
  void Stats::fill(double x) {
    ++numEntries_;
    sum_ += x;
    if(x < min_) min_ = x;
    if(x > max_) max_ = x;
  }

  //================================================================
}
