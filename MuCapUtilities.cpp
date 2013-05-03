#include "MuCapUtilities.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

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
}
