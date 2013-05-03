// Andrei Gaponenko, 2013

#ifndef MuCapUtilities_h
#define MuCapUtilities_h

namespace MuCapUtilities {
  // TWIST units are MeV/c^2
  double mass(int pdgId);

  double kineticEnergy(int pdgId, double ptot);
}

#endif/*MuCapUtilities_h*/
