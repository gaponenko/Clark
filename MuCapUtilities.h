// Andrei Gaponenko, 2013

#ifndef MuCapUtilities_h
#define MuCapUtilities_h

#include <vector>
#include "TDCHitWP.h"
#include "WireCluster.h"

class EventClass;

namespace MuCapUtilities {

  // TWIST units are MeV/c^2
  double mass(int pdgId, const EventClass& evt); // evt is to decide between G3 and G4

  double kineticEnergy(int pdgId, double ptot, const EventClass& evt); // evt is to decide between G3 and G4

  unsigned numPlanes(const TDCHitWPPtrCollection& hits);

  unsigned numWires(const WireClusterCollection& clusters);

  //================================================================
  class Stats {
  public:
    void fill(double x);

    double mean() const;
    double median() const;
    double min() const;
    double max() const;

    double minOr(double defaultValue) const { return numEntries_ ? min_ : defaultValue; }
    double maxOr(double defaultValue) const { return numEntries_ ? max_ : defaultValue; }

    int numEntries() const { return numEntries_; }

    Stats();

  private:
    std::vector<double> values_;
    double sum_;
    double min_;
    double max_;
    int numEntries_;
  };

  //================================================================

}

#endif/*MuCapUtilities_h*/
