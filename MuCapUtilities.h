// Andrei Gaponenko, 2013

#ifndef MuCapUtilities_h
#define MuCapUtilities_h

#include <vector>
#include "TDCHitWP.h"

namespace MuCapUtilities {

  // TWIST units are MeV/c^2
  double mass(int pdgId);

  double kineticEnergy(int pdgId, double ptot);

  unsigned numPlanes(const TDCHitWPPtrCollection& hits);

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
