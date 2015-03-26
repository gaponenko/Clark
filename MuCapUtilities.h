// Andrei Gaponenko, 2013

#ifndef MuCapUtilities_h
#define MuCapUtilities_h

#include <vector>
#include "TDCHitWP.h"
#include "WireCluster.h"
#include "DetectorGeo.h"

class EventClass;

namespace MuCapUtilities {

  static int const PID_G3_PROTON = 14;
  static int const PID_G3_DEUTERON = 45;
  static int const PID_G3_MUMINUS = 66;

  static int const PID_G4_MUMINUS = +13;
  static int const PID_G4_EMINUS = +11;
  static int const PID_G4_PROTON = 2212;
  static int const PID_G4_NEUTRON = 2112;

  static int const PROC_G4_PRIMARY = 56; // ProcessCode::mu2ePrimary

  // TWIST units are MeV/c^2
  double mass(int pdgId, const EventClass& evt); // evt is to decide between G3 and G4

  double kineticEnergy(int pdgId, double ptot, const EventClass& evt); // evt is to decide between G3 and G4

  unsigned numPlanes(const TDCHitWPPtrCollection& hits);

  unsigned numWires(const WireClusterCollection& clusters);

  // Find the last plane in a hit plane range contiguous with the track
  int findExtendedLastPlane(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters);

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
  double computeHitRMax(const TDCHitWPPtrCollection& hits, WirePlane::DetType pt, const DetectorGeo& geo);

  //================================================================

}

#endif/*MuCapUtilities_h*/
