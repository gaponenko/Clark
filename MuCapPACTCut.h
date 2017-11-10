// Andrei Gaponenko, 2013

#ifndef MuCapPACTCut_h
#define MuCapPACTCut_h

#include <vector>
#include "MuCapPACTQuadrant.h"

class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class WireCluster;
class EventClass;

//================================================================
class MuCapPACTCut {
public:

  // Does Clark allow book histograms in the constructor?  Perhaps not...
  void init(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile &conf);

  //----------------------------------------------------------------
  // Useful plot for understanding the PACT quadrants
  //
  // y (pc6)  b       a
  // ^         \  2  /
  // |          \   /
  // |           \ /
  // |        1   X  3
  // |           / \                 .
  // |          /   \                .
  // |         /  4  \               .
  // |
  //  ------------------------------->x (pc5)
  //
  // This method returns a quadrant 1-4 as defined above, or 0 if
  // can't determine (e.g. too many hit wires per plane).

  int quadrant(const WireCluster& pc5cluster, const WireCluster& pc6cluster, const EventClass& evt);

  MuCapPACTCut() : hClusterSize_() {}

private :
  TH2 *hClusterSize_;
  typedef std::vector<std::vector<MuCapPACTQuadrant> > Quadrants;
  Quadrants qq_;
};

#endif/*MuCapPACTCut_h*/
