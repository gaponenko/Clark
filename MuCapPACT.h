// Andrei Gaponenko, 2013

#ifndef MuCapPACT_h
#define MuCapPACT_h

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;

//================================================================
class MuCapPACT {
public:

  // Does Clark allow book histograms in the constructor?  Perhaps not...
  void init(HistogramFactory &hf, const ConfigFile &conf);

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

  int quadrant(const WireCluster& pc5cluster, const WireCluster& pc6cluster);

  MuCapPACT()
    : slopea_()
    , intercepta_()
    , slopeb_()
    , interceptb_()
    , hClusterSize_()
    , hpc6vs5widthAll_()
    , hpc6vs5widthQ1_()
  {}

private :
  double slopea_;
  double intercepta_;
  double slopeb_;
  double interceptb_;

  TH2 *hClusterSize_;
  TH2 *hpc6vs5widthAll_;
  TH2 *hpc6vs5widthQ1_;
};

#endif/*MuCapPACT_h*/
