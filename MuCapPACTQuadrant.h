// Andrei Gaponenko, 2013

#ifndef MuCapPACTQuadrant_h
#define MuCapPACTQuadrant_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class WireCluster;

//================================================================
class MuCapPACTQuadrant {
public:

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

  MuCapPACTQuadrant(HistogramFactory &hf, const ConfigFile& conf,
                    const std::string& hdir, const std::string& suffix);

private :
  double slopea_;
  double intercepta_;
  double slopeb_;
  double interceptb_;

  TH2 *hpc6vs5widthAll_;
  TH2 *hpc6vs5widthQ1_;
};

#endif/*MuCapPACTQuadrant_h*/
