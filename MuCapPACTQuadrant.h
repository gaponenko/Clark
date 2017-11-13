// Andrei Gaponenko, 2013

#ifndef MuCapPACTQuadrant_h
#define MuCapPACTQuadrant_h

#include <string>
#include <random>

class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class WireCluster;
class EventClass;

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
  int quadrant(const WireCluster& pc5cluster, const WireCluster& pc6cluster, const EventClass& evt);

  MuCapPACTQuadrant(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile& conf,
                    const std::string& hdir, const std::string& suffix);

private :
  double slopea_;
  double intercepta_;
  double slopeb_;
  double interceptb_;

  double extrasmearing_; // for systematic studies

  TH2 *hpc6vs5widthAll_;
  TH2 *hpc6vs5widthQ1_;

  // transformed coordinates:
  // (x,y) => (dib, dia) == (y - sb*x - ib, y - sa*x - ia)

  TH2 *dia_vs_dib_all_;
  TH2 *dia_vs_dib_q1_;

  //----------------------------------------------------------------
  bool doMCTruth_;
  double targetCenterZ_;
  double targetThickness_;
  double pc6CenterZ_;
  double pc6wireRadius_;

  TH2 *mctruthTargetStops_;
  TH2 *mctruthWireStops_;
  TH2 *mctruthOtherStops_;

  TH2 *mctruthTargetStopsi_;
  TH2 *mctruthWireStopsi_;
  TH2 *mctruthOtherStopsi_;

  enum class MuStopRegion { TARGET, PC6WIRE, OTHER };
  MuStopRegion muStopKind(const EventClass& evt) const;

  //----------------------------------------------------------------
  int    lastSeededRun_; // for reproducibility of random numbers
  std::mt19937 eng_;
  std::normal_distribution<double> gaus_;

  double smearPCWidth(const EventClass& evt, double origWidth);
  //----------------------------------------------------------------
};

#endif/*MuCapPACTQuadrant_h*/
