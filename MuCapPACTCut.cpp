#include "MuCapPACTCut.h"

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "WireCluster.h"

void MuCapPACTCut::init(HistogramFactory &hf, const ConfigFile& conf) {
  hClusterSize_ = hf.DefineTH2D("MuCapture/PACT", "MuStopClusterSize", "Muon stop cluster V vs U width (cell units)", 6, 0.5, 6.5,  6, 0.5, 6.5);
  hClusterSize_->SetOption("colz");

  qq_.resize(2);

  qq_[0].push_back(MuCapPACTQuadrant(hf, conf, "MuCapture/PACT/q11", "_11"));
  qq_[0].push_back(MuCapPACTQuadrant(hf, conf, "MuCapture/PACT/q12", "_12"));
  qq_[1].push_back(MuCapPACTQuadrant(hf, conf, "MuCapture/PACT/q21", "_21"));
  qq_[1].push_back(MuCapPACTQuadrant(hf, conf, "MuCapture/PACT/q22", "_22"));
}

int MuCapPACTCut::quadrant(const WireCluster& pc5cluster, const WireCluster& pc6cluster) {
  int res = 0;

  hClusterSize_->Fill(pc5cluster.numCells(), pc6cluster.numCells());

  if( (pc5cluster.numCells() <= 2) && (pc6cluster.numCells() <= 2) ) {
    res = qq_.at(pc5cluster.numCells()-1).at(pc6cluster.numCells()-1).quadrant(pc5cluster, pc6cluster);
  }

  return res;
}
