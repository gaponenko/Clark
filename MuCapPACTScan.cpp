#include "MuCapPACTScan.h"

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "WireCluster.h"

#include "MuCapPACTScanQuadrant.h"
#include "MuCapPACTScanSlope.h"

template<class T>
void MuCapPACTScan<T>::init(HistogramFactory &hf, const std::string& hdir, const DetectorGeo& geom, const ConfigFile& conf) {
  hClusterSize_ = hf.DefineTH2D(hdir, "MuStopClusterSize", "Muon stop cluster V vs U width (cell units)", 6, 0.5, 6.5,  6, 0.5, 6.5);
  hClusterSize_->SetOption("colz");

  qq_.resize(2);

  qq_[0].push_back(T(hf, hdir, geom, conf, "11"));
  qq_[0].push_back(T(hf, hdir, geom, conf, "12"));
  qq_[1].push_back(T(hf, hdir, geom, conf, "21"));
  qq_[1].push_back(T(hf, hdir, geom, conf, "22"));
}

template<class T>
void MuCapPACTScan<T>::fill(const EventClass& evt, const WireCluster& pc5cluster, const WireCluster& pc6cluster) {
  hClusterSize_->Fill(pc5cluster.numCells(), pc6cluster.numCells());

  if( (pc5cluster.numCells() <= 2) && (pc6cluster.numCells() <= 2) ) {
    qq_.at(pc5cluster.numCells()-1).at(pc6cluster.numCells()-1).fill(evt, pc5cluster, pc6cluster);
  }
}

template class MuCapPACTScan<MuCapPACTScanQuadrant>;
template class MuCapPACTScan<MuCapPACTScanSlope>;
