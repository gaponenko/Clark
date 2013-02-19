#include "MuCapPACT.h"

#include <algorithm>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

void MuCapPACT::init(HistogramFactory &hf, const ConfigFile& conf) {
  slopea_ = conf.read<double>("MuCapture/PACT/slopea");
  slopeb_ = conf.read<double>("MuCapture/PACT/slopeb");
  intercepta_ = conf.read<double>("MuCapture/PACT/intercepta");
  interceptb_ = conf.read<double>("MuCapture/PACT/interceptb");

  hClusterSize_ = hf.DefineTH2D("MuCapture/PACT", "MuStopClusterSize", "Muon stop cluster V vs U width (cell units)", 5, 0.5, 5.5,  5, 0.5, 5.5);

  hpc6vs5widthAll_ = hf.DefineTH2D("MuCapture/PACT", "pc6vs5widthAll", "PC6 vs 5 TDC width, all", 500, 0., 500.,  500, 0., 500.);
  hpc6vs5widthQ1_ = hf.DefineTH2D("MuCapture/PACT", "pc6vs5widthQ1", "PC6 vs 5 TDC width, quadrant 1", 500, 0., 500.,  500, 0., 500.);
}

int MuCapPACT::quadrant(const WireCluster& pc5cluster, const WireCluster& pc6cluster) {
  int res = 0;

  hClusterSize_->Fill(pc5cluster.numCells(), pc6cluster.numCells());

  if( (pc5cluster.numCells() <= 2) && (pc6cluster.numCells() <= 2) ) {
    const double pc5width = pc5cluster.totalTDCWidth();
    const double pc6width = pc6cluster.totalTDCWidth();
    hpc6vs5widthAll_->Fill(pc5width, pc6width);

    const double linea = intercepta_ + slopea_ * pc5width;
    const double lineb = interceptb_ + slopeb_ * pc5width;

    if (pc6width > linea && pc6width > lineb) res = 2;
    if (pc6width > linea && pc6width < lineb) res = 1;
    if (pc6width < linea && pc6width > lineb) res = 3;
    if (pc6width < linea && pc6width < lineb) res = 4;

    if(res == 1) {
      hpc6vs5widthQ1_->Fill(pc5width, pc6width);
    }
  }

  return res;
}
