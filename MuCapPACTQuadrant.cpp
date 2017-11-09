#include "MuCapPACTQuadrant.h"

#include <algorithm>

#include "TH2.h"

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"
#include "WireCluster.h"

MuCapPACTQuadrant::MuCapPACTQuadrant(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile& conf,
                                     const std::string& hdir, const std::string& suffix)
  : slopea_(conf.read<double>("MuCapture/PACT/slopea"+suffix))
  , intercepta_(conf.read<double>("MuCapture/PACT/intercepta"+suffix))
  , slopeb_(conf.read<double>("MuCapture/PACT/slopeb"+suffix))
  , interceptb_(conf.read<double>("MuCapture/PACT/interceptb"+suffix))
  , hpc6vs5widthAll_()
  , hpc6vs5widthQ1_()
{
  hpc6vs5widthAll_ = hf.DefineTH2D(hdir, "pc6vs5widthAll"+suffix, "PC6 vs 5 TDC width, all"+suffix, 500, 0., 500.,  500, 0., 500.);
  hpc6vs5widthAll_->SetOption("colz");

  hpc6vs5widthQ1_ = hf.DefineTH2D(hdir, "pc6vs5widthQ1"+suffix, "PC6 vs 5 TDC width, quadrant 1"+suffix, 500, 0., 500.,  500, 0., 500.);
  hpc6vs5widthQ1_->SetOption("colz");
}

int MuCapPACTQuadrant::quadrant(const WireCluster& pc5cluster, const WireCluster& pc6cluster) {

  const double pc5width = pc5cluster.totalTDCWidth();
  const double pc6width = pc6cluster.totalTDCWidth();
  hpc6vs5widthAll_->Fill(pc5width, pc6width);

  const double linea = intercepta_ + slopea_ * pc5width;
  const double lineb = interceptb_ + slopeb_ * pc5width;

  int res = (pc6width < linea) ?
    ((pc6width < lineb)? 4 : 3):
    ((pc6width < lineb)? 1 : 2);

  if(res == 1) {
    hpc6vs5widthQ1_->Fill(pc5width, pc6width);
  }

  return res;
}
