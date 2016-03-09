#include "MuCapPACTScanQuadrant.h"

#include <algorithm>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "WireCluster.h"



MuCapPACTScanQuadrant::MuCapPACTScanQuadrant(HistogramFactory& hf,
                                             const std::string& topdir,
                                             const DetectorGeo& geom,
                                             const ConfigFile& conf,
                                             const std::string& suffix)
  : slopea_(conf.read<double>("MuCapture/PACT/slopea_"+suffix))
  , intercepta_(conf.read<double>("MuCapture/PACT/intercepta_"+suffix))

  , slopeb_(conf.read<double>("MuCapture/PACT/slopeb_"+suffix))
  , interceptb_min_(conf.read<double>("MuCapture/PACT/interceptb_min_"+suffix))
  , interceptb_max_(conf.read<double>("MuCapture/PACT/interceptb_max_"+suffix))
  , interceptb_npoints_(conf.read<double>("MuCapture/PACT/interceptb_npoints_"+suffix))

  , hpc6vs5widthAll_()
{
  std::string hdir = topdir + "/q"+suffix;

  hpc6vs5widthAll_ = hf.DefineTH2D(hdir, "pc6vs5widthAll"+suffix, "PC6 vs 5 TDC width, all "+suffix, 500, 0., 500.,  500, 0., 500.);
  hpc6vs5widthAll_->SetOption("colz");

  if(interceptb_npoints_ < 2) {
    throw std::runtime_error("MuCapPACTScanQuadrant: bad config interceptb_npoints < 2\n");
  }

  for(int i=0; i<interceptb_npoints_; ++i) {
    const double ib = interceptb_min_ + i*(interceptb_max_ - interceptb_min_)/(interceptb_npoints_ - 1);
    interceptb_.push_back(ib);

    std::ostringstream osname;
    osname<<"pc6vs5widthQ1_"<<suffix<<"_pt"<<i;

    std::ostringstream ostitle;
    ostitle<<"Accepted PC6 vs 5 TDC width "<<suffix<<", ib="<<ib;

    hpc6vs5widthQ1_scanb_.push_back(hf.DefineTH2D(hdir, osname.str(), ostitle.str(), 500, 0., 500.,  500, 0., 500.));
    hpc6vs5widthQ1_scanb_.back()->SetOption("colz");

    std::ostringstream mcdir;
    mcdir << hdir <<"/"<<i;
    hMcMuStops_scanb_.push_back(new HistMuCapMCTgtStops());
    hMcMuStops_scanb_.back()->init(hf, mcdir.str(), geom, conf);
  }
}

void MuCapPACTScanQuadrant::fill(const EventClass& evt, const WireCluster& pc5cluster, const WireCluster& pc6cluster) {

  const double pc5width = pc5cluster.totalTDCWidth();
  const double pc6width = pc6cluster.totalTDCWidth();
  hpc6vs5widthAll_->Fill(pc5width, pc6width);

  for(unsigned i=0; i<interceptb_.size(); ++i) {
    const double linea = intercepta_ + slopea_ * pc5width;
    const double lineb = interceptb_[i] + slopeb_ * pc5width;

    if (pc6width > linea && pc6width < lineb) {
      // accepted - quadrant 1
      hpc6vs5widthQ1_scanb_[i]->Fill(pc5width, pc6width);
      hMcMuStops_scanb_[i]->fill(evt);
    }
  }

}
