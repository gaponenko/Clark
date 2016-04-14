#include "MuCapPACTScanSlope.h"

#include <algorithm>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "WireCluster.h"

MuCapPACTScanSlope::MuCapPACTScanSlope(HistogramFactory& hf,
                                       const std::string& topdir,
                                       const DetectorGeo& geom,
                                       const ConfigFile& conf,
                                       const std::string& suffix)
  : slopea_(conf.read<double>("MuCapture/PACTss/slopea_"+suffix))
{
  double slopeb_min = conf.read<double>("MuCapture/PACTss/slopeb_min_"+suffix);
  double slopeb_max = conf.read<double>("MuCapture/PACTss/slopeb_max_"+suffix);
  double slopeb_npoints = conf.read<double>("MuCapture/PACTss/slopeb_npoints_"+suffix);

  double ia_min = conf.read<double>("MuCapture/PACTss/ia_min_"+suffix);
  double ia_max = conf.read<double>("MuCapture/PACTss/ia_max_"+suffix);
  double ia_nbins = conf.read<double>("MuCapture/PACTss/ia_nbins_"+suffix);

  double ib_min = conf.read<double>("MuCapture/PACTss/ib_min_"+suffix);
  double ib_max = conf.read<double>("MuCapture/PACTss/ib_max_"+suffix);
  double ib_nbins = conf.read<double>("MuCapture/PACTss/ib_nbins_"+suffix);

  if(slopeb_npoints < 2) {
    throw std::runtime_error("MuCapPACTScanSlope: bad config slopeb_npoints < 2\n");
  }


  std::string hdir = topdir + "/q"+suffix;

  for(int i=0; i<slopeb_npoints; ++i) {
    const double sb = slopeb_min + i*(slopeb_max - slopeb_min)/(slopeb_npoints - 1);
    slopeb_.push_back(sb);

    std::ostringstream osname;
    osname<<"ia_vs_ib_"<<suffix<<"_pt"<<i;

    std::ostringstream ostitle;
    ostitle<<"PACT ia vs ib "<<suffix<<", sa="<<slopea_<<", sb="<<sb;

    hia_vs_ib_scan_.push_back(hf.DefineTH2D(hdir, osname.str(), ostitle.str(),
                                            ib_nbins, ib_min, ib_max,
                                            ia_nbins, ia_min, ia_max
                                            ));
    hia_vs_ib_scan_.back()->SetOption("colz");
  }
}

void MuCapPACTScanSlope::fill(const EventClass& evt, const WireCluster& pc5cluster, const WireCluster& pc6cluster) {
  const double pc5width = pc5cluster.totalTDCWidth();
  const double pc6width = pc6cluster.totalTDCWidth();

  for(unsigned i=0; i<slopeb_.size(); ++i) {
    const double ia = pc6width - slopea_ * pc5width;
    const double ib = pc6width - slopeb_[i] * pc5width;
    hia_vs_ib_scan_[i]->Fill(ib, ia);
  }
}
