// Andrei Gaponenko, 2014

#include "HistTDCSinglePlanePID.h"

#include <limits>
#include <sstream>
#include <cmath>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
TDCPlanePIDResult::TDCPlanePIDResult()
  : analyzed(false)
  , nwires(0)
  , raw(std::numeric_limits<double>::quiet_NaN())
  , calibrated(std::numeric_limits<double>::quiet_NaN())
 {}

//================================================================
void HistTDCSinglePlanePID::init(const std::string& hdir,
                                 int globalPlaneNumber,
                                 HistogramFactory &hf,
                                 const ConfigFile& conf)
{
  globalPlaneNumber_ = globalPlaneNumber;

  std::ostringstream osplane;
  osplane<<globalPlaneNumber;

  hNumClusters_ = hf.DefineTH1D(hdir, "numClusters", "num clusters in plane "+osplane.str(), 9, -0.5, 8.5);
  hNumClusters_->SetOption("colz");

  hSingleClusterSize_ = hf.DefineTH1D(hdir, "singleClusterSize", "cluster size in plane "+osplane.str()+" for single cluster events", 9, -0.5, 8.5);
  hSingleClusterSize_->SetOption("colz");

  hsumwcos_vs_p_1_ = hf.DefineTH2D(hdir, "sumwcos1", "Plane "+osplane.str()+" TDC mean width x cos(theta) vs p, 1", 150, 0., 300., 160, 0, 800.);
  hsumwcos_vs_p_1_->SetOption("colz");
  hsumwcos_vs_p_1_->GetXaxis()->SetTitle("momentum [MeV/c]");
  hsumwcos_vs_p_1_->GetYaxis()->SetTitle("cos(theta) * mean TDC width");

  hsumwcos_vs_p_2_ = hf.DefineTH2D(hdir, "sumwcos2", "Plane "+osplane.str()+" TDC mean width s cos(theta) vs p, 2", 150, 0., 300., 160, 0, 800.);
  hsumwcos_vs_p_2_->SetOption("colz");
  hsumwcos_vs_p_2_->GetXaxis()->SetTitle("momentum [MeV/c]");
  hsumwcos_vs_p_2_->GetYaxis()->SetTitle("cos(theta) * mean TDC width");

  hsumwcos_vs_p_3_ = hf.DefineTH2D(hdir, "sumwcos3", "Plane "+osplane.str()+" TDC mean width s cos(theta) vs p, 3", 150, 0., 300., 160, 0, 800.);
  hsumwcos_vs_p_3_->SetOption("colz");
  hsumwcos_vs_p_3_->GetXaxis()->SetTitle("momentum [MeV/c]");
  hsumwcos_vs_p_3_->GetYaxis()->SetTitle("cos(theta) * mean TDC width");
}

//================================================================
TDCPlanePIDResult HistTDCSinglePlanePID::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
  TDCPlanePIDResult res;

  hNumClusters_->Fill(protonGlobalClusters[globalPlaneNumber_].size());
  const WireClusterCollection& clusters = protonGlobalClusters[globalPlaneNumber_];

  if(clusters.size() == 1) {
    hSingleClusterSize_->Fill(clusters[0].numCells());

    // Use mean value per cell, not per hit
    double meanWidth = clusters[0].totalTDCWidth()/clusters[0].numCells();
    const double rawPIDVar = meanWidth * evt.costh[itrack];

    res.analyzed = true;
    res.nwires = clusters[0].numCells();
    res.raw = rawPIDVar;

    if(clusters[0].numCells() == 1) {
      hsumwcos_vs_p_1_->Fill(evt.ptot[itrack], rawPIDVar);
    }
    else if(clusters[0].numCells() == 2) {
      hsumwcos_vs_p_2_->Fill(evt.ptot[itrack], rawPIDVar);
    }
    else if(clusters[0].numCells() == 3) {
      hsumwcos_vs_p_3_->Fill(evt.ptot[itrack], rawPIDVar);
    }
  }

  return res;
}

//================================================================
