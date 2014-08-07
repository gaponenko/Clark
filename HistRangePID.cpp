// Andrei Gaponenko, 2014

#include "HistRangePID.h"

#include <cmath>

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
void HistRangePID::init(const std::string& hdir,
                        HistogramFactory &hf,
                        const ConfigFile& conf)
{
  planeRangeVsPz_ = hf.DefineTH2D(hdir, "planeRangeVsPz", "Last plane hit vs pz",150, 0., 300., 30, 0.5, 30.5);
  planeRangeVsPz_->SetOption("colz");
  planeRangeVsPz_->GetXaxis()->SetTitle("pz [MeV/c]");
  planeRangeVsPz_->GetYaxis()->SetTitle("plane-28");

  planeRangecosVsP_ = hf.DefineTH2D(hdir, "planeRangecosVsP", "Last plane hit/|cos(theta)| vs p",150, 0., 300., 50, 0., 50.);
  planeRangecosVsP_->SetOption("colz");
  planeRangecosVsP_->GetXaxis()->SetTitle("p [MeV/c]");
  planeRangecosVsP_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

  trackRangeVsPz_ = hf.DefineTH2D(hdir, "trackRangeVsPz", "Last track hit vs pz",150, 0., 300., 30, 0.5, 30.5);
  trackRangeVsPz_->SetOption("colz");
  trackRangeVsPz_->GetXaxis()->SetTitle("pz [MeV/c]");
  trackRangeVsPz_->GetYaxis()->SetTitle("track end-28");

  trackRangecosVsP_ = hf.DefineTH2D(hdir, "trackRangecosVsP", "Last track hit/|cos(theta)| vs p",150, 0., 300., 50, 0., 50.);
  trackRangecosVsP_->SetOption("colz");
  trackRangecosVsP_->GetXaxis()->SetTitle("p [MeV/c]");
  trackRangecosVsP_->GetYaxis()->SetTitle("(track end-28)/|cos(theta)|");
}

//================================================================
void HistRangePID::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
  int trackEnd = evt.hefit_pstop[itrack];
  const double pz = evt.ptot[itrack]*evt.costh[itrack];

  trackRangeVsPz_->Fill(pz, trackEnd-28);
  trackRangecosVsP_->Fill(evt.ptot[itrack], (trackEnd-28)/std::abs(evt.costh[itrack]));

  // Find the last plane contiguous with the track
  int lastPlane = trackEnd;
  while(++lastPlane < protonGlobalClusters.size() && !protonGlobalClusters[lastPlane].empty()) {}
  --lastPlane;
  //std::cout<<"lastPlane - trackEnd = "<<lastPlane - trackEnd<<std::endl;

  planeRangeVsPz_->Fill(pz, lastPlane-28);
  planeRangecosVsP_->Fill(evt.ptot[itrack], (lastPlane-28)/std::abs(evt.costh[itrack]));
}

//================================================================
