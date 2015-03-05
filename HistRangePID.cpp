// Andrei Gaponenko, 2014

#include "HistRangePID.h"
#include "MuCapUtilities.h"

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

  planeRangeVsP_ = hf.DefineTH2D(hdir, "planeRangeVsP", "Last plane hit vs p",150, 0., 300., 30, 0.5, 30.5);
  planeRangeVsP_->SetOption("colz");
  planeRangeVsP_->GetXaxis()->SetTitle("p [MeV/c]");
  planeRangeVsP_->GetYaxis()->SetTitle("plane-28");

  planeRangecosVsP_ = hf.DefineTH2D(hdir, "planeRangecosVsP", "Last plane hit/|cos(theta)| vs p",150, 0., 300., 50, 0., 50.);
  planeRangecosVsP_->SetOption("colz");
  planeRangecosVsP_->GetXaxis()->SetTitle("p [MeV/c]");
  planeRangecosVsP_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

  trackRangeVsPz_ = hf.DefineTH2D(hdir, "trackRangeVsPz", "Last track hit vs pz",150, 0., 300., 30, 0.5, 30.5);
  trackRangeVsPz_->SetOption("colz");
  trackRangeVsPz_->GetXaxis()->SetTitle("pz [MeV/c]");
  trackRangeVsPz_->GetYaxis()->SetTitle("track end-28");

  trackRangeVsP_ = hf.DefineTH2D(hdir, "trackRangeVsP", "Last track hit vs p",150, 0., 300., 30, 0.5, 30.5);
  trackRangeVsP_->SetOption("colz");
  trackRangeVsP_->GetXaxis()->SetTitle("p [MeV/c]");
  trackRangeVsP_->GetYaxis()->SetTitle("track end-28");

  trackRangecosVsP_ = hf.DefineTH2D(hdir, "trackRangecosVsP", "Last track hit/|cos(theta)| vs p",150, 0., 300., 50, 0., 50.);
  trackRangecosVsP_->SetOption("colz");
  trackRangecosVsP_->GetXaxis()->SetTitle("p [MeV/c]");
  trackRangecosVsP_->GetYaxis()->SetTitle("(track end-28)/|cos(theta)|");
}

//================================================================
double HistRangePID::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
  int trackEnd = evt.hefit_pstop[itrack];
  const double ptot = evt.ptot[itrack];
  const double pz = ptot*evt.costh[itrack];

  trackRangeVsPz_->Fill(pz, trackEnd-28);
  trackRangeVsP_->Fill(ptot, trackEnd-28);
  trackRangecosVsP_->Fill(evt.ptot[itrack], (trackEnd-28)/std::abs(evt.costh[itrack]));

  // Find the last plane contiguous with the track
  int lastPlane = MuCapUtilities::findExtendedLastPlane(evt, itrack, protonGlobalClusters);
  planeRangeVsPz_->Fill(pz, lastPlane-28);
  planeRangeVsP_->Fill(ptot, lastPlane-28);
  const double planeRangecosPIDVar = (lastPlane-28)/std::abs(evt.costh[itrack]);
  planeRangecosVsP_->Fill(evt.ptot[itrack], planeRangecosPIDVar);

  return planeRangecosPIDVar;
}

//================================================================
