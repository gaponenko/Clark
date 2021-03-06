// Andrei Gaponenko, 2015

#include "HistRangeStudies.h"
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include "HitBasedObservables.h"
#include "EventClass.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

//================================================================
void HistRangeStudies::init(HistogramFactory& hf,
                            const std::string& hdir,
                            const DetectorGeo& geom,
                            const ConfigFile& conf)
{
  //----------------------------------------------------------------
  extendedVsTrackRange_ = hf.DefineTH2D(hdir, "extendedVsTrackRange", "Extended vs track range", 28, 0.5, 28.5, 28, 0.5, 28.5);
  extendedVsTrackRange_->SetOption("colz");
  extendedVsTrackRange_->GetXaxis()->SetTitle("track range");
  extendedVsTrackRange_->GetYaxis()->SetTitle("extended range");

  //----------------------------------------------------------------
  rangeDiffVsPrec_ = hf.DefineTH2D(hdir, "rangeDiff", "Range diff vs momentum", 88, 30., 250., 28, -0.5, 27.5);
  rangeDiffVsPrec_->SetOption("colz");
  rangeDiffVsPrec_->GetXaxis()->SetTitle("p [MeV/c]");
  rangeDiffVsPrec_->GetYaxis()->SetTitle("range diff");

  //----------------------------------------------------------------
  clusterMultiplicity_ = hf.DefineTProfile2D(hdir, "clusterMultiplicity", "cluster multiplicity - 1  in track range profile", 88, 30., 250., 28, 28.5, 56.5);
  clusterMultiplicity_->SetOption("colz");
  clusterMultiplicity_->GetXaxis()->SetTitle("p [MeV/c]");
  clusterMultiplicity_->GetYaxis()->SetTitle("plane");

  //----------------------------------------------------------------
  tdcWidthTrack_ = hf.DefineTH1D(hdir, "tdcWidthTrack", "DC TDC width for hits in track range", 400, 0., 10000.);
  tdcWidthExtended_ = hf.DefineTH1D(hdir, "tdcWidthExtended", "DC TDC width for hits adjacent to track range", 400, 0., 10000.);

  tdcMaxWidthTrack_ = hf.DefineTH1D(hdir, "tdcMaxWidthTrack", "DC TDC max  width for hits in track range", 400, 0., 10000.);
  tdcMaxWidthExtended_ = hf.DefineTH1D(hdir, "tdcMaxWidthExtended", "DC TDC max width for hits adjacent to track range", 400, 0., 10000.);

  tdcMaxWidthDenseTrack_ = hf.DefineTH1D(hdir, "tdcMaxWidthDenseTrack", "DC TDC max width for hits in track range, DS", 400, 0., 10000.);
  tdcMaxWidthDenseExtended_ = hf.DefineTH1D(hdir, "tdcMaxWidthDenseExtended", "DC TDC max width for hits adjacent to track range, DS", 400, 0., 10000.);

  hocc_.init(hdir, "extendedHitMapDC", 44, 80, hf, conf);

}

//================================================================
void HistRangeStudies::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
  const double ptot = evt.ptot[itrack];

  const int extended = MuCapUtilities::findExtendedLastPlane(evt, itrack, protonGlobalClusters);
  const int trange = MuCapUtilities::findTrackLastPlane(evt, itrack, protonGlobalClusters);

  extendedVsTrackRange_->Fill(trange-28, extended-28);
  rangeDiffVsPrec_->Fill(ptot, extended-trange);

  for(int iplane=29; iplane<=trange; ++iplane) {
    const int numClusters= protonGlobalClusters.at(iplane).size();
    clusterMultiplicity_->Fill(ptot, iplane, numClusters - 1);
  }

  //----------------------------------------------------------------
  // Process hits in range
  for(int iplane = evt.hefit_pstart[itrack]; iplane <= evt.hefit_pstop[itrack]; ++iplane) {
    const WireClusterCollection& planeClusters = protonGlobalClusters[iplane];
    double maxWidthInPlane = 0.;
    for(int icluster=0; icluster<planeClusters.size(); ++icluster) {
      const TDCHitWPPtrCollection& hits = planeClusters[icluster].hits();
      for(int ihit=0; ihit<hits.size(); ++ihit) {
        tdcWidthTrack_->Fill(hits[ihit]->width());
        maxWidthInPlane = std::max<double>(hits[ihit]->width(), maxWidthInPlane);
      }
    }
    tdcMaxWidthTrack_->Fill(maxWidthInPlane);
    // 56 - 4(pc) - 8(dense) = 44
    if((44 <= iplane) && (iplane <= 52)) {
      tdcMaxWidthDenseTrack_->Fill(maxWidthInPlane);
    }
  }

  // Process extended hits
  for(int iplane=1+evt.hefit_pstop[itrack];
      (iplane<protonGlobalClusters.size()) && !protonGlobalClusters[iplane].empty();
      ++iplane) {

    const WireClusterCollection& planeClusters = protonGlobalClusters[iplane];
    double maxWidthInPlane = 0.;
    for(int icluster=0; icluster<planeClusters.size(); ++icluster) {
      const TDCHitWPPtrCollection& hits = planeClusters[icluster].hits();
      for(int ihit=0; ihit<hits.size(); ++ihit) {
        tdcWidthExtended_->Fill(hits[ihit]->width());
        maxWidthInPlane = std::max<double>(hits[ihit]->width(), maxWidthInPlane);
        hocc_.fill(*hits[ihit]);
      }
    }

    tdcMaxWidthExtended_->Fill(maxWidthInPlane);
    // 56 - 4(pc) - 8(dense) = 44
    if((44 <= iplane) && (iplane <= 52)) {
      tdcMaxWidthDenseExtended_->Fill(maxWidthInPlane);
    }
  }

}

//================================================================
