// Andrei Gaponenko, 2015

#include "HistRangeStudies.h"
#include <limits>
#include <cmath>
#include <algorithm>

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
  clusterMultiplicity_ = hf.DefineTProfile2D(hdir, "clusterMultiplicity", "cluster multiplicity  in track range profile", 88, 30., 250., 28, 28.5, 56.5);
  clusterMultiplicity_->SetOption("colz");
  clusterMultiplicity_->GetXaxis()->SetTitle("p [MeV/c]");
  clusterMultiplicity_->GetYaxis()->SetTitle("plane");

  //----------------------------------------------------------------
}

//================================================================
void HistRangeStudies::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
  const double ptot = evt.ptot[itrack];

  const int extended = MuCapUtilities::findExtendedLastPlane(evt, itrack, protonGlobalClusters);
  const int trange = MuCapUtilities::findTrackLastPlane(evt, itrack, protonGlobalClusters);

  extendedVsTrackRange_->Fill(trange-28, extended-28);
  rangeDiffVsPrec_->Fill(ptot, extended-trange);

  for(int iplane=29; iplane<=trange; ++iplane) {
    const unsigned numClusters= protonGlobalClusters.at(iplane).size();
    clusterMultiplicity_->Fill(ptot, iplane, numClusters);
  }
}

//================================================================
