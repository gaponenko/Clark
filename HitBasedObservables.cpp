// Andrei Gaponenko, 2014

#include "HitBasedObservables.h"

#include <sstream>
#include <iomanip>

#include "TH2.h"
#include "HistogramFactory.h"

//================================================================
void HistHitBasedAmbiguities::init(HistogramFactory& hf,
                                   const std::string& hdir,
                                   const DetectorGeo& geom,
                                   const ConfigFile& conf)
{
  hClusterMultiplicity_ = hf.DefineTH2D(hdir, "clusterMultiplicity", "Plane vs cluster multiplicity (in range, for accepted events)", 20, -0.5, 19.5, 56, 0.5, 56.5);
  hClusterMultiplicity_->SetOption("colz");
  hClusterMultiplicity_->GetXaxis()->SetTitle("num clusters");
  hClusterMultiplicity_->GetYaxis()->SetTitle("plane number");

  hDiffNWires_ = hf.DefineTH2D(hdir, "diffNWires", " Plane vs numWires(best) - numWires(other)", 161, -80.5, 80.5, 56, 0.5, 56.5);
  hDiffNWires_->SetOption("colz");
  hDiffNWires_->GetXaxis()->SetTitle("nbest - nother");
  hDiffNWires_->GetYaxis()->SetTitle("plane number");
}

//================================================================
template<class ClusterCmp>
HitBasedObservables<ClusterCmp>::HitBasedObservables(const ClustersByPlane& protonGlobalClusters, HistHitBasedAmbiguities *hh)
  : dnCPlanes_(0)
  , dnCWires_(0)
{
  // Start at PC7
  for(int iplane = 29; (iplane < protonGlobalClusters.size()) && !protonGlobalClusters.at(iplane).empty(); ++iplane) {
    const WireClusterCollection& clusters = protonGlobalClusters[iplane];

    if(hh) { hh->hClusterMultiplicity_->Fill(clusters.size(), iplane); }

    // Find the best cluster
    WireClusterCollection::const_iterator best = clusters.begin();
    for(WireClusterCollection::const_iterator current = ++clusters.begin(); current != clusters.end(); ++current) {
      if(cmp_(*best, *current)) {
        best = current;
      }
    }

    if(hh) {
      for(WireClusterCollection::const_iterator current = ++clusters.begin(); current != clusters.end(); ++current) {
        if(current != best) {
          hh->hDiffNWires_->Fill(best->numCells() - current->numCells(), iplane);
        }
      }
    }

    clusters_.push_back(*best);
  }

  dnCPlanes_ = clusters_.size();

  dnCWires_ = 0;
  for(int i=0; i<clusters_.size(); ++i) {
    dnCWires_ += clusters_[i].numCells();
  }
}

template HitBasedObservables<WireClusterCmpNumCells>::HitBasedObservables(const ClustersByPlane&, HistHitBasedAmbiguities*);
template HitBasedObservables<WireClusterCmpTDCWidth>::HitBasedObservables(const ClustersByPlane&, HistHitBasedAmbiguities*);
