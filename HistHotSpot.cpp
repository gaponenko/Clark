// Andrei Gaponenko, 2014

#include "HistHotSpot.h"
#include <limits>
#include <cmath>
#include <algorithm>

#include "HitBasedObservables.h"
#include "EventClass.h"

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

#include "TimeWindow.h"

extern TimeWindowingResults *gwres;

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

//================================================================
void HistHotSpot::init(HistogramFactory& hf,
                       const std::string& hdir,
                       const DetectorGeo& geom,
                       const ConfigFile& conf)
{
  hitTimeWidest_ = hf.DefineTH1D(hdir, "hitTimeWidest", "Time of the widest hit (78)", 400, 0., 10000.);

  hitTimeAll78_ = hf.DefineTH1D(hdir, "hitTimeAll78", "Time of PC7,8 hits", 400, 0., 10000.);

  hitTimeAllPlane_.resize(2);
  hitTimeAllPlane_[0] = hf.DefineTH1D(hdir, "hitTimeAll7", "Time of PC7 hits", 400, 0., 10000.);
  hitTimeAllPlane_[1] = hf.DefineTH1D(hdir, "hitTimeAll8", "Time of PC8 hits", 400, 0., 10000.);

  //----------------------------------------------------------------
  tdcWidthWidestHit8vs7_ = hf.DefineTH2D(hdir, "widestWidth87", "TDC width 8 vs 7, widest hits", 50, 0., 1000., 50, 0., 1000.);

  tdcWidthWidestHit8vs7_->SetOption("colz");
  tdcWidthWidestHit8vs7_->GetXaxis()->SetTitle("PC7 hit width");
  tdcWidthWidestHit8vs7_->GetYaxis()->SetTitle("PC8 hit width");

  //----------------------------------------------------------------
  posWidestHit8vs7_ = hf.DefineTH2D(hdir, "widestPosition8vs7", "TDC widest hit wire, PC8 vs PC7",
                                    52, 54.5, 106.5, 52, 54.5, 106.5);

  posWidestHit8vs7_->SetOption("colz");
  posWidestHit8vs7_->GetXaxis()->SetTitle("PC7 hit width");
  posWidestHit8vs7_->GetYaxis()->SetTitle("PC8 hit width");

  //----------------------------------------------------------------
  clusterSize8vs7_ = hf.DefineTH2D(hdir, "clusterSize8vs7", "TDC cluster size, PC8 vs PC7",
                                    50, 0.5, 50.5, 50, 0.5, 50.5);

  clusterSize8vs7_->SetOption("colz");
  clusterSize8vs7_->GetXaxis()->SetTitle("PC7 hit width");
  clusterSize8vs7_->GetYaxis()->SetTitle("PC8 hit width");

  //----------------------------------------------------------------
  numUnassignedDCHits_ = hf.DefineTH1D(hdir, "numUnassignedDCHits", "Number of unassigned DC hits", 80, -0.5, 79.5);

  tdcWidthAll78_ = hf.DefineTH1D(hdir, "tdcWidthAll78", "TDC width PC7,8", 1000, 0., 1000.);

  tdcWidthAllPlane_.resize(2);
  tdcWidthAllPlane_[0] = hf.DefineTH1D(hdir, "tdcWidthAll7", "TDC width PC7", 1000, 0., 1000.);
  tdcWidthAllPlane_[1] = hf.DefineTH1D(hdir, "tdcWidthAll8", "TDC width PC8", 1000, 0., 1000.);

  dt_.resize(2);
  dt_[0] = hf.DefineTH1D(hdir, "dt7", "t(hit)-t(earliest hit), PC7", 200, 0., 200.);
  dt_[1] = hf.DefineTH1D(hdir, "dt8", "t(hit)-t(earliest hit), PC8", 200, 0., 200.);

  //----------------------------------------------------------------
}

//================================================================
void HistHotSpot::fill(const EventClass& evt, const ClustersByPlane& protonGlobalClusters) {
  HitBasedObservables obs(protonGlobalClusters);

  // Require hits in both PC7 and PC8
  if(obs.dnCPlanes() < 2) {
    return;
  }

  TDCHitWPPtr hit7 = maxTDCWidthHit(protonGlobalClusters.at(29));
  TDCHitWPPtr hit8 = maxTDCWidthHit(protonGlobalClusters.at(30));
  TDCHitWPPtr hit78 = (hit7->width() < hit8->width()) ? hit8 : hit7;

  hitTimeWidest_->Fill(hit78->time());
  tdcWidthWidestHit8vs7_->Fill(hit7->width(), hit8->width());
  posWidestHit8vs7_->Fill(hit7->cell(), hit8->cell());
  clusterSize8vs7_->Fill(obs.clusterSize().at(0), obs.clusterSize().at(1));

  numUnassignedDCHits_->Fill(gwres->unassignedDCHits.size());

  for(int iplane=29; iplane<=30; ++iplane) {

    // per plane
    TDCHitWPPtr earliestHit;

    for(int icluster=0; icluster<protonGlobalClusters.at(iplane).size(); ++icluster) {
      for(int ihit = 0; ihit<protonGlobalClusters[iplane].at(icluster).hits().size(); ++ihit) {
        TDCHitWPPtr hit = protonGlobalClusters[iplane][icluster].hits()[ihit];
        hitTimeAll78_->Fill(hit->time());
        hitTimeAllPlane_.at(iplane - 29)->Fill(hit->time());

        tdcWidthAll78_->Fill(hit->width());
        tdcWidthAllPlane_.at(iplane - 29)->Fill(hit->width());

        if((earliestHit == TDCHitWPPtr()) || (hit->time() < earliestHit->time())) {
          earliestHit = hit;
        }
      }
    }

    // Fill per-plane dt histgrams
    for(int icluster=0; icluster<protonGlobalClusters.at(iplane).size(); ++icluster) {
      for(int ihit = 0; ihit<protonGlobalClusters[iplane].at(icluster).hits().size(); ++ihit) {
        TDCHitWPPtr hit = protonGlobalClusters[iplane][icluster].hits()[ihit];
        if(hit != earliestHit) {
          dt_.at(iplane-29)->Fill(hit->time() - earliestHit->time());
        }
      }
    }
  }

}

//================================================================
