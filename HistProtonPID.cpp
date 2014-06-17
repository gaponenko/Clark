// Andrei Gaponenko, 2014

#include "HistProtonPID.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"
#include "EventClass.h"

//================================================================
void HistProtonPID::init(const std::string& hdir,
                         HistogramFactory& hf,
                         const ConfigFile& conf)
{
  hNumClusters78_ = hf.DefineTH2D(hdir, "numClusters78", "PC8 vs PC7 num clusters", 5, -0.5, 4.5, 5, -0.5, 4.5);
  hNumClusters78_->SetOption("colz");

  hSingleClusterSize78_ = hf.DefineTH2D(hdir, "signleClusterSize78", "PC8 vs PC7 cluster size for single cluster events", 5, -0.5, 4.5, 5, -0.5, 4.5);
  hSingleClusterSize78_->SetOption("colz");

  hsum78cos_vs_p_all_ = hf.DefineTH2D(hdir, "sum78cosall", "PC78 mean width cos(theta) vs p, all", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_all_->SetOption("colz");

  hsum78cos_vs_p_11_ = hf.DefineTH2D(hdir, "sum78cos11", "PC78 mean width cos(theta) vs p, 11", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_11_->SetOption("colz");

  hsum78cos_vs_p_12_ = hf.DefineTH2D(hdir, "sum78cos12", "PC78 mean width cos(theta) vs p, 12", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_12_->SetOption("colz");

  hsum78cos_vs_p_21_ = hf.DefineTH2D(hdir, "sum78cos21", "PC78 mean width cos(theta) vs p, 21", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_21_->SetOption("colz");

  hsum78cos_vs_p_22_ = hf.DefineTH2D(hdir, "sum78cos22", "PC78 mean width cos(theta) vs p, 22", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_22_->SetOption("colz");
}

//================================================================
void HistProtonPID::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters)
{
  hNumClusters78_->Fill(protonGlobalClusters[29].size(), protonGlobalClusters[30].size());

  const WireClusterCollection& pc7clusters = protonGlobalClusters[29];
  const WireClusterCollection& pc8clusters = protonGlobalClusters[30];

  if((pc7clusters.size() == 1) && (pc8clusters.size() == 1)) {

    hSingleClusterSize78_->Fill(pc7clusters[0].numCells(),
                                pc8clusters[0].numCells()
                                );

    // Use mean value per cell, not per hit
    double pc78meanWidth = (pc7clusters[0].totalTDCWidth() + pc8clusters[0].totalTDCWidth())
      /(pc7clusters[0].numCells() + pc8clusters[0].numCells());
    hsum78cos_vs_p_all_->Fill(evt.ptot[itrack], pc78meanWidth * evt.costh[itrack]);

    TH1 *hver = 0;
    if(pc7clusters[0].numCells() == 1) {
      if(pc8clusters[0].numCells() == 1) {
        hver = hsum78cos_vs_p_11_;
      }
      else if(pc8clusters[0].numCells() == 2) {
        hver = hsum78cos_vs_p_12_;
      }
    }
    else if(pc7clusters[0].numCells() == 2) {
      if(pc8clusters[0].numCells() == 1) {
        hver = hsum78cos_vs_p_21_;
      }
      else if(pc8clusters[0].numCells() == 2) {
        hver = hsum78cos_vs_p_22_;
      }
    }
    if(hver) {
      hver->Fill(evt.ptot[itrack], pc78meanWidth * evt.costh[itrack]);
    }

  }
}

//================================================================
