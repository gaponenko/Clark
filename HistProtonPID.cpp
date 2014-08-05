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
  pidPC7_.init(hdir+"/pidPC7", 22+7, hf, conf);
  pidPC8_.init(hdir+"/pidPC8", 22+8, hf, conf);
  pidDC23_.init(hdir+"/pidDC23", 8+23, hf, conf);
  pidDC24_.init(hdir+"/pidDC24", 8+24, hf, conf);

  //----------------------------------------------------------------
  hNumClusters78_ = hf.DefineTH2D(hdir, "numClusters78", "PC8 vs PC7 num clusters", 5, -0.5, 4.5, 5, -0.5, 4.5);
  hNumClusters78_->SetOption("colz");

  hSingleClusterSize78_ = hf.DefineTH2D(hdir, "singleClusterSize78", "PC8 vs PC7 cluster size for single cluster events", 5, -0.5, 4.5, 5, -0.5, 4.5);
  hSingleClusterSize78_->SetOption("colz");

  hsum78cos_vs_p_11_ = hf.DefineTH2D(hdir, "sum78cos11", "PC78 mean width cos(theta) vs p, 11", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_11_->SetOption("colz");

  hsum78cos_vs_p_12_ = hf.DefineTH2D(hdir, "sum78cos12", "PC78 mean width cos(theta) vs p, 12", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_12_->SetOption("colz");

  hsum78cos_vs_p_21_ = hf.DefineTH2D(hdir, "sum78cos21", "PC7 mean width cos(theta) vs p, 21", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_21_->SetOption("colz");

  hsum78cos_vs_p_22_ = hf.DefineTH2D(hdir, "sum78cos22", "PC7 mean width cos(theta) vs p, 22", 150, 0., 300., 120, 0, 600.);
  hsum78cos_vs_p_22_->SetOption("colz");

  //----------------------------------------------------------------
  // calibrated PC distributions
  hcs78cos_vs_p_11_ = hf.DefineTH2D(hdir, "calib78cos11", "PC78 calib width cos(theta) vs p, 11", 150, 0., 300., 60, -150., 150.);
  hcs78cos_vs_p_11_->SetOption("colz");

  hcs78cos_vs_p_12_ = hf.DefineTH2D(hdir, "calib78cos12", "PC78 calib width cos(theta) vs p, 12", 150, 0., 300., 60, -150., 150.);
  hcs78cos_vs_p_12_->SetOption("colz");

  hcs78cos_vs_p_21_ = hf.DefineTH2D(hdir, "calib78cos21", "PC7 calib width cos(theta) vs p, 21", 150, 0., 300., 60, -150., 150.);
  hcs78cos_vs_p_21_->SetOption("colz");

  hcs78cos_vs_p_22_ = hf.DefineTH2D(hdir, "calib78cos22", "PC7 calib width cos(theta) vs p, 22", 150, 0., 300., 60, -150., 150.);
  hcs78cos_vs_p_22_->SetOption("colz");

  //----------------------------------------------------------------
  // DC hists
  hNumClusters2324_ = hf.DefineTH2D(hdir, "numClusters2324", "DC24 vs DC23 num clusters", 5, -0.5, 4.5, 5, -0.5, 4.5);
  hNumClusters2324_->SetOption("colz");

  hSingleClusterSize2324_ = hf.DefineTH2D(hdir, "singleClusterSize2324", "DC24 vs DC23 cluster size for single cluster events", 5, -0.5, 4.5, 5, -0.5, 4.5);
  hSingleClusterSize2324_->SetOption("colz");

  hsum2324cos_vs_p_11_ = hf.DefineTH2D(hdir, "sum2324cos11", "DC2324 mean width cos(theta) vs p, 11", 200, 0., 300., 160, 0, 800.);
  hsum2324cos_vs_p_11_->SetOption("colz");

  hsum2324cos_vs_p_12_ = hf.DefineTH2D(hdir, "sum2324cos12", "DC2324 mean width cos(theta) vs p, 12", 200, 0., 300., 160, 0, 800.);
  hsum2324cos_vs_p_12_->SetOption("colz");

  hsum2324cos_vs_p_21_ = hf.DefineTH2D(hdir, "sum2324cos21", "DC23 mean width cos(theta) vs p, 21", 200, 0., 300., 160, 0, 800.);
  hsum2324cos_vs_p_21_->SetOption("colz");

  hsum2324cos_vs_p_22_ = hf.DefineTH2D(hdir, "sum2324cos22", "DC23 mean width cos(theta) vs p, 22", 200, 0., 300., 160, 0, 800.);
  hsum2324cos_vs_p_22_->SetOption("colz");
}

//================================================================
void HistProtonPID::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters)
{
  pidPC7_.fill(evt, itrack, protonGlobalClusters);
  pidPC8_.fill(evt, itrack, protonGlobalClusters);
  pidDC23_.fill(evt, itrack, protonGlobalClusters);
  pidDC24_.fill(evt, itrack, protonGlobalClusters);

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

    if(pc7clusters[0].numCells() == 1) {
      if(pc8clusters[0].numCells() == 1) {
        hsum78cos_vs_p_11_->Fill(evt.ptot[itrack], pc78meanWidth * evt.costh[itrack]);
        //hcs78cos_vs_p_11_->Fill(evt.ptot[itrack], getCalib(11, pc8meanWidth * evt.costh[itrack]));
      }
      else if(pc8clusters[0].numCells() == 2) {
        hsum78cos_vs_p_12_->Fill(evt.ptot[itrack], pc78meanWidth * evt.costh[itrack]);
      }
    }
    else if(pc7clusters[0].numCells() == 2) {
      if(pc8clusters[0].numCells() == 1) {
        hsum78cos_vs_p_21_->Fill(evt.ptot[itrack], pc78meanWidth * evt.costh[itrack]);
      }
      else if(pc8clusters[0].numCells() == 2) {
        hsum78cos_vs_p_22_->Fill(evt.ptot[itrack], pc78meanWidth * evt.costh[itrack]);
      }
    }
  }


  //----------------------------------------------------------------
  // DC variables
  hNumClusters2324_->Fill(protonGlobalClusters[31].size(), protonGlobalClusters[32].size());

  const WireClusterCollection& dc23clusters = protonGlobalClusters[31];
  const WireClusterCollection& dc24clusters = protonGlobalClusters[32];

  if((dc23clusters.size() == 1) && (dc24clusters.size() == 1)) {
    hSingleClusterSize2324_->Fill(dc23clusters[0].numCells(),
                                  dc23clusters[0].numCells()
                                  );

    // Use mean value per cell, not per hit
    double dc2324meanWidth = (dc23clusters[0].totalTDCWidth() + dc24clusters[0].totalTDCWidth())
      /(dc23clusters[0].numCells() + dc24clusters[0].numCells());

    if(dc23clusters[0].numCells() == 1) {
      if(dc24clusters[0].numCells() == 1) {
        hsum2324cos_vs_p_11_->Fill(evt.ptot[itrack], dc2324meanWidth * evt.costh[itrack]);
      }
      else if(dc24clusters[0].numCells() == 2) {
        hsum2324cos_vs_p_12_->Fill(evt.ptot[itrack], dc2324meanWidth * evt.costh[itrack]);
      }
    }
    else if(dc23clusters[0].numCells() == 2) {
      if(dc24clusters[0].numCells() == 1) {
        hsum2324cos_vs_p_21_->Fill(evt.ptot[itrack], dc2324meanWidth * evt.costh[itrack]);
      }
      else if(dc24clusters[0].numCells() == 2) {
        hsum2324cos_vs_p_22_->Fill(evt.ptot[itrack], dc2324meanWidth * evt.costh[itrack]);
      }
    }
  }

}

//================================================================
