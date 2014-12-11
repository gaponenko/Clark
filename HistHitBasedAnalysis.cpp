// Andrei Gaponenko, 2014

#include "HistHitBasedAnalysis.h"
#include "HitBasedObservables.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistHitBasedAnalysis::init(HistogramFactory& hf,
                                const std::string& hdir,
                                const DetectorGeo& geom,
                                const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  const int recoCWiresNBins = 200;
  const double recoCWiresXMin = -0.5;
  const double recoCWiresXMax = 199.5;
  const int recoCPlanesNBins = geom.numGlobal()/2;
  const double recoCPlanesXMin = 0.5;
  const double recoCPlanesXMax =recoCPlanesNBins + 0.5;

  //----------------------------------------------------------------
  lastconPlaneVsCWires_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires", "Last contiguous plane vs sum(largest cluster size)",
                                        recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

  lastconPlaneVsCWires_->SetOption("colz");
  lastconPlaneVsCWires_->GetXaxis()->SetTitle("cluster wires");
  lastconPlaneVsCWires_->GetYaxis()->SetTitle("cplane");

  if(doMCTruth_) {
    lastconPlaneVsCWires_mcproton_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires_mcproton", "Last contiguous plane vs sum(largest cluster size), mcproton",
                                                   recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);
    lastconPlaneVsCWires_mcproton_->SetOption("colz");
    lastconPlaneVsCWires_mcproton_->GetXaxis()->SetTitle("cluster wires");
    lastconPlaneVsCWires_mcproton_->GetYaxis()->SetTitle("cplane");

    lastconPlaneVsCWires_mcdeuteron_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires_mcdeuteron", "Last contiguous plane vs sum(largest cluster size), mcdeuteron",
                                                   recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);
    lastconPlaneVsCWires_mcdeuteron_->SetOption("colz");
    lastconPlaneVsCWires_mcdeuteron_->GetXaxis()->SetTitle("cluster wires");
    lastconPlaneVsCWires_mcdeuteron_->GetYaxis()->SetTitle("cplane");

    lastconPlaneVsCWires_mcdio_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires_mcdio", "Last contiguous plane vs sum(largest cluster size), mcdio",
                                                   recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);
    lastconPlaneVsCWires_mcdio_->SetOption("colz");
    lastconPlaneVsCWires_mcdio_->GetXaxis()->SetTitle("cluster wires");
    lastconPlaneVsCWires_mcdio_->GetYaxis()->SetTitle("cplane");
  }

  noncontiguous_.init(hf, hdir+"/noncontiguous", geom, conf);

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning
    const int gen1nbins = 400; // 2.5 MeV/c bins
    const double gen1pmin = 0.;
    const double gen1pmax = 400.;

    // Migration matrices for the contained channel, two generator binnings
    migration_ = hf.DefineTH3D(hdir, "migration", "Hit based channel migration",
                               gen1nbins, gen1pmin, gen1pmax,
                               recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                               recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    migration_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_->GetYaxis()->SetTitle("cluster wires");
    migration_->GetZaxis()->SetTitle("cplane");

    hTruth_in_.init(hf, hdir+"/MCTruth_in", conf);
    hTruth_accepted_.init(hf, hdir+"/MCTruth_accepted", conf);
  }

  hOuterVetoNumHitPlanes_ = hf.DefineTH1D(hdir, "outerVetoNumHitPlanes", "outerVetoNumHitPlanes", 6, -0.5, 5.5);
  hNumPC7Clusters_ = hf.DefineTH1D(hdir, "numPC7Clusters", "numPC7Clusters", 20, -0.5, 19.5);
}

//================================================================
bool HistHitBasedAnalysis::accepted(const EventClass& evt, const ClustersByPlane& protonGlobalClusters) {
  if(doMCTruth_) {
    hTruth_in_.fill(evt);
  }

  // Z containment cut
  int numOuterVetoHitPlanes(0);
  for(int i=0; i<4; ++i) {
    if(!protonGlobalClusters[56-i].empty()) {
      ++numOuterVetoHitPlanes;
    }
  }

  hOuterVetoNumHitPlanes_->Fill(numOuterVetoHitPlanes);
  if(numOuterVetoHitPlanes > 1) {
    return false;
  }

  const unsigned numPC7Clusters = protonGlobalClusters.at(29).size();
  hNumPC7Clusters_->Fill(numPC7Clusters);
  if(!numPC7Clusters) {
    return false;
  }

  //----------------------------------------------------------------
  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  // Simulated DIO have no easily accessible MC truth.  We'll tread PID=zero as DIO down in this code.
  const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

  HitBasedObservables obs(protonGlobalClusters);

  lastconPlaneVsCWires_->Fill(obs.dnCWires(), obs.dnCPlanes());
  noncontiguous_.fill(protonGlobalClusters, evt);

  if(doMCTruth_) {
    if(imcvtxStart  != -1) {
      migration_->Fill(evt.mcvertex_ptot[imcvtxStart], obs.dnCWires(), obs.dnCPlanes());
    }

    hTruth_accepted_.fill(evt);

    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      lastconPlaneVsCWires_mcproton_->Fill(obs.dnCWires(), obs.dnCPlanes());
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      lastconPlaneVsCWires_mcdeuteron_->Fill(obs.dnCWires(), obs.dnCPlanes());
      break;
    case 0:
      lastconPlaneVsCWires_mcdio_->Fill(obs.dnCWires(), obs.dnCPlanes());
      break;
    }
  }

  return true;
}

//================================================================
