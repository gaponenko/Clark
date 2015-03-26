// Andrei Gaponenko, 2015

#include "HistMuCapMuonRange.h"

#include "TH1.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapMuonRange::init(HistogramFactory& hf,
                              const std::string& hdir,
                              const DetectorGeo& geom,
                              const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  hMuonLastPlane_ = hf.DefineTH1D(hdir, "muonLastPlane", "Muon last plane", 56, 0.5, 56.5);

  hMuonRangeGaps_ = hf.DefineTH2D(hdir, "muonRangeGaps", "Muon range gap end vs start", 57, -0.5, 56.5, 57, -0.5, 56.5);
  hMuonRangeGaps_->SetOption("colz");

  hMuonMissingPlanes_ = hf.DefineTH1D(hdir, "muonMissingPlanes", "Muon missing planes", 56, 0.5, 56.5);

  if(doMCTruth_) {
    hTruePProtonVsLastPlane_ = hf.DefineTH2D(hdir, "muonLastPlaneVsTruePCapture_proton", "Muon last plane vs true p proton",
                                             56, 0.5, 56.5, 500, 0., 500.);
    hTruePProtonVsLastPlane_->SetOption("colz");
    hTruePProtonVsLastPlane_->GetXaxis()->SetTitle("Muon last plane");
    hTruePProtonVsLastPlane_->GetYaxis()->SetTitle("p proton, MeV/c");

    hTruePDeuteronVsLastPlane_ = hf.DefineTH2D(hdir, "muonLastPlaneVsTruePCapture_deuteron", "Muon last plane vs true p deuteron",
                                               56, 0.5, 56.5, 500, 0., 500.);
    hTruePDeuteronVsLastPlane_->SetOption("colz");
    hTruePDeuteronVsLastPlane_->GetXaxis()->SetTitle("Muon last plane");
    hTruePDeuteronVsLastPlane_->GetYaxis()->SetTitle("p deuteron, MeV/c");
  }

}

//================================================================
void HistMuCapMuonRange::fill(const PlaneRange& muonRange,
                              const ClustersByPlane& muonGlobalClusters,
                              const EventClass& evt)
{

  hMuonLastPlane_->Fill(muonRange.max());

  for(int i = 0; i + 1 < muonRange.segments().size(); ++i) {
    hMuonRangeGaps_->Fill(muonRange.segments()[i].max + 1, muonRange.segments()[i+1].min - 1);
  }
  for(int iplane = 1; iplane <= muonRange.max(); ++iplane) {
    if(muonGlobalClusters[iplane].empty()) {
      hMuonMissingPlanes_->Fill(iplane);
    }
  }

  //----------------
  if(doMCTruth_) {
    const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
    const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;
    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      hTruePProtonVsLastPlane_->Fill(muonRange.max(), evt.mcvertex_ptot[imcvtxStart]);
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      hTruePDeuteronVsLastPlane_->Fill(muonRange.max(), evt.mcvertex_ptot[imcvtxStart]);
      break;
    default:
      break;
    }
  }

}

//================================================================
