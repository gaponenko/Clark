// Andrei Gaponenko, 2013

#include "HistMuStopTruth.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "TH2.h"

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuStopTruth::init(HistogramFactory &hf,
                           const std::string& hdir,
                           const DetectorGeo& geom,
                           const ConfigFile &conf)
{
  zTargetCenter_ = geom.zTargetCenter();
  targetHalfThickness_ = geom.targetThickness();

  hMcMuonTotalMultiplicity_ = hf.DefineTH1D(hdir, "mcMuonTotalMultiplicity", "MC mu- total multiplicity", 10, -0.5, 9.5);
  hMcMuonTrigCandidateMultiplicity_ = hf.DefineTH1D(hdir, "mcMuonTrigCandidateMultiplicity", "MC mu- trig candidate multiplicity", 10, -0.5, 9.5);

  hMcMuonStopTime_ = hf.DefineTH1D(hdir, "mcMuonStopTime", "MC muon stop time", 300, -50., 2950.);

  hstopZ1_ = hf.DefineTH1D(hdir, "stopz1", "true Z stop", 1200, -60., 60.);

  // Same binning as in_zstop in the reference sample handling
  hstopZ2_ = hf.DefineTH1D(hdir, "stopz2", "MC muon Z stop position", 120, -0.0200, -0.0080);

  // Align bin boundaries with the target.
  const int numBinsInTarget = 35; // in the 71um target
  const double binSize = geom.targetThickness()/numBinsInTarget;

  const double approximateHalfRange = 0.5; // cm, half histo scale
  const int numBins =  1 + 2*int(approximateHalfRange/binSize);
  const double fullRange = binSize * numBins;
  hstopdz_ = hf.DefineTH1D(hdir, "stopdz", "true Z stop - Z tgt", numBins, -0.5*fullRange, +0.5*fullRange);

  hpastTargetdz_ = hf.DefineTH1D(hdir, "pastTargetdz", "true Z stop - Z tgt for past target stops", numBins, -0.5*fullRange, +0.5*fullRange);
  hpastTargetUV_ = hf.DefineTH2D(hdir, "pastTargetUV", "true V vs U muon stop position",  200, -4.0, +4.0, 200, -4.0, +4.0);
  hpastTargetUV_->SetOption("colz");
}

//================================================================
void HistMuStopTruth::fill(const EventClass& evt) {
  hMcMuonTotalMultiplicity_->Fill(evt.mcMuonTotalMultiplicity);
  hMcMuonTrigCandidateMultiplicity_->Fill(evt.mcMuonTrigCandidateMultiplicity);

  if(evt.iMuStopMcVtxEnd != -1) {
    hMcMuonStopTime_->Fill(evt.mcvertex_time[evt.iMuStopMcVtxEnd]);
    const double zstop = evt.mcvertex_vz[evt.iMuStopMcVtxEnd];

    hstopZ1_->Fill(zstop);
    hstopZ2_->Fill(zstop);

    const double dz = zstop - zTargetCenter_;
    hstopdz_->Fill(dz);

    if(dz > targetHalfThickness_) {
      hpastTargetdz_->Fill(dz);
      hpastTargetUV_->Fill(evt.mcvertex_vu[evt.iMuStopMcVtxEnd], evt.mcvertex_vv[evt.iMuStopMcVtxEnd]);
    }
  }
}

//================================================================
