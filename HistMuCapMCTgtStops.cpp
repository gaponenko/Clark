// Andrei Gaponenko, 2014

#include "HistMuCapMCTgtStops.h"

#include <stdexcept>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapMCTgtStops::init(HistogramFactory& hf,
                               const std::string& hdir,
                               const DetectorGeo& geom,
                               const ConfigFile& conf)
{
  geom_ = &geom;

  hMcMuonTotalMultiplicity_ = hf.DefineTH1D(hdir, "mcMuonTotalMultiplicity", "MC mu- total multiplicity", 10, -0.5, 9.5);
  hMcMuonTrigCandidateMultiplicity_ = hf.DefineTH1D(hdir, "mcMuonTrigCandidateMultiplicity", "MC mu- trig candidate multiplicity", 10, -0.5, 9.5);

  hMcMuonStopTime_ = hf.DefineTH1D(hdir, "mcMuonStopTime", "MC muon stop time", 300, -50., 2950.);

  hMcMuonZStopFine_ = hf.DefineTH1D(hdir, "mcMuonZStopFine", "MC muon Z stop position", 120, -0.0200, -0.0080);
  hMcMuonZStopCoarse_ = hf.DefineTH1D(hdir, "mcMuonZStopCoarse", "MC muon Z stop position", 200, -1.0, +1.0);
}

//================================================================
void HistMuCapMCTgtStops::fill(const EventClass& evt) {
  hMcMuonTotalMultiplicity_->Fill(evt.mcMuonTotalMultiplicity);
  hMcMuonTrigCandidateMultiplicity_->Fill(evt.mcMuonTrigCandidateMultiplicity);

  if(evt.iMuStopMcVtxEnd != -1) {
    hMcMuonStopTime_->Fill(evt.mcvertex_time[evt.iMuStopMcVtxEnd]);
    const double zstop = evt.mcvertex_vz[evt.iMuStopMcVtxEnd];
    hMcMuonZStopFine_->Fill(zstop);
    hMcMuonZStopCoarse_->Fill(zstop);
  }
}

//================================================================
