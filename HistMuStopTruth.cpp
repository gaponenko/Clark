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
  hnumPrimaries_ = hf.DefineTH1D(hdir, "numPrimaries", "number of primary tracks", 10, -0.5, 9.5);
  hstopZ1_ = hf.DefineTH1D(hdir, "stopz1", "true Z stop", 1200, -60., 60.);
  hstopZ2_ = hf.DefineTH1D(hdir, "stopz2", "true Z stop", 200, -1., 1.);
  hstopZ3_ = hf.DefineTH1D(hdir, "stopz3", "true Z stop - Z tgt", 300, -0.03, 0.03);
}

//================================================================
void HistMuStopTruth::fill(const EventClass& evt) {
  hnumPrimaries_->Fill(evt.numPrimaryMcTrkCandidates);
  if(evt.iPrimaryMcTrk != -1) {
    const double zstop = evt.mcvertex_vz[evt.iPrimaryMcVtxEnd];
    hstopZ1_->Fill(zstop);
    hstopZ2_->Fill(zstop);
    hstopZ3_->Fill(zstop - zTargetCenter_);
  }
}

//================================================================
