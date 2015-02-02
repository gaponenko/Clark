// Andrei Gaponenko, 2013

#include "HistMuCapTruth.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapTruth::init(HistogramFactory &hf,
                          const std::string& hdir,
                          const ConfigFile &conf)
{
  static const double PI = 4.*atan(1.);

  hNumMCCaptureTracks_ = hf.DefineTH1D(hdir, "numCaptureMcTrkCandidates", "numCaptureMcTrkCandidates", 10, -0.5, 9.5);
  hCaptureTime_ = hf.DefineTH1D(hdir, "time", "MC time", 640, -100., 1500.);
  hptot_ = hf.DefineTH1D(hdir, "ptot", "MC ptot", 500, 0., 500.);
  hptot_proton_ = hf.DefineTH1D(hdir, "ptot_proton", "MC ptot proton", 500, 0., 500.);
  hptot_deuteron_ = hf.DefineTH1D(hdir, "ptot_deuteron", "MC ptot deuteron", 500, 0., 500.);
  hek_ = hf.DefineTH1D(hdir, "ek", "MC Ek", 500, 0., 50.);
  hphi_ = hf.DefineTH1D(hdir, "phi", "MC momentum phi", 100, -PI, +PI);
  hpcos_ = hf.DefineTH2D(hdir, "cosVsPtos", "MC cosTheta vs ptot", 500, 0., 500., 100, -1., 1.);
  hpcos_->SetOption("colz");

  hVUend_ = hf.DefineTH2D(hdir, "uvend", "V vs U primary end vtx for |Zend|<60 cm", 600, -30., 30., 600, -30., 30.);
  hVUend_->SetOption("box");

  hRZend_ = hf.DefineTH2D(hdir, "rzend", "Primary end vtx R vs Z", 1200, -60., 60., 300, 0., 30.);
  hRZend_->SetOption("colz");

  hRendVsPstart_ = hf.DefineTH2D(hdir, "rpstart", "End R vs start ptot for |Zend|<60 cm", 500, 0., 500., 300, 0., 30.);
  hRendVsPstart_->SetOption("colz");

  hZendVsPstart_ = hf.DefineTH2D(hdir, "zpstart", "End Z vs start ptot", 500, 0., 500., 1200, -60., 60.);
  hZendVsPstart_->SetOption("colz");

  // Same binning as G3 H4
  hZStart1_ = hf.DefineTH1D(hdir, "zstart1", "Capture track start Z, coarse", 100, -50., 50.);

  // Same bin size as G3 H5, but different range because of the target shift
  hZStart2_ = hf.DefineTH1D(hdir, "zstart2", "Capture track start Z", 300, -0.030, 0.030);
}

//================================================================
void HistMuCapTruth::fill(const EventClass& evt) {
  hNumMCCaptureTracks_->Fill(evt.numCaptureMcTrkCandidates);
  if(evt.iCaptureMcTrk != -1) {
    const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
    const unsigned imcvtxEnd = evt.iCaptureMcVtxEnd;
    const unsigned imctrk = evt.iCaptureMcTrk;

    hCaptureTime_->Fill(evt.mcvertex_time[imcvtxStart]);

    hptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
    const int mcParticle = evt.mctrack_pid[imctrk];
    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      hptot_proton_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      hptot_deuteron_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    }

    hek_->Fill(MuCapUtilities::kineticEnergy(evt.mctrack_pid[imctrk], evt.mcvertex_ptot[imcvtxStart], evt));
    hphi_->Fill(evt.mcvertex_phimuv[imcvtxStart]);
    hpcos_->Fill(evt.mcvertex_ptot[imcvtxStart], evt.mcvertex_costh[imcvtxStart]);

    const double endR = std::sqrt(std::pow(evt.mcvertex_vu[imcvtxEnd], 2) + std::pow(evt.mcvertex_vv[imcvtxEnd], 2));
    hRZend_->Fill(evt.mcvertex_vz[imcvtxEnd], endR);

    if(std::abs(evt.mcvertex_vz[imcvtxEnd]) < 60.) {
      hVUend_->Fill(evt.mcvertex_vu[imcvtxEnd], evt.mcvertex_vv[imcvtxEnd]);
      hRendVsPstart_->Fill(evt.mcvertex_ptot[imcvtxStart], endR);
    }

    hZendVsPstart_->Fill(evt.mcvertex_ptot[imcvtxStart], evt.mcvertex_vz[imcvtxEnd]);

    hZStart1_->Fill(evt.mcvertex_vz[imcvtxStart]);
    hZStart2_->Fill(evt.mcvertex_vz[imcvtxStart]);
  }
}

//================================================================
