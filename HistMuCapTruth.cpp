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

  helectrons_.init(hf, hdir+"/electrons", conf);

  hNumMCCaptureTracks_ = hf.DefineTH1D(hdir, "numCaptureMcTrkCandidates", "numCaptureMcTrkCandidates", 10, -0.5, 9.5);
  hCaptureTime_ = hf.DefineTH1D(hdir, "time", "MC time", 640, -100., 1500.);
  hptot_ = hf.DefineTH1D(hdir, "ptot", "MC ptot", 650, 0., 650.);
  hptot_proton_ = hf.DefineTH1D(hdir, "ptot_proton", "MC ptot proton", 650, 0., 650.);
  hptot_deuteron_ = hf.DefineTH1D(hdir, "ptot_deuteron", "MC ptot deuteron", 650, 0., 650.);
  hptot_triton_ = hf.DefineTH1D(hdir, "ptot_triton", "MC ptot triton", 650, 0., 650.);
  hptot_alpha_ = hf.DefineTH1D(hdir, "ptot_alpha", "MC ptot alpha", 650, 0., 650.);
  hek_ = hf.DefineTH1D(hdir, "ek", "MC Ek", 500, 0., 50.);
  hbeta_ = hf.DefineTH1D(hdir, "beta", "MC beta", 600, 0., 0.6);
  hphi_ = hf.DefineTH1D(hdir, "phi", "MC momentum phi", 100, -PI, +PI);
  hpcos_ = hf.DefineTH2D(hdir, "cosVsPtos", "MC cosTheta vs ptot", 650, 0., 650., 100, -1., 1.);
  hpcos_->SetOption("colz");

  hVUend_ = hf.DefineTH2D(hdir, "uvend", "V vs U primary end vtx for |Zend|<60 cm", 600, -30., 30., 600, -30., 30.);
  hVUend_->SetOption("box");

  hRZend_ = hf.DefineTH2D(hdir, "rzend", "Primary end vtx R vs Z", 1200, -60., 60., 300, 0., 30.);
  hRZend_->SetOption("colz");

  hRendVsPstart_ = hf.DefineTH2D(hdir, "rpstart", "End R vs start ptot for |Zend|<60 cm", 650, 0., 650., 300, 0., 30.);
  hRendVsPstart_->SetOption("colz");

  hZendVsPstart_ = hf.DefineTH2D(hdir, "zpstart", "End Z vs start ptot", 650, 0., 650., 1200, -60., 60.);
  hZendVsPstart_->SetOption("colz");

  // Same binning as G3 H4
  hZStart1_ = hf.DefineTH1D(hdir, "zstart1", "Capture track start Z, coarse", 100, -50., 50.);
  hZStop1_ = hf.DefineTH1D(hdir, "zstop1", "Muon stop Z, coarse", 100, -50., 50.);

  // Same bin size as G3 H5, but different range because of the target shift
  hZStart2_ = hf.DefineTH1D(hdir, "zstart2", "Capture track start Z", 300, -0.030, 0.030);
  hZStop2_ = hf.DefineTH1D(hdir, "zstop2", "Muon stop Z", 300, -0.030, 0.030);

  // Same binning as in_zstop in the reference sample handling
  hZStop3_ = hf.DefineTH1D(hdir, "zstop3", "Muon stop Z", 120, -0.0200, -0.0080);

  // dt_geo.00061: foil positions around the center: -0.41318 -0.00904 0.38671
  // foil thickness: 0.0005960 cm
  hZStop4_ = hf.DefineTH1D(hdir, "zstop4", "Muon stop Z", 1000, -0.5, +0.5);
}

//================================================================
void HistMuCapTruth::fill(const EventClass& evt) {
  helectrons_.fill(evt);
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
    case MuCapUtilities::PID_G3_TRITON:
      hptot_triton_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    case MuCapUtilities::PID_G3_ALPHA:
      hptot_alpha_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    default:
      std::ostringstream os;
      os<<"HistMuCapTruth: unknown capture PID="<<mcParticle;
      throw std::runtime_error(os.str());
    }

    const double Ek = MuCapUtilities::kineticEnergy(evt.mctrack_pid[imctrk], evt.mcvertex_ptot[imcvtxStart], evt);
    const double mass = MuCapUtilities::mass(evt.mctrack_pid[imctrk], evt);
    const double beta = evt.mcvertex_ptot[imcvtxStart] / (Ek + mass);

    hek_->Fill(Ek);
    hbeta_->Fill(beta);
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

  // Muon stop z is filled regardless of the presence of a capture track
  if(true) {
    // This code is stolen from HistMuCapAnalysisChannels::fillReferenceSample()
    // Look for the trigger muon stop
    int iMuStopTrack = -1;
    for(unsigned i=0; i<evt.nmctr; ++i) {
      if((evt.mctrack_pid[i] == MuCapUtilities::PID_G3_MUMINUS) ||
         (evt.mctrack_pid[i] == MuCapUtilities::PID_G4_MUMINUS)) {

        // Look at the end vertex of the muon track
        const int itmpvtxstart = evt.getFirstMCVertexIndexForTrack(i);
        const int itmpvtxend = itmpvtxstart + evt.mctrack_nv[i] - 1;
        const double stoptime = evt.mcvertex_time[itmpvtxend];
        if(std::abs(stoptime) < 100.) {
          if(iMuStopTrack == -1) {
            iMuStopTrack = i;
          }
        }
      }
    }
    if(iMuStopTrack != -1) {
      // Set vertex indexes
      const int iMuStopVtxStart = evt.getFirstMCVertexIndexForTrack(iMuStopTrack);
      const int iMuStopVtxEnd = iMuStopVtxStart + evt.mctrack_nv[iMuStopTrack] - 1;
      const double zstop = evt.mcvertex_vz[iMuStopVtxEnd];
      hZStop1_->Fill(zstop);
      hZStop2_->Fill(zstop);
      hZStop3_->Fill(zstop);
      hZStop4_->Fill(zstop);
    }
  }

}

//================================================================
