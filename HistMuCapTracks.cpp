// Andrei Gaponenko, 2013

#include "HistMuCapTracks.h"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "EventClass.h"

//================================================================
void HistMuCapTracks::init(const std::string& hdir,
                           const std::string& namePrefix,
                           int chargeSign,
                           HistogramFactory &hf,
                           const ConfigFile& conf)
{
  cutCharge_ = chargeSign;
  cutTrackWinTimeDiff_ = conf.read<double>("MuCapture/ProtonWindow/cutTrackTimeDiff");
  cutTrackRmax_ = conf.read<double>("MuCapture/ProtonWindow/cutTrackRmax");
  cutTrackMuonOffset_ = conf.read<double>("MuCapture/ProtonWindow/cutTrackMuonOffset");

  trackwintimeLargeScale_ = hf.DefineTH1D(hdir, namePrefix+"trackWinTimeLS",
                                          namePrefix + "t(track) - t(win)",
                                          1100, -10000., 1000.);

  trackwintime_ = hf.DefineTH1D(hdir, namePrefix+"trackWinTime",
                                namePrefix + "t(track) - t(win)",
                                700, -200., 500.);

  hStartStop_ = hf.DefineTH2D(hdir,
                              namePrefix+"startStopPlane",
                              namePrefix + "startStopPlane",
                              56, 0.5, 56.5, 56, 0.5, 56.5);

  hStartStop_->SetOption("colz");

  trackz_ = hf.DefineTH1D(hdir, namePrefix+"trackz",
                          namePrefix + "trackz",
                          700, -89.5, 610.5);

  trackMuonOffset_ = hf.DefineTH2D(hdir,
                                   namePrefix+"trackMuonOffset",
                                   namePrefix + "track-muon UV",
                                   601, -3.05, 3.05, 601, -3.05, +3.05);

  trackMuonOffset_->SetOption("colz");

  trackMuondr_ = hf.DefineTH1D(hdir,
                               namePrefix+"trackMuondr",
                               namePrefix + "track-muon UV distance",
                               100, 0., 5.);

  costhVsPtot_ = hf.DefineTH2D(hdir,
                               namePrefix+"ptotVsCosth",
                               namePrefix + "ptot vs costh",
                               140, 0., 70., 100, -1., +1.);

  costhVsPtot_->SetOption("colz");

  u0v0_ = hf.DefineTH2D(hdir,
                        namePrefix+"u0v0",
                        namePrefix + "u0v0",
                        640, -16., 16, 640, -16., 16.);

  u0v0_->SetOption("colz");

  trackRL_ = hf.DefineTH2D(hdir,
                           namePrefix+"trackRL",
                           namePrefix + "track wavelengh vs radius",
                           200, 0., 20, 1600, -80., 80.);

  trackRL_->SetOption("colz");

  helixCenterUV_ = hf.DefineTH2D(hdir,
                                 namePrefix+"helixCenterUV",
                                 namePrefix + "helix center UV",
                                 320, -16., 16., 320, -16., 16.);

  helixCenterUV_->SetOption("colz");


  trackROut_ = hf.DefineTH1D(hdir, namePrefix+"rout",
                             namePrefix+"max track distance from the Z axis",
                             400, 0, 40.);

  hNumTracks_ = hf.DefineTH1D(hdir, namePrefix+"numSelectedTracks",
                              namePrefix+"numSelectedTracks",
                              10, -0.5, 9.5);
}

//================================================================
void HistMuCapTracks::fill(const EventClass& evt,
                           double timeWinStart,
                           double muStopU,
                           double muStopV)
{
  int numSelected(0);

  for(int i = 0; i < evt.ntr; ++i) {
    if(evt.hefit_ierror[i] == 0) {
      if(evt.hefit_q[i] == cutCharge_) {

        const double dt = evt.hefit_time[i] - timeWinStart;
        trackwintimeLargeScale_->Fill(dt);
        trackwintime_->Fill(dt);
        if(std::abs(dt) < cutTrackWinTimeDiff_) {

          trackRL_->Fill(evt.radius[i], evt.wavelen[i]);
          if(evt.radius[i] < cutTrackRmax_) {

            const double du = evt.hefit_u0[i] - muStopU;
            const double dv = evt.hefit_v0[i] - muStopV;
            const double dr = sqrt(std::pow(du,2)+std::pow(dv,2));
            trackMuonOffset_->Fill(du, dv);
            trackMuondr_->Fill(dr);
            if(dr < cutTrackMuonOffset_) {

              hStartStop_->Fill(evt.hefit_pstart[i], evt.hefit_pstop[i]);
              // Select tracks that go from tgt to the end of tracker
              if( (31== evt.hefit_pstart[i]) && (52 == evt.hefit_pstop[i])) {

                ++numSelected;
                trackz_->Fill(evt.hefit_z[i]);

                helixCenterUV_->Fill(evt.hefit_ucenter[i], evt.hefit_vcenter[i]);

                const double rout =
                  sqrt(std::pow(evt.hefit_ucenter[i], 2) + std::pow(evt.hefit_vcenter[i], 2))
                  + evt.radius[i];

                trackROut_->Fill(rout);

                costhVsPtot_->Fill(evt.ptot[i], evt.costh[i]);
                u0v0_->Fill(evt.hefit_u0[i], evt.hefit_v0[i]);

              } // pstart,pstop

            } // cutTrackMuonOffset

          } // cutTrackRmax
        } // cutTrackWinTimeDiff
      } // cutCharge
    } // ierror
  }

  hNumTracks_->Fill(numSelected);
}

//================================================================
