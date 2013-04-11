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

  trackwintimeLargeScale_ = hf.DefineTH1D(hdir, namePrefix+"trackWinTimeLS",
                                          namePrefix + "t(track) - t(win)",
                                          1100, -10000., 1000.);

  trackwintime_ = hf.DefineTH1D(hdir, namePrefix+"trackWinTime",
                                namePrefix + "t(track) - t(win)",
                                700, -200., 500.);

  costhVsPtot_ = hf.DefineTH2D(hdir,
                               namePrefix+"ptotVsCosth",
                               namePrefix + "ptot vs costh",
                               140, 0., 70., 100, -1., +1.);

  costhVsPtot_->SetOption("colz");

  hStartStop_ = hf.DefineTH2D(hdir,
                              namePrefix+"startStopPlane",
                              namePrefix + "startStopPlane",
                              56, 0.5, 56.5, 56, 0.5, 56.5);

  hStartStop_->SetOption("colz");

  trackz_ = hf.DefineTH1D(hdir, namePrefix+"trackz",
                          namePrefix + "trackz",
                          700, -89.5, 610.5);

  u0v0_ = hf.DefineTH2D(hdir,
                        namePrefix+"u0v0",
                        namePrefix + "u0v0",
                        640, -160., 160, 640, -160., 160.);

  u0v0_->SetOption("colz");

  hNumTracks_ = hf.DefineTH1D(hdir, namePrefix+"numSelectedTracks",
                              namePrefix+"numSelectedTracks",
                              10, -0.5, 9.5);
}

//================================================================
void HistMuCapTracks::fill(const EventClass& evt, double timeWinStart) {
  int numSelected(0);

  for(int i = 0; i < evt.ntr; ++i) {
    if(evt.hefit_ierror[i] == 0) {
      if(evt.hefit_q[i] == cutCharge_) {
        const double dt = evt.hefit_time[i] - timeWinStart;
        trackwintimeLargeScale_->Fill(dt);
        trackwintime_->Fill(dt);
        if(std::abs(dt) < cutTrackWinTimeDiff_) {
          hStartStop_->Fill(evt.hefit_pstart[i], evt.hefit_pstop[i]);
          // Select tracks that go from tgt to the end of tracker
          if( (31== evt.hefit_pstart[i]) && (52 == evt.hefit_pstop[i])) {
            ++numSelected;
            trackz_->Fill(evt.hefit_z[i]);
            costhVsPtot_->Fill(evt.ptot[i], evt.costh[i]);
            u0v0_->Fill(evt.hefit_u0[i], evt.hefit_v0[i]);
          }
        }
      }
    }
  }

  hNumTracks_->Fill(numSelected);
}

//================================================================
