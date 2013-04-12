// Andrei Gaponenko, 2013

#include "MuCapUVAnalysis.h"

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "EventClass.h"

#include "PlaneRange.h"

//================================================================
void MuCapUVAnalysis::init(const std::string& hdir,
                           HistogramFactory &hf,
                           const ConfigFile& conf)
{
  cutCharge_ = -1;
  cutTrackWinTimeDiff_ = conf.read<double>("MuCapture/ProtonWindow/cutTrackTimeDiff");
  cutTrackRmax_ = conf.read<double>("MuCapture/ProtonWindow/cutTrackRmax");
  cutTrackMuonOffset_ = conf.read<double>("MuCapture/ProtonWindow/cutTrackMuonOffset");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D(hdir, "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  //----------------------------------------------------------------
  trackwintimeLargeScale_ = hf.DefineTH1D(hdir, "trackWinTimeLS",
                                          "t(track) - t(win)",
                                          1100, -10000., 1000.);

  trackwintime_ = hf.DefineTH1D(hdir, "trackWinTime",
                                "t(track) - t(win)",
                                700, -200., 500.);

  hStartStop_ = hf.DefineTH2D(hdir,
                              "startStopPlane",
                              "Track stop vs start plane",
                              56, 0.5, 56.5, 56, 0.5, 56.5);

  hStartStop_->SetOption("colz");

  hHitRange_ = hf.DefineTH2D(hdir,
                             "hitRange",
                             "time window hit range stop vs start",
                             56, 0.5, 56.5, 56, 0.5, 56.5);

  hHitRange_->SetOption("colz");

  trackz_ = hf.DefineTH1D(hdir, "trackz",
                          "trackz",
                          700, -89.5, 610.5);

  trackMuonOffset_ = hf.DefineTH2D(hdir,
                                   "trackMuonOffset",
                                   "track-muon UV",
                                   601, -3.05, 3.05, 601, -3.05, +3.05);

  trackMuonOffset_->SetOption("colz");

  trackMuondr_ = hf.DefineTH1D(hdir,
                               "trackMuondr",
                               "track-muon UV distance",
                               100, 0., 5.);

  costhVsPtot_ = hf.DefineTH2D(hdir,
                               "ptotVsCosth",
                               "ptot vs costh",
                               140, 0., 70., 100, -1., +1.);

  costhVsPtot_->SetOption("colz");

  u0v0_ = hf.DefineTH2D(hdir,
                        "u0v0",
                        "u0v0",
                        640, -16., 16, 640, -16., 16.);

  u0v0_->SetOption("colz");

  trackRL_ = hf.DefineTH2D(hdir,
                           "trackRL",
                           "track wavelengh vs radius",
                           200, 0., 20, 1600, -80., 80.);

  trackRL_->SetOption("colz");

  helixCenterUV_ = hf.DefineTH2D(hdir,
                                 "helixCenterUV",
                                 "helix center UV",
                                 320, -16., 16., 320, -16., 16.);

  helixCenterUV_->SetOption("colz");


  trackROut_ = hf.DefineTH1D(hdir, "rout",
                             "max track distance from the Z axis",
                             400, 0, 40.);

  hNumTracks_ = hf.DefineTH1D(hdir, "numSelectedTracks",
                              "numSelectedTracks",
                              10, -0.5, 9.5);
}

//================================================================
void MuCapUVAnalysis::process(const EventClass& evt,
                              double timeWinStart,
                              const ClustersByPlane& globalClusters,
                              double muStopU, double muStopV
                              )
{
  int numSelected(0);
  for(int i = 0; i < evt.ntr; ++i) {
    if(processTrack(i, evt, timeWinStart, globalClusters, muStopU, muStopV)) {
      ++numSelected;
    }
  }
  hNumTracks_->Fill(numSelected);
}

//================================================================
bool MuCapUVAnalysis::processTrack(int itrack,
                                   const EventClass& evt,
                                   double timeWinStart,
                                   const ClustersByPlane& globalClusters,
                                   double muStopU, double muStopV
                                   )
{
  CutNumber c = analyzeTrack(itrack, evt, timeWinStart, globalClusters, muStopU, muStopV);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c==CUTS_ACCEPTED;
}

//================================================================
MuCapUVAnalysis::CutNumber MuCapUVAnalysis::
analyzeTrack(int i, const EventClass& evt,
             double timeWinStart,
             const ClustersByPlane& globalClusters,
             double muStopU, double muStopV
             )
{
  //----------------------------------------------------------------
  if(evt.hefit_ierror[i] != 0) {
    return CUT_IERROR;
  }

  //----------------------------------------------------------------
  if(evt.hefit_q[i] != cutCharge_) {
    return CUT_CHARGE;
  }

  //----------------------------------------------------------------
  const double dt = evt.hefit_time[i] - timeWinStart;
  trackwintimeLargeScale_->Fill(dt);
  trackwintime_->Fill(dt);
  if(std::abs(dt) > cutTrackWinTimeDiff_) {
    return CUT_TIME;
  }

  //----------------------------------------------------------------
  trackRL_->Fill(evt.radius[i], evt.wavelen[i]);
  if(evt.radius[i] > cutTrackRmax_) {
    return CUT_RADIUS;
  }

  //----------------------------------------------------------------
  const double du = evt.hefit_u0[i] - muStopU;
  const double dv = evt.hefit_v0[i] - muStopV;
  const double dr = sqrt(std::pow(du,2)+std::pow(dv,2));
  trackMuonOffset_->Fill(du, dv);
  trackMuondr_->Fill(dr);
  if(dr > cutTrackMuonOffset_) {
    return CUT_TRACK_MUON_OFFSET;
  }

  //----------------------------------------------------------------
  const PlaneRange gr = findPlaneRange(globalClusters);
  hHitRange_->Fill(gr.min, gr.max);

  hStartStop_->Fill(evt.hefit_pstart[i], evt.hefit_pstop[i]);
  // Select tracks that go from tgt to the end of tracker
  // if( (31== evt.hefit_pstart[i]) && (52 == evt.hefit_pstop[i])) {}

  trackz_->Fill(evt.hefit_z[i]);

  helixCenterUV_->Fill(evt.hefit_ucenter[i], evt.hefit_vcenter[i]);

  const double rout =
    sqrt(std::pow(evt.hefit_ucenter[i], 2) + std::pow(evt.hefit_vcenter[i], 2))
    + evt.radius[i];

  trackROut_->Fill(rout);

  costhVsPtot_->Fill(evt.ptot[i], evt.costh[i]);
  u0v0_->Fill(evt.hefit_u0[i], evt.hefit_v0[i]);

  //----------------------------------------------------------------
  return CUTS_ACCEPTED;
}
