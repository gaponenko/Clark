// Andrei Gaponenko, 2013

#include "MuCapUVAnalysis.h"

#include <cmath>
#include <sstream>

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "EventClass.h"

#include "PlaneRange.h"

//================================================================
namespace {
  std::string formatTrack(const EventClass& evt, unsigned i) {
    std::ostringstream os;
    os<<"track["<<i<<"]: ptot="<<evt.ptot[i]
      <<", costh="<<evt.costh[i]
      <<", time="<<evt.hefit_time[i]
      <<", z="<<evt.hefit_z[i]
      <<", start="<<evt.hefit_pstart[i]
      <<", stop="<<evt.hefit_pstop[i]
      <<", u0="<<evt.hefit_u0[i]
      <<", v0="<<evt.hefit_v0[i]
      ;
    return os.str();
  }

  int distanceToTarget(int globalPlane) {
    return int(std::abs((globalPlane - 28.5)));
  }
}

//================================================================
void MuCapUVAnalysis::init(const std::string& hdir,
                           HistogramFactory &hf,
                           const ConfigFile& conf,
                           TimeWindow::StreamType cutStream,
                           double cutMinTime)
{
  cutCharge_ = -1;
  cutStream_ = cutStream;
  cutTrackMinTime_ = cutMinTime;
  cutTrackRmax_ = conf.read<double>("MuCapture/UVAnalysis/cutTrackRmax");
  cutCosThetaMin_ = conf.read<double>("MuCapture/UVAnalysis/cutCosThetaMin");
  cutCosThetaMax_ = conf.read<double>("MuCapture/UVAnalysis/cutCosThetaMax");
  cutPtMin_ = conf.read<double>("MuCapture/UVAnalysis/cutPtMin");
  cutPzMin_ = conf.read<double>("MuCapture/UVAnalysis/cutPzMin");
  cutPtotMin_ = conf.read<double>("MuCapture/UVAnalysis/cutPtotMin");
  cutPtotMax_ = conf.read<double>("MuCapture/UVAnalysis/cutPtotMax");
  cutTrackMuonOffset_ = conf.read<double>("MuCapture/UVAnalysis/cutTrackMuonOffset");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "cuts_r", "Tracks rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D(hdir, "cuts_p", "Tracks before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  //----------------------------------------------------------------
  hStartStop_ = hf.DefineTH2D(hdir,
                              "startStopPlane",
                              "Track stop vs start plane",
                              56, 0.5, 56.5, 56, 0.5, 56.5);

  hStartStop_->SetOption("colz");

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

  final_costhVsPtot_ = hf.DefineTH2D(hdir,
                                     "final_ptotVsCosth",
                                     "ptot vs costh",
                                     140, 0., 70., 100, -1., +1.);

  final_costhVsPtot_->SetOption("colz");

  final_u0v0_ = hf.DefineTH2D(hdir,
                              "final_u0v0",
                              "final_u0v0",
                              640, -16., 16, 640, -16., 16.);

  final_u0v0_->SetOption("colz");

  trackRL_ = hf.DefineTH2D(hdir,
                           "trackRL",
                           "track wavelengh vs radius",
                           200, 0., 20, 1600, -80., 80.);

  trackRL_->SetOption("colz");

  final_trackRL_ = hf.DefineTH2D(hdir,
                                 "final trackRL",
                                 "final track wavelengh vs radius",
                                 200, 0., 20, 1600, -80., 80.);

  final_trackRL_->SetOption("colz");

  helixCenterUV_ = hf.DefineTH2D(hdir,
                                 "helixCenterUV",
                                 "helix center UV",
                                 320, -16., 16., 320, -16., 16.);

  helixCenterUV_->SetOption("colz");


  final_trackROut_ = hf.DefineTH1D(hdir, "final_rout",
                                   "max track distance from the Z axis",
                                   400, 0, 40.);

  trackTime_ = hf.DefineTH1D(hdir, "trackTime", "track time",
                             1100, -1000., 10000.);

  final_trackTime_ = hf.DefineTH1D(hdir, "final_trackTime",
                                   "final track time",
                                   1100, -1000., 10000.);

  hNumTracks_ = hf.DefineTH1D(hdir, "numAcceptedTracks",
                              "numAcceptedTracks",
                              10, -0.5, 9.5);
}

//================================================================
int MuCapUVAnalysis::process(const EventClass& evt,
                             const ROOT::Math::XYPoint& muStopUV
                             )
{
  std::vector<unsigned> accepted;
  for(int i = 0; i < evt.ntr; ++i) {
    if(processTrack(i, evt, muStopUV)) {
      accepted.push_back(i);
    }
  }
  hNumTracks_->Fill(accepted.size());

  //// Look at the multiple-track events
  // if(accepted.size() > 1) {
  //   std::cout<<"Multiple accepted tracks in run "<<evt.nrun<<", event "<<evt.nevt<<std::endl;
  //   for(unsigned i=0; i < accepted.size(); ++i) {
  //     std::cout<<"\t"<<formatTrack(evt, accepted[i])<<std::endl;
  //   }
  // }

  int selected = accepted.empty() ? -1 : accepted[0];
  // find a track that starts closer to the target
  for(int i=1; i<accepted.size(); ++i) {
    if(distanceToTarget(evt.hefit_pstart[i]) < distanceToTarget(evt.hefit_pstart[selected])) {
      selected = i;
    }
  }

  return selected;
}

//================================================================
bool MuCapUVAnalysis::processTrack(int itrack,
                                   const EventClass& evt,
                                   const ROOT::Math::XYPoint& muStopUV
                                   )
{
  CutNumber c = analyzeTrack(itrack, evt, muStopUV);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c==CUTS_ACCEPTED;
}

//================================================================
MuCapUVAnalysis::CutNumber MuCapUVAnalysis::
analyzeTrack(int i, const EventClass& evt,
             const ROOT::Math::XYPoint& muStopUV
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
  if((cutStream_==TimeWindow::DOWNSTREAM)&&(evt.costh[i] < 0.) ||
     (cutStream_==TimeWindow::UPSTREAM)&&(evt.costh[i] > 0.) ) {
    return CUT_STREAM;
  }

  //----------------------------------------------------------------
  trackTime_->Fill(evt.hefit_time[i]);
  if(evt.hefit_time[i] < cutTrackMinTime_) {
    return CUT_TIME;
  }

  //----------------------------------------------------------------
  // Make sure the tracks are fully contained in the transverse direction
  trackRL_->Fill(evt.radius[i], evt.wavelen[i]);
  costhVsPtot_->Fill(evt.ptot[i], evt.costh[i]);
  if(evt.radius[i] > cutTrackRmax_) {
    return CUT_RADIUS;
  }

  //----------------------------------------------------------------
  // Other kinematic cuts
  if(std::abs(evt.costh[i]) < cutCosThetaMin_) {
    return CUT_COSTHETAMIN;
  }

  if(std::abs(evt.costh[i]) > cutCosThetaMax_) {
    return CUT_COSTHETAMAX;
  }

  if(evt.pt[i] < cutPtMin_) {
    return CUT_PTMIN;
  }

  if(std::abs(evt.pz[i]) < cutPzMin_) {
    return CUT_PZMIN;
  }

  if(evt.ptot[i] < cutPtotMin_) {
    return CUT_PTOTMIN;
  }

  if(evt.ptot[i] > cutPtotMax_) {
    return CUT_PTOTMAX;
  }

  //----------------------------------------------------------------
  const double du = evt.hefit_u0[i] - muStopUV.x();
  const double dv = evt.hefit_v0[i] - muStopUV.y();
  const double dr = sqrt(std::pow(du,2)+std::pow(dv,2));
  trackMuonOffset_->Fill(du, dv);
  trackMuondr_->Fill(dr);
  if(dr > cutTrackMuonOffset_) {
    return CUT_TRACK_MUON_OFFSET;
  }

  //----------------------------------------------------------------
  hStartStop_->Fill(evt.hefit_pstart[i], evt.hefit_pstop[i]);
  // Select tracks that go from tgt to the end of tracker
  // if( (31== evt.hefit_pstart[i]) && (52 == evt.hefit_pstop[i])) {}
  trackz_->Fill(evt.hefit_z[i]);

  helixCenterUV_->Fill(evt.hefit_ucenter[i], evt.hefit_vcenter[i]);

  const double rout =
    sqrt(std::pow(evt.hefit_ucenter[i], 2) + std::pow(evt.hefit_vcenter[i], 2))
    + evt.radius[i];

  final_trackROut_->Fill(rout);
  final_costhVsPtot_->Fill(evt.ptot[i], evt.costh[i]);
  final_trackRL_->Fill(evt.radius[i], evt.wavelen[i]);
  final_u0v0_->Fill(evt.hefit_u0[i], evt.hefit_v0[i]);
  final_trackTime_->Fill(evt.hefit_time[i]);

  //----------------------------------------------------------------
  return CUTS_ACCEPTED;
}
