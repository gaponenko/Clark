// Andrei Gaponenko, 2013

#include "MuCapTrkAnalysisHF.h"

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

  double distanceToMuon(const EventClass& evt, int itrack, const ROOT::Math::XYPoint& muStopUV) {
    const double du = evt.hefit_u0[itrack] - muStopUV.x();
    const double dv = evt.hefit_v0[itrack] - muStopUV.y();
    const double dr = sqrt(std::pow(du,2)+std::pow(dv,2));
    return dr;
  }

}

//================================================================
void MuCapTrkAnalysisHF::init(const std::string& hdir,
                              HistogramFactory &hf,
                              const ConfigFile& conf,
                              int cutCharge,
                              TimeWindow::StreamType cutStream,
                              RecoResMuCapTrk *result)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");
  result_ = result;
  cutCharge_ = cutCharge;
  cutStream_ = cutStream;
  cutTrackWinTimedt_ = conf.read<int>("MuCapture/TrkAnalysisHF/cutTrackWinTimedt");
  cutTrackRmax_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutTrackRmax");
  cutCosThetaMin_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutCosThetaMin");
  cutCosThetaMax_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutCosThetaMax");
  cutPtMin_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutPtMin");
  cutPzMin_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutPzMin");
  cutPtotMin_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutPtotMin");
  cutPtotMax_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutPtotMax");
  // cutChi2_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutChi2");
  cutTrackMuonOffset_ = conf.read<double>("MuCapture/TrkAnalysisHF/cutTrackMuonOffset");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "cuts_r", "Tracks rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());
  h_cuts_r->SetOption("hist text");

  h_cuts_p = hf.DefineTH1D(hdir, "cuts_p", "Tracks before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);
  h_cuts_p->SetOption("hist text");

  //----------------------------------------------------------------
  trackWinTime_ = hf.DefineTH1D(hdir, "trackWinTime", "track time - win time",
                                201, -100.5, 100.5);

  hCharge_ = hf.DefineTH1D(hdir, "charge", "track charge", 3, -1.5, 1.5);

  hStartStop_ = hf.DefineTH2D(hdir,
                              "startStopPlane",
                              "Track stop vs start plane",
                              56, 0.5, 56.5, 56, 0.5, 56.5);

  hStartStop_->SetOption("colz");

  final_StartStop_ = hf.DefineTH2D(hdir,
                                   "final_startStopPlane",
                                   "Final track stop vs start plane",
                                   56, 0.5, 56.5, 56, 0.5, 56.5);

  final_StartStop_->SetOption("colz");

  final_trackz_ = hf.DefineTH1D(hdir, "trackz",
                                "trackz",
                                700, -89.5, 610.5);

  hTrackQuality_.init(hf, hdir+"/TrackQuality", conf);

  trackChi2_ = hf.DefineTH1D(hdir, "trackchi2", "trackchi2", 1000, 0, 500.);

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
                               300, 0., 300., 100, -1., +1.);

  costhVsPtot_->SetOption("colz");

  final_costhVsPtot_ = hf.DefineTH2D(hdir,
                                     "final_ptotVsCosth",
                                     "final ptot vs costh",
                                     300, 0., 300., 100, -1., +1.);

  final_costhVsPtot_->SetOption("colz");

  final_u0v0_ = hf.DefineTH2D(hdir,
                              "final_u0v0",
                              "final_u0v0",
                              640, -16., 16, 640, -16., 16.);

  final_u0v0_->SetOption("colz");

  trackRL_ = hf.DefineTH2D(hdir,
                           "trackRL",
                           "track wavelengh vs radius",
                           200, 0., 20, 800, -400., 400.);

  trackRL_->SetOption("colz");

  final_trackRL_ = hf.DefineTH2D(hdir,
                                 "final_trackRL",
                                 "final track wavelengh vs radius",
                                 200, 0., 20, 800, -400., 400.);

  final_trackRL_->SetOption("colz");

  final_helixCenterUV_ = hf.DefineTH2D(hdir,
                                       "final_helixCenterUV",
                                       "final helix center UV",
                                       320, -16., 16., 320, -16., 16.);

  final_helixCenterUV_->SetOption("colz");


  final_trackROut_ = hf.DefineTH1D(hdir, "final_rout",
                                   "max track distance from the Z axis",
                                   400, 0, 40.);

  trackTime_ = hf.DefineTH1D(hdir, "trackTime", "track time",
                             1100, -1000., 10000.);

  final_trackTime_ = hf.DefineTH1D(hdir, "final_trackTime",
                                   "final track time",
                                   1100, -1000., 10000.);

  final_trackWinTime_ = hf.DefineTH1D(hdir, "final_trackWinTime", "final track time - win time",
                                      201, -100.5, 100.5);

  hNumTracks_ = hf.DefineTH1D(hdir, "numAcceptedTracks",
                              "numAcceptedTracks",
                              10, -0.5, 9.5);

  hSelectedTrackQuality_.init(hf, hdir+"/SelectedTrackQuality", conf);

  hPerEventMomentum_ = hf.DefineTH1D(hdir, "perEventMomentum",
                                     "selected track momentum",
                                     300, 0., 300.);

  hSelectorPlane_ = hf.DefineTH2D(hdir, "selectorPlane", "track to target dplane max vs min",
                                  30, -0.5, 29.5, 30, -0.5, 29.5);

  hSelectorPlane_->SetOption("colz");

  hSelectordrmu_ = hf.DefineTH2D(hdir, "selectordrmu", "drmu max vs min",
                                 50, 0., 5., 50, 0., 5.);

  hSelectordrmu_->SetOption("colz");

  if(doMCTruth_) {
    htruthAccepted_.init(hf, hdir+"/truthAccepted", conf);
  }
}

//================================================================
int MuCapTrkAnalysisHF::process(const EventClass& evt,
                                const ROOT::Math::XYPoint& muStopUV,
                                const TimeWindow& protonWin)
{
  std::vector<unsigned> accepted;
  for(int i = 0; i < evt.ntr; ++i) {
    if(processTrack(i, evt, muStopUV, protonWin)) {
      accepted.push_back(i);
    }
  }
  hNumTracks_->Fill(accepted.size());

  // Look at the multiple-track events
  if(false && (accepted.size() > 1)) {
    std::cout<<"Multiple accepted tracks in run "<<evt.nrun<<", event "<<evt.nevt<<std::endl;
    for(unsigned i=0; i < accepted.size(); ++i) {
      std::cout<<"\t"<<formatTrack(evt, accepted[i])<<std::endl;
    }
  }

  int selected = accepted.empty() ? -1 : accepted[0];
  // find a track that starts closer to the target
  for(int i=1; i<accepted.size(); ++i) {
    const int distBest = distanceToTarget(evt.hefit_pstart[selected]);
    const int distCurrent = distanceToTarget(evt.hefit_pstart[accepted[i]]) ;
    hSelectorPlane_->Fill(std::min(distBest,distCurrent), std::max(distBest,distCurrent));
    if(distCurrent < distBest) {
      selected = accepted[i];
    }
    else if(distCurrent == distBest) {
      const double drmuBest = distanceToMuon(evt, selected, muStopUV);
      const double drmuCurrent = distanceToMuon(evt, accepted[i], muStopUV);
      hSelectordrmu_->Fill(std::min(drmuBest,drmuCurrent), std::max(drmuBest, drmuCurrent));
      if(drmuCurrent < drmuBest) {
        selected = accepted[i];
      }
    }
  }

  if(selected >= 0) {
    // Note: resul_->accepted must be set to false higher in the call chain.
    // This class does not see *all* of the incoming events, thus
    // it can't reset the result.
    if(result_) {
      result_->accepted = true;
      result_->momentum = evt.ptot[selected];
    }

    hPerEventMomentum_->Fill(evt.ptot[selected]);
    hSelectedTrackQuality_.fill(evt, selected, distanceToMuon(evt, selected, muStopUV));

    if(doMCTruth_) {
      htruthAccepted_.fill(evt);
    }
  }

  return selected;
}

//================================================================
bool MuCapTrkAnalysisHF::processTrack(int itrack,
                                      const EventClass& evt,
                                      const ROOT::Math::XYPoint& muStopUV,
                                      const TimeWindow& protonWin)
{
  CutNumber c = analyzeTrack(itrack, evt, muStopUV, protonWin);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c==CUTS_ACCEPTED;
}

//================================================================
MuCapTrkAnalysisHF::CutNumber MuCapTrkAnalysisHF::
analyzeTrack(int i, const EventClass& evt,
             const ROOT::Math::XYPoint& muStopUV,
             const TimeWindow& protonWin)
{
  //----------------------------------------------------------------
  if(evt.hefit_ierror[i] != 0) {
    return CUT_IERROR;
  }

  //----------------------------------------------------------------
  trackTime_->Fill(evt.hefit_time[i]);
  trackWinTime_->Fill(evt.hefit_time[i] - protonWin.tstart);
  if(std::abs(evt.hefit_time[i]  - protonWin.tstart) > cutTrackWinTimedt_) {
    return CUT_TRKWINTIME;
  }

  //----------------------------------------------------------------
  hStartStop_->Fill(evt.hefit_pstart[i], evt.hefit_pstop[i]);
  if((cutStream_==TimeWindow::DOWNSTREAM)
     &&((evt.hefit_pstart[i]<=28)||(evt.hefit_pstop[i]<=28))
     ||
     (cutStream_==TimeWindow::UPSTREAM)
     &&((evt.hefit_pstart[i]>28)||(evt.hefit_pstop[i]>28))
     ) {
    return CUT_STREAM;
  }

  // debug:
  if((cutStream_==TimeWindow::DOWNSTREAM) && evt.costh[i] < 0.) {
    std::cout<<"Downstream event with "
             <<" costh = " <<evt.costh[i]
             <<", charge = "<<evt.hefit_q[i]
             <<": run = "<<evt.nrun
             <<" event "<<evt.nevt
             <<std::endl;
  }

  //----------------------------------------------------------------
  hCharge_->Fill(evt.hefit_q[i]);
  if(evt.hefit_q[i] != cutCharge_) {
    return CUT_CHARGE;
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
  hTrackQuality_.fill(evt, i, dr);

  if(dr > cutTrackMuonOffset_) {
    return CUT_TRACK_MUON_OFFSET;
  }

  //----------------------------------------------------------------
  trackChi2_->Fill(evt.hefit_chi2[i]);
//  if(evt.hefit_chi2[i] > cutChi2_) {
//    return CUT_CHI2;
//  }

  //----------------------------------------------------------------
  final_helixCenterUV_->Fill(evt.hefit_ucenter[i], evt.hefit_vcenter[i]);

  const double rout =
    sqrt(std::pow(evt.hefit_ucenter[i], 2) + std::pow(evt.hefit_vcenter[i], 2))
    + evt.radius[i];

  final_trackz_->Fill(evt.hefit_z[i]);
  final_trackROut_->Fill(rout);
  final_costhVsPtot_->Fill(evt.ptot[i], evt.costh[i]);
  final_trackRL_->Fill(evt.radius[i], evt.wavelen[i]);
  final_u0v0_->Fill(evt.hefit_u0[i], evt.hefit_v0[i]);
  final_trackTime_->Fill(evt.hefit_time[i]);
  final_trackWinTime_->Fill(evt.hefit_time[i] - protonWin.tstart);
  final_StartStop_->Fill(evt.hefit_pstart[i], evt.hefit_pstop[i]);

  return CUTS_ACCEPTED;
}
