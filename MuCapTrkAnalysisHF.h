// A track-based capture analysis that used helix fit results as inputs
//
// Andrei Gaponenko, 2013

#ifndef MuCapTrkAnalysisHF_h
#define MuCapTrkAnalysisHF_h

#include <string>

#include "TimeWindow.h"
#include "WireCluster.h"
#include "HistMuCapTruth.h"
#include "RecoResMuCapTrk.h"
#include "HistTrkQuality.h"
#include "HistMuCapTrkResolution.h"

#include "TAxis.h"
#include "Math/Point2D.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class MuCapTrkAnalysisHF {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_IERROR, "ierror");
    ax->SetBinLabel(1+CUT_CHARGE, "charge");
    ax->SetBinLabel(1+CUT_TRKWINTIME, "trkwintime");
    ax->SetBinLabel(1+CUT_STREAM, "stream");
    ax->SetBinLabel(1+CUT_RADIUS, "radius");

    ax->SetBinLabel(1+CUT_COSTHETAMIN, "cos(theta) min");
    ax->SetBinLabel(1+CUT_COSTHETAMAX, "cos(theta) max");
    ax->SetBinLabel(1+CUT_PTMIN, "pt min");
    ax->SetBinLabel(1+CUT_PZMIN, "pz min");
    ax->SetBinLabel(1+CUT_PTOTMIN, "ptot min");
    ax->SetBinLabel(1+CUT_PTOTMAX, "ptot max");
    ax->SetBinLabel(1+CUT_TRACK_MUON_OFFSET, "drmu");
    // ax->SetBinLabel(1+CUT_CHI2,     "chi2");

    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_IERROR,
    CUT_CHARGE,
    CUT_TRKWINTIME,
    CUT_STREAM,
    CUT_RADIUS,
    CUT_COSTHETAMIN,
    CUT_COSTHETAMAX,
    CUT_PTMIN,
    CUT_PZMIN,
    CUT_PTOTMIN,
    CUT_PTOTMAX,
    CUT_TRACK_MUON_OFFSET,
    // CUT_CHI2,
    CUTS_ACCEPTED,
    CUTS_END
  };

  void init(const std::string& hdir,
            HistogramFactory &hf,
            const DetectorGeo& geom,
            const ConfigFile &conf,
            const std::string& cutSet,
            TimeWindow::StreamType cutStream,
            RecoResMuCapTrk *result=0);

  // Returns the index of an accepted track in the event, or -1
  int process(const EventClass& evt,
              const ROOT::Math::XYPoint& muStopUV,
              const TimeWindow& protonWin);

  MuCapTrkAnalysisHF()
    : doMCTruth_(false)
    , result_(nullptr)
    , cutCharge_()
    , cutTrackWinTimedt_()
    , cutTrackRmax_()
    , cutCosThetaMin_()
    , cutCosThetaMax_()
    , cutPtMin_()
    , cutPzMin_()
    , cutPtotMin_()
    , cutPtotMax_()
      // , cutChi2_()
    , cutTrackMuonOffset_()

    , h_cuts_r()
    , h_cuts_p()

    , trackTime_()
    , trackWinTime_()
    , trackRL_()
    , costhVsPtot_()
    , pfine_()

    , hStartStop_()

    , trackMuonOffset_()
    , trackMuondr_()

    , final_helixCenterUV_()
    , final_trackz_()
    , final_trackRL_()
    , final_costhVsPtot_()
    , final_pfine_()
    , final_u0v0_()
    , final_trackROut_()
    , final_trackTime_()
    , final_trackWinTime_()
    , final_StartStop_()

    , hNumTracks_()
    , hPerEventMomentum_()
    , hPerEventTime_()

    , hSelectorPlane_()
    , hSelectordrmu_()
  {}

private :
  bool doMCTruth_;
  RecoResMuCapTrk *result_;
  int cutCharge_;
  TimeWindow::StreamType cutStream_;
  double cutTrackWinTimedt_;
  double cutTrackRmax_;
  double cutCosThetaMin_;
  double cutCosThetaMax_;
  double cutPtMin_;
  double cutPzMin_;
  double cutPtotMin_;
  double cutPtotMax_;
  // double cutChi2_;
  double cutTrackMuonOffset_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1 *hCharge_;
  TH1 *trackTime_;
  TH1 *trackWinTime_;

  TH2 *trackRL_;
  TH2 *costhVsPtot_;
  TH1 *pfine_;

  TH2 *hStartStop_;
  TH1 *trackChi2_;
  HistTrkQuality hTrackQuality_;
  TH2 *trackMuonOffset_;
  TH1 *trackMuondr_;

  TH2 *final_helixCenterUV_;
  TH1 *final_trackz_;
  TH2 *final_trackRL_; // radius and wavelength
  TH2 *final_costhVsPtot_;
  TH1 *final_pfine_;
  TH2 *final_u0v0_;
  TH1 *final_trackROut_;
  TH1 *final_trackTime_;
  TH1 *final_trackWinTime_;
  TH2 *final_StartStop_;

  // Unlike the histos above, which are per track,
  // the following distributions are filled once per event,
  TH1 *hNumTracks_;
  TH1 *hPerEventMomentum_;
  TH1 *hPerEventTime_;
  HistMuCapTruth htruthAccepted_;
  HistTrkQuality hSelectedTrackQuality_;
  HistMuCapTrkResolution hSelectedTrkResolution_;

  // Selection among multiple track passing the cuts
  TH2 *hSelectorPlane_;
  TH2 *hSelectordrmu_;

  // true iff the track passed the cuts
  bool processTrack(int itrack,
                    const EventClass& evt,
                    const ROOT::Math::XYPoint& muStopUV,
                    const TimeWindow& protonWin);

  CutNumber analyzeTrack(int itrack,
                         const EventClass& evt,
                         const ROOT::Math::XYPoint& muStopUV,
                         const TimeWindow& protonWin);
};

#endif/*MuCapTrkAnalysisHF_h*/
