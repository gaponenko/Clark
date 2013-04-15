// Extract (uv) distribution of selected muons using reconstructed DIO tracks.
//
// Andrei Gaponenko, 2013

#ifndef MuCapUVAnalysis_h
#define MuCapUVAnalysis_h

#include <string>
#include <fstream>

#include "WireCluster.h"

#include "TAxis.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class MuCapUVAnalysis {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_IERROR, "ierror");
    ax->SetBinLabel(1+CUT_CHARGE, "charge");
    ax->SetBinLabel(1+CUT_TIME, "time");
    ax->SetBinLabel(1+CUT_RADIUS, "radius");

    ax->SetBinLabel(1+CUT_COSTHETAMIN, "cos(theta) min");
    ax->SetBinLabel(1+CUT_COSTHETAMAX, "cos(theta) max");
    ax->SetBinLabel(1+CUT_PTMIN, "pt min");
    ax->SetBinLabel(1+CUT_PZMIN, "pz min");
    ax->SetBinLabel(1+CUT_PTOTMIN, "ptot min");
    ax->SetBinLabel(1+CUT_PTOTMAX, "ptot max");

    ax->SetBinLabel(1+CUT_TRACK_MUON_OFFSET, "drmu");
    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_IERROR,
    CUT_CHARGE,
    CUT_TIME,
    CUT_RADIUS,
    CUT_COSTHETAMIN,
    CUT_COSTHETAMAX,
    CUT_PTMIN,
    CUT_PZMIN,
    CUT_PTOTMIN,
    CUT_PTOTMAX,
    CUT_TRACK_MUON_OFFSET,
    CUTS_ACCEPTED,
    CUTS_END
  };

  void init(const std::string& hdir, HistogramFactory &hf, const ConfigFile &conf);

  void process(const EventClass& evt,
               double timeWinStart,
               const ClustersByPlane& globalClusters,
               double muStopU, double muStopV
               );

  MuCapUVAnalysis()
    : cutCharge_()
    , cutTrackWinTimeDiff_()
    , cutTrackRmax_()
    , cutCosThetaMin_()
    , cutCosThetaMax_()
    , cutPtMin_()
    , cutPzMin_()
    , cutPtotMin_()
    , cutPtotMax_()
    , cutTrackMuonOffset_()

    , h_cuts_r()
    , h_cuts_p()

    , trackwintimeLargeScale_()
    , trackwintime_()

    , trackRL_()
    , costhVsPtot_()

    , hStartStop_()
    , hHitRange_()
    , trackz_()

    , trackMuonOffset_()
    , trackMuondr_()

    , helixCenterUV_()

    , final_trackRL_()
    , final_costhVsPtot_()
    , final_u0v0_()
    , final_trackROut_()

    , hNumTracks_()
  {}

private :
  int cutCharge_;
  double cutTrackWinTimeDiff_;
  double cutTrackRmax_;
  double cutCosThetaMin_;
  double cutCosThetaMax_;
  double cutPtMin_;
  double cutPzMin_;
  double cutPtotMin_;
  double cutPtotMax_;
  double cutTrackMuonOffset_;

  std::string uvOutFileName_;
  std::ofstream uvOutFile_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1 *trackwintimeLargeScale_;
  TH1 *trackwintime_;

  TH2 *trackRL_;
  TH2 *costhVsPtot_;

  TH2 *hStartStop_;
  TH2 *hHitRange_;
  TH1 *trackz_;
  TH2 *trackMuonOffset_;
  TH1 *trackMuondr_;

  TH2 *helixCenterUV_;

  TH2 *final_trackRL_; // radius and wavelength
  TH2 *final_costhVsPtot_;
  TH2 *final_u0v0_;
  TH1 *final_trackROut_;

  TH1 *hNumTracks_;

  // true iff the track passed the cuts
  bool processTrack(int itrack,
                    const EventClass& evt,
                    double timeWinStart,
                    const ClustersByPlane& globalClusters,
                    double muStopU, double muStopV
                    );

  CutNumber analyzeTrack(int itrack,
                         const EventClass& evt,
                         double timeWinStart,
                         const ClustersByPlane& globalClusters,
                         double muStopU, double muStopV
                         );
};

#endif/*MuCapUVAnalysis_h*/
