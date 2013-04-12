// Extract (uv) distribution of selected muons using reconstructed DIO tracks.
//
// Andrei Gaponenko, 2013

#ifndef MuCapUVAnalysis_h
#define MuCapUVAnalysis_h

#include <string>

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
    ax->SetBinLabel(1+CUT_TRACK_MUON_OFFSET, "drmu");
    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_IERROR,
    CUT_CHARGE,
    CUT_TIME,
    CUT_RADIUS,
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
    , cutTrackMuonOffset_()

    , h_cuts_r()
    , h_cuts_p()

    , trackwintimeLargeScale_()
    , trackwintime_()
    , hStartStop_()
    , hHitRange_()
    , trackz_()
    , trackMuonOffset_()
    , trackMuondr_()
    , costhVsPtot_()
    , u0v0_()
    , trackRL_()
    , helixCenterUV_()
    , trackROut_()
    , hNumTracks_()
  {}

private :
  int cutCharge_;
  double cutTrackWinTimeDiff_;
  double cutTrackRmax_;
  double cutTrackMuonOffset_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1 *trackwintimeLargeScale_;
  TH1 *trackwintime_;
  TH2 *hStartStop_;
  TH2 *hHitRange_;
  TH1 *trackz_;
  TH2 *trackMuonOffset_;
  TH1 *trackMuondr_;

  TH2 *costhVsPtot_;
  TH2 *u0v0_;
  TH2 *trackRL_; // radius and wavelength
  TH2 *helixCenterUV_;
  TH1 *trackROut_;

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
