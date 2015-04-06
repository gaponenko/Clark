// Check that a track does not exit the detector but ranges out inside.
// The code assumes downstream tracks.
//
// Andrei Gaponenko, 2014

#ifndef MuCapTrkContainmentCut_h
#define MuCapTrkContainmentCut_h

#include <string>

#include "TAxis.h"

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class EventClass;

//================================================================
class MuCapTrkContainmentCut {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_NO_TRACK, "no track");
    ax->SetBinLabel(1+CUT_ZRANGE, "z range");
    ax->SetBinLabel(1+CUT_ROUT, "r out");
    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_NO_TRACK,
    CUT_ZRANGE,
    CUT_ROUT,
    CUTS_ACCEPTED,
    CUTS_END
  };

  void init(const std::string& hdir,
            HistogramFactory &hf,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  bool contained(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters, int *outLastPlane=0);

  MuCapTrkContainmentCut()
    : cutMaxPlane_()
    , cutMaxRout_()

    , h_cuts_r()
    , h_cuts_p()

    , hExtendedLastPlane_()
    , hExtendedLastPlaneFinal_()
    , hRout_()
    , hRoutFinal_()
    , rangeFinder_()
  {}

private :
  bool doMCTruth_;
  int cutMaxPlane_;
  double cutMaxRout_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1 *hExtendedLastPlane_;
  TH1 *hExtendedLastPlaneFinal_;
  TH1 *hRout_;
  TH1 *hRoutFinal_;

  CutNumber analyzeTrack(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters, int *outLastPlane);
  int (*rangeFinder_)(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters);
};

#endif/*MuCapTrkContainmentCut_h*/
