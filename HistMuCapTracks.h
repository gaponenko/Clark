// Andrei Gaponenko, 2013

#ifndef HistMuCapTracks_h
#define HistMuCapTracks_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistMuCapTracks {
public:
  void init(const std::string& hdir,
            const std::string& namePrefix,
            int chargeSign,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const EventClass& evt, double timeWinStart, double muu, double muv);

  HistMuCapTracks()
    : cutCharge_()
    , cutTrackWinTimeDiff_()
    , cutTrackRmax_()
    , cutTrackMuonOffset_()
    , trackwintimeLargeScale_()
    , trackwintime_()
    , hStartStop_()
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

  TH1 *trackwintimeLargeScale_;
  TH1 *trackwintime_;
  TH2 *hStartStop_;
  TH1 *trackz_;
  TH2 *trackMuonOffset_;
  TH1 *trackMuondr_;

  TH2 *costhVsPtot_;
  TH2 *u0v0_;
  TH2 *trackRL_; // radius and wavelength
  TH2 *helixCenterUV_;
  TH1 *trackROut_;

  TH1 *hNumTracks_;
};

#endif/*HistMuCapTracks_h*/
