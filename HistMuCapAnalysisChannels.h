// Andrei Gaponenko, 2014

#ifndef HistMuCapAnalysisChannels_h
#define HistMuCapAnalysisChannels_h

#include <string>

#include "HistMuCapTruth.h"
#include "HistMuCapTrkResolution.h"

class TH1;
class TH2;
class TH3;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistMuCapAnalysisChannels {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const EventClass& evt,
            int iPosTrack,
            int iNegTrack,
            bool isPosTrackContained,
            double rangePIDVar);

  HistMuCapAnalysisChannels()
    : doMCTruth_(false)
    , contained_prange_()
    , uncontained_p_()
    , lost1_ptot_()
    , noncapture_lost_()
    , mcproton_ptot_()
    , mcdeuteron_ptot_()
    , noncapture_count_() // "channels" input reco events that are not true captures (e.g. DIO)
    , containedMigration1_()
    , uncontainedMigration1_()

  {}

private :
  bool doMCTruth_;

  // "channel" analysis histograms
  TH2* contained_prange_;
  TH1* uncontained_p_;
  TH1* lost1_ptot_;
  TH1* noncapture_lost_;
  TH1* mcproton_ptot_;
  TH1* mcdeuteron_ptot_;
  TH1* noncapture_count_; // "channels" input reco events that are not true captures (e.g. DIO)
  TH3* containedMigration1_;
  TH2* uncontainedMigration1_;

  HistMuCapTruth hTruthTrkContained_;
  HistMuCapTruth hTruthTrkUncontained_;

  HistMuCapTrkResolution hResolutionContained_;
  HistMuCapTrkResolution hResolutionUncontained_;
};

#endif/*HistMuCapAnalysisChannels_h*/
