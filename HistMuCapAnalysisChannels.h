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

    , mclost2_ptot_()
    , mclost2_count_()

    , mcin_proton_ptot_()
    , mcin_deuteron_ptot_()
    , mcin_dio_count_()

    , containedMigration_()
    , uncontainedMigration_()

  {}

private :
  bool doMCTruth_;

  // "channel" analysis histograms
  TH2* contained_prange_;
  TH1* uncontained_p_;

  // lost with 2 channels
  TH1* mclost2_ptot_;
  TH1* mclost2_count_; // for DIO, where we don't have mcptot easily available

  TH1* mcin_proton_ptot_;
  TH1* mcin_deuteron_ptot_;
  TH1* mcin_dio_count_;

  TH3* containedMigration_;
  TH2* uncontainedMigration_;

  HistMuCapTruth hTruthTrkContained_;
  HistMuCapTruth hTruthTrkUncontained_;

  HistMuCapTrkResolution hResolutionContained_;
  HistMuCapTrkResolution hResolutionUncontained_;
};

#endif/*HistMuCapAnalysisChannels_h*/
