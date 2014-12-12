// Andrei Gaponenko, 2014

#ifndef HistMuCapAnalysisChannels_h
#define HistMuCapAnalysisChannels_h

#include <string>

#include "HistHitBasedAnalysis.h"
#include "HistMuCapTruth.h"
#include "HistMuCapTrkResolution.h"
#include "WireCluster.h"
#include "HistTDCBCSWidth.h"

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
            double rangePIDVar,
            const ClustersByPlane& globalPlaneClusters );

  HistMuCapAnalysisChannels() : doMCTruth_(false) {}

private :
  bool doMCTruth_;

  //----------------------------------------------------------------
  // Essential "channel" analysis histograms
  TH2* contained_prange_;
  TH1* uncontained_p_;

  // lost with 2 channels
  TH1* mclost2_ptot_;
  TH1* mclost2_count_; // for DIO, where we don't have mcptot easily available

  TH3* containedMigration_;
  TH2* uncontainedMigration_;

  HistHitBasedAnalysis hitbased_;

  // lost with 3 channels
  TH1* mclost3_ptot_;
  TH1* mclost3_count_;

  //----------------------------------------------------------------
  // Extra histograms to plot efficiencies etc.  Not essential for the
  // unfolding.

  TH1* mcin_proton_ptot_;
  TH1* mcin_deuteron_ptot_;
  TH1* mcin_dio_count_;

  // same content as the reco hists, but which one is filled depends on MC PID
  TH2* contained_prange_mcproton_;
  TH2* contained_prange_mcdeuteron_;
  TH2* contained_prange_mcdio_;
  TH1* uncontained_p_mcproton_;
  TH1* uncontained_p_mcdeuteron_;
  TH1* uncontained_p_mcdio_;

  HistMuCapTruth hTruthTrkContained_;
  HistMuCapTruth hTruthTrkUncontained_;

  HistMuCapTrkResolution hResolutionContained_;
  HistMuCapTrkResolution hResolutionUncontained_;

  HistTDCBCSWidth hTDCWidthContained_;
  HistTDCBCSWidth hTDCWidthUncontained_;
};

#endif/*HistMuCapAnalysisChannels_h*/
