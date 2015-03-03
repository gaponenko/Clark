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

  // Must be called for all events.  Must be called before fill() on the same event.
  void fillReferenceSample(const EventClass& evt);

  bool referenceSampleAccepted() const { return referenceSampleAccepted_; }

  void fill(const EventClass& evt,
            int iPosTrack,
            int iNegTrack,
            bool isPosTrackContained,
            double rangePIDVar,
            const ClustersByPlane& globalPlaneClusters );

  HistMuCapAnalysisChannels() : doMCTruth_(false), referenceSample_nrun_(0), referenceSample_nevt_(0) {}

private :
  bool doMCTruth_;

  int referenceSample_nrun_;
  int referenceSample_nevt_;
  bool referenceSampleAccepted_;

  TH1* refsample_muminus_multiplicity_;
  TH1* refsample_endvtx_time_;
  TH1* refsample_num_stops_;

  TH1* refsample_in_zstop_;
  TH1* refsample_accepted_count_;

  // Reference sample spectrum: common for all reco channels
  TH1* mcin_proton_ptot_;
  TH1* mcin_deuteron_ptot_;
  TH1* mcin_dio_count_;


  //----------------------------------------------------------------
  // Essential "channel" analysis histograms
  TH2* contained_prange_;
  TH1* uncontained_p_;

  TH3* containedMigration_;
  TH3* containedMigration_mcproton_;
  TH3* containedMigration_mcdeuteron_;

  TH2* uncontainedMigration_;
  TH2* uncontainedMigration_mcproton_;
  TH2* uncontainedMigration_mcdeuteron_;

  // Contamination: like migration, but for events not in the reference sample: no true tgt stop.
  TH3* containedContamination_;
  TH3* containedContamination_mcproton_;
  TH3* containedContamination_mcdeuteron_;

  TH2* uncontainedContamination_;
  TH2* uncontainedContamination_mcproton_;
  TH2* uncontainedContamination_mcdeuteron_;

  HistHitBasedAnalysis hitbased_;

  // same content as the reco hists, but which one is filled depends on MC PID
  TH2* contained_prange_mcproton_;
  TH2* contained_prange_mcdeuteron_;
  TH2* contained_prange_mcdio_;
  TH1* uncontained_p_mcproton_;
  TH1* uncontained_p_mcdeuteron_;
  TH1* uncontained_p_mcdio_;

  //----------------

  HistMuCapTruth hTruthTrkContained_;
  HistMuCapTruth hTruthTrkUncontained_;

  HistMuCapTrkResolution hResolutionContained_;
  HistMuCapTrkResolution hResolutionUncontained_;

  HistTDCBCSWidth hTDCWidthContained_;
  HistTDCBCSWidth hTDCWidthUncontained_;
};

#endif/*HistMuCapAnalysisChannels_h*/
