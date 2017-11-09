// Andrei Gaponenko, 2014

#ifndef HistMuCapAnalysisChannels_h
#define HistMuCapAnalysisChannels_h

#include <string>

#include "WireCluster.h"

#include "HistMuCapContainedChannel.h"
#include "HistMuCapUncontainedChannel.h"
#include "HistMuCapHitbasedChannel.h"

#include "HistTDCBCSWidth.h"
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
            const std::string& htopdir,
            const std::string& channelsetname,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  // Must be called for all events.  Must be called before fill() on the same event.
  void fillReferenceSample(const EventClass& evt);

  bool referenceSampleAccepted() const { return referenceSampleAccepted_; }

  void fill(const EventClass& evt,
            int iPosTrack,
            int iNegTrack,
            const ClustersByPlane& globalPlaneClusters);

  HistMuCapAnalysisChannels()
    : doMCTruth_(false)
    , targetCenterZ_(0.)
    , targetThickness_(0.)
    , referenceSample_nrun_(0)
    , referenceSample_nevt_(0)
  {}

private :
  bool doMCTruth_;

  double targetCenterZ_;
  double targetThickness_;

  int referenceSample_nrun_;
  int referenceSample_nevt_;
  bool referenceSampleAccepted_;

  TH1* refsample_in_zstop_;
  TH1* refsample_accepted_count_;

  //----------------------------------------------------------------
  // Essential "channel" analysis histograms
  HistMuCapContainedChannel contained_;
  HistMuCapUncontainedChannel uncontained_;
  HistMuCapHitbasedChannel hitbased_;

  // Reference sample spectrum: common for all reco channels
  // Essential for normalization of the unfolding
  TH1* mcin_proton_ptot_;
  TH1* mcin_deuteron_ptot_;
  TH1* mcin_triton_ptot_;
  TH1* mcin_alpha_ptot_;
  TH1* mcin_dio_count_;

  //----------------
  // Extra distributions
  bool fillExtras_TDC_;
  HistTDCBCSWidth hTDCWidthContained_;
  HistTDCBCSWidth hTDCWidthUncontained_;
  HistTDCBCSWidth hTDCWidthHitbased_;
  HistTDCBCSWidth hTDCWidthNone_;

  HistMuCapTruth hTruthContained_;
  HistMuCapTruth hTruthUncontained_;
  HistMuCapTruth hTruthHitbased_;
  HistMuCapTruth hTruthNone_;

  HistMuCapTrkResolution hResolutionContained_;
  HistMuCapTrkResolution hResolutionUncontained_;
};

#endif/*HistMuCapAnalysisChannels_h*/
