// Andrei Gaponenko, 2015

#ifndef HistMuCapContainedChannel_h
#define HistMuCapContainedChannel_h

#include <string>

#include "MuCapContainedVars.h"
#include "WireCluster.h"

#include "TAxis.h"

class TH1;
class TH2;
class TH3;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistMuCapContainedChannel {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_POSTRK, "No pos trk");
    ax->SetBinLabel(1+CUT_DIO, "DIO veto");
    ax->SetBinLabel(1+CUT_CONTAINED, "Containment cut");
    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_POSTRK,
    CUT_DIO,
    CUT_CONTAINED,
    CUTS_ACCEPTED,
    CUTS_END
  };

  void init(HistogramFactory& hf,
            const std::string& hgrandtopdir,
            const std::string& channelsetname,
            const DetectorGeo& geom,
            const ConfigFile& conf,
            MuCapContainedVars::IVarProcessor& cvp);

  bool accepted(const EventClass& evt,
                bool referenceSampleAccepted,
                int iPosTrack,
                int iNegTrack,
                const ClustersByPlane& globalPlaneClusters );

  HistMuCapContainedChannel() : doMCTruth_(false), cvp_(0) {}

private :
  bool doMCTruth_;

  MuCapContainedVars::IVarProcessor* cvp_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH2* reco_;
  TH2* reco_mcproton_;
  TH2* reco_mcdeuteron_;
  TH2* reco_mctriton_;
  TH2* reco_mcalpha_;
  TH2* reco_mcdio_;

  TH3* migration_;
  TH3* migration_mcproton_;
  TH3* migration_mcdeuteron_;
  TH3* migration_mctriton_;
  TH3* migration_mcalpha_;

  // Contamination: like migration, but for events not in the reference sample: no true tgt stop.
  TH3* contamination_;
  TH3* contamination_mcproton_;
  TH3* contamination_mcdeuteron_;
  TH3* contamination_mctriton_;
  TH3* contamination_mcalpha_;

  CutNumber analyzeEvent(const EventClass& evt,
                         bool referenceSampleAccepted,
                         int iPosTrack,
                         int iNegTrack,
                         const ClustersByPlane& globalPlaneClusters );
};

#endif/*HistMuCapContainedChannel_h*/
