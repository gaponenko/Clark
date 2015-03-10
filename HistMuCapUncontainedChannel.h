// Andrei Gaponenko, 2015

#ifndef HistMuCapUncontainedChannel_h
#define HistMuCapUncontainedChannel_h

#include <string>

#include "TAxis.h"

class TH1;
class TH2;
class TH3;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistMuCapUncontainedChannel {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_POSTRK, "No pos trk");
    ax->SetBinLabel(1+CUT_DIO, "DIO veto");
    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_POSTRK,
    CUT_DIO,
    CUTS_ACCEPTED,
    CUTS_END
  };

  void init(HistogramFactory& hf,
            const std::string& hgrandtopdir,
            const std::string& channelsetname,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  bool accepted(const EventClass& evt,
                bool referenceSampleAccepted,
                int iPosTrack,
                int iNegTrack);

  HistMuCapUncontainedChannel() : doMCTruth_(false) {}

private :
  bool doMCTruth_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1* reco_;
  TH1* reco_mcproton_;
  TH1* reco_mcdeuteron_;
  TH1* reco_mcdio_;

  TH2* migration_;
  TH2* migration_mcproton_;
  TH2* migration_mcdeuteron_;

  // Contamination: like migration, but for events not in the reference sample: no true tgt stop.
  TH2* contamination_;
  TH2* contamination_mcproton_;
  TH2* contamination_mcdeuteron_;

  CutNumber analyzeEvent(const EventClass& evt,
                         bool referenceSampleAccepted,
                         int iPosTrack,
                         int iNegTrack );
};

#endif/*HistMuCapUncontainedChannel_h*/
