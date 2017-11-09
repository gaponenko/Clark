// These are the truth-only plots that can be made without any event
// selection.  That is, the do not use reconstructed quantities.
//
// Andrei Gaponenko, 2013

#ifndef HistMuStopTruth_h
#define HistMuStopTruth_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class EventClass;

//================================================================
class HistMuStopTruth {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const EventClass& evt);

  HistMuStopTruth()
    : zTargetCenter_()
    , hMcMuonTotalMultiplicity_()
    , hMcMuonTrigCandidateMultiplicity_()
    , hMcMuonStopTime_()
    , hstopZ1_()
    , hstopZ2_()
    , hstopdz_()
  {}

private :
  double zTargetCenter_;

  TH1* hMcMuonTotalMultiplicity_;
  TH1* hMcMuonTrigCandidateMultiplicity_;
  TH1* hMcMuonStopTime_;

  TH1 *hstopZ1_;
  TH1 *hstopZ2_;
  TH1 *hstopdz_;
};

#endif/*HistMuStopTruth_h*/
