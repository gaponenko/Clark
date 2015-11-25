// Andrei Gaponenko, 2014

#ifndef HistMuCapMCTgtStops_h
#define HistMuCapMCTgtStops_h

#include <string>

class TH1;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistMuCapMCTgtStops {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const EventClass& evt);

private :
  const DetectorGeo *geom_;

  TH1* hMcMuonTotalMultiplicity_;
  TH1* hMcMuonTrigCandidateMultiplicity_;

  TH1* hMcMuonStopTime_;

  TH1* hMcMuonZStopFine_;
  TH1* hMcMuonZStopCoarse_;
};

#endif/*HistMuCapMCTgtStops_h*/
