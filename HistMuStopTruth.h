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
class ConfigFile;
class EventClass;

//================================================================
class HistMuStopTruth {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const ConfigFile &conf);

  void fill(const EventClass& evt);

  HistMuStopTruth()
    : hstopZ1_()
    , hstopZ2_()
    , hstopZ3_()
  {}

private :
  TH1 *hnumPrimaries_;
  TH1 *hstopZ1_;
  TH1 *hstopZ2_;
  TH1 *hstopZ3_;
};

#endif/*HistMuStopTruth_h*/
