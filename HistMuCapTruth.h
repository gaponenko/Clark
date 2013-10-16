// These are the truth-only plots that can be made without any event
// selection.  That is, the do not use reconstructed quantities.
//
// Andrei Gaponenko, 2013

#ifndef HistMuCapTruth_h
#define HistMuCapTruth_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistMuCapTruth {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const ConfigFile &conf);

  void fill(const EventClass& evt);

  HistMuCapTruth()
    : hNumMCCaptureTracks_()
    , hptot_()
    , hek_()
    , hphi_()
    , hpcos_()
    , hVUend_()
    , hRZend_()
    , hRendVsPstart_()
    , hZendVsPstart_()
  {}

private :
  TH1 *hNumMCCaptureTracks_;
  TH1 *hptot_;
  TH1 *hek_;
  TH1 *hphi_;
  TH2 *hpcos_;
  TH2 *hVUend_;
  TH2 *hRZend_;
  TH2 *hRendVsPstart_;
  TH2 *hZendVsPstart_;
};

#endif/*HistMuCapTruth_h*/
