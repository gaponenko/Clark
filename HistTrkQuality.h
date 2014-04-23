// Andrei Gaponenko, 2014

#ifndef HistTrkQuality_h
#define HistTrkQuality_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistTrkQuality {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const ConfigFile& conf);

  void fill(const EventClass& evt, int itrack);

  HistTrkQuality()
    : hchi2_()
    , hNDF_vs_chi2_()
    , hprob_()
    , hkineprob_()
  {}

private :
  TH1 *hchi2_;
  TH2 *hNDF_vs_chi2_;
  TH1 *hprob_;
  TProfile2D *hkineprob_;
};

#endif/*HistTrkQuality_h*/
