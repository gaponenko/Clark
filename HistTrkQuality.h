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

  void fill(const EventClass& evt, int itrack, double drmu);

  HistTrkQuality()
    : hchi2_()
    , hchi2overndf_()
    , hNDF_vs_chi2_()
    , hprob_()
    , hkineprob_()
    , hchi2ndfdrmu_()
  {}

private :
  TH1 *hchi2_;
  TH1 *hchi2overndf_;
  TH2 *hNDF_vs_chi2_;
  TH1 *hprob_;
  TProfile2D *hkineprob_;
  TProfile2D *hchi2ndfdrmu_;
};

#endif/*HistTrkQuality_h*/
