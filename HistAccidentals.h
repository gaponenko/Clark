// Andrei Gaponenko, 2013

#ifndef HistAccidentals_h
#define HistAccidentals_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class TimeWindowingResults;

//================================================================
class HistAccidentals {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const TimeWindowingResults& wres);

  HistAccidentals()
    : cutPreTrigTimeMin_(), cutPreTrigTimeMax_()
    , htstartAll_(), htstartDn_()
    , hnumwinAll_(), hnumwinDn_()
  {}

private :
  // Only consider time windows starting in this range
  double cutPreTrigTimeMin_;
  double cutPreTrigTimeMax_;

  TH1 *htstartAll_;
  TH1 *htstartDn_;
  TH1 *hnumwinAll_;
  TH1 *hnumwinDn_;
};

#endif/*HistAccidentals_h*/
