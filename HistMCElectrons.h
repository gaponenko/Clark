// Details about truth electrons, to look at G3 delta ray production.
//
// Andrei Gaponenko, 2015

#ifndef HistMCElectrons_h
#define HistMCElectrons_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistMCElectrons {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const ConfigFile &conf);

  void fill(const EventClass& evt);

  HistMCElectrons()
    : hNumAllElectrons_()
    , hNumElectrons_()
    , hEkAll_()
    , hEk_()
    , hZStartAll_()
    , hZStart_()
  {}

private :
  TH1 *hNumAllElectrons_;
  TH1 *hNumElectrons_;
  TH1 *hEkAll_;
  TH1 *hEk_;
  TH1 *hZStartAll_;
  TH1 *hZStart_;
};

#endif/*HistMCElectrons_h*/
