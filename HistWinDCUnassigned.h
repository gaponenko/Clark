// Andrei Gaponenko, 2014

#ifndef HistWinDCUnassigned_h
#define HistWinDCUnassigned_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class TimeWindowingResults;

//================================================================
class HistWinDCUnassigned {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const TimeWindowingResults& wres);

  HistWinDCUnassigned() : hWinDCUnassignedAll_(), hWinDCUnassignedAfterTrig_(), hWinDCNumUnassignedAfterTrig_() {}

private :
  TH2 *hWinDCUnassignedAll_;
  TH2 *hWinDCUnassignedAfterTrig_;
  TH1 *hWinDCNumUnassignedAfterTrig_;
};

#endif/*HistWinDCUnassigned_h*/
