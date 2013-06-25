// Andrei Gaponenko, 2013

#ifndef HistWinTime_h
#define HistWinTime_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class TimeWindowingResults;

//================================================================
class HistWinTime {
public:
  void init(const std::string& hdir,
            const std::string& nameSuffix,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const TimeWindowingResults& wres);

  HistWinTime() : hWinTime_()  {}

private :
  TH1 *hWinTime_;
};

#endif/*HistWinTime_h*/
