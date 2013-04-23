// These are the plots that combine MC and reconstructed quantities.
//
// Andrei Gaponenko, 2013

#ifndef HistMuCapRTruth_h
#define HistMuCapRTruth_h

#include <string>

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistMuCapRTruth {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const ConfigFile &conf);

  void fill(const EventClass& evt, int lastPlane, double extrapolatedRmax);

  HistMuCapRTruth()
    : hLastPlaneVsMCPstart_()
    , hRmaxContained_()
    , hRmaxUncontained_()
  {}

private :
  TH2 *hLastPlaneVsMCPstart_;
  TH1 *hRmaxContained_;
  TH1 *hRmaxUncontained_;
};

#endif/*HistMuCapRTruth_h*/
