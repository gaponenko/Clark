// Andrei Gaponenko, 2014

#ifndef HistMuCapTrkResolution_h
#define HistMuCapTrkResolution_h

#include <string>

class TH1;
class TH2;
class TProfile;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistMuCapTrkResolution {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const ConfigFile& conf);

  void fill(const EventClass& evt, int iRecoTrack);

  HistMuCapTrkResolution()
    : hGlobalResolution_()
    , hMomResVsMom_()
    , hMomResVsCosth_()
  {}

private :
  TH1 *hGlobalResolution_;
  TProfile *hMomResVsMom_;
  TProfile *hMomResVsCosth_;
};

#endif/*HistMuCapTrkResolution_h*/
