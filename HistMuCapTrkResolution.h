// Andrei Gaponenko, 2014

#ifndef HistMuCapTrkResolution_h
#define HistMuCapTrkResolution_h

#include <string>

class TH1;
class TH2;
class TProfile;
class TProfile2D;

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
    , hMomBias2DMC_()
    , hMomBias2DReco_()
  {}

private :
  TH1 *hGlobalResolution_;
  TProfile *hMomResVsMom_;
  TProfile *hMomResVsCosth_;
  TProfile2D *hMomBias2DMC_;
  TProfile2D *hMomBias2DReco_;
};

#endif/*HistMuCapTrkResolution_h*/
