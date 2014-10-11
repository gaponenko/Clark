// Andrei Gaponenko, 2014

#ifndef HistMuCapRooUnfold_h
#define HistMuCapRooUnfold_h

#include <string>

//#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfoldDummy.h"
typedef RooUnfoldDummy RooUnfoldResponse;
class TH1D;
class TH2D;
class TH3D;


#define NBUNFOLD1D 11

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;



//================================================================
class HistMuCapRooUnfold {
public:
  class HistUnfold1D {
    private:
      bool Selected;
      RooUnfoldResponse Response_;
      TH1D* hTruthMomentum_;
      TH1D* hTruthMomentumReco_;
      TH2D* hMeasVsTruthMomentum_;
      TH1D* hTruthMomentumNotReco_;
      TH1D* hMeasuredMomentum_;

    public:
      HistUnfold1D(HistogramFactory &H, std::string Dir, std::string Name, bool MCTruth, int NBinP, double MaxP);
      void Reset();
      void FillMeasured(double mom);
      void FillTruth(double reco, double tru);
      void MissTruth(double tru);

  };

  class HistUnfold2D {
    private:
      bool Selected;
      RooUnfoldResponse Response_;
      TH2D* hTruthMomentum_;
      TH2D* hTruthMomentumReco_;
      TH2D* hMeasVsTruthMomentumTruProtons_;
      TH2D* hMeasVsTruthMomentumTruDeuterons_;
      TH2D* hTruthMomentumNotReco_;
      TH2D* hMeasuredMomentum_;
      TH2D* hPlaneVsMomentum_;
      TH2D* hPlaneVsMomentumProt_;
      TH2D* hPlaneVsMomentumDeut_;
      TH2D* hPlaneCosThVsMomentum_;
      TH2D* hPlaneCosThVsMomentumProt_;
      TH2D* hPlaneCosThVsMomentumDeut_;

    public:
      HistUnfold2D(HistogramFactory &H, std::string Dir, std::string Name, bool MCTruth, int NBinP, double MaxP);
      void Reset();
      void FillMeasured(int PID, double mom, int lastPlane, double lastPlaneOvCosTh);
      void FillTruth(int recoPID, int truPID, double reco, double tru);
      void MissTruth(int truPID, double tru);

  };

  void init(HistogramFactory& H,
            const std::string& hdir,
            const ConfigFile& conf);

  // Must be called for all events.  Must be called BEFORE analyze() in process
  void SaveEventVariables(const EventClass& evt);

  void FillMeasured(const EventClass& evt,
            int iPosTrack,
            int iNegTrack,
            bool isPosTrackContained,
            double rangePIDVar);

  void Fill(const EventClass& evt,
            int iPosTrack,
            int iNegTrack,
            bool isPosTrackContained,
            double rangePIDVar);

  // Must be called for all events.  Must be called AFTER analyze() in process
  void FillAndMiss();

  HistMuCapRooUnfold() : doMCTruth_(false) {}

private :
  bool doMCTruth_;
  double p_true_;
  int TruePIDProton_;
  bool truCaptEvt_;

  HistUnfold1D *FullSpectrum_;
  HistUnfold1D *Contained_;
  HistUnfold1D *PlnRngCutPln_;
  HistUnfold1D *PlnVsPCutZone1_;
  HistUnfold1D *PlnVsPCutZone2_;
  HistUnfold1D *PlnVsPCutZone3_;
  HistUnfold1D *PlnVsPCutZone4_;
  HistUnfold1D *PlnVsPCutZone5_;
  HistUnfold1D *PlnVsPCutZone6_;
  HistUnfold1D *PlnVsPCutZone1AllTrks_;
  HistUnfold1D *PlnVsPCutZone2AllTrks_;

  HistUnfold1D *AllUnfold1D_[NBUNFOLD1D];

  HistUnfold2D *WithPID_;
  HistUnfold2D *WithPIDRestrictPlane_;
  HistUnfold2D *WithPIDAllTrks_;

};

#endif/*HistMuCapRooUnfold_h*/
