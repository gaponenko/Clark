// Andrei Gaponenko, 2014

#ifndef HistMuCapRooUnfold_h
#define HistMuCapRooUnfold_h

#include <string>

#include "RooUnfold/RooUnfoldResponse.h"

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
	
  void init(HistogramFactory& H,
            const std::string& hdir,
            const ConfigFile& conf);

  // Must be called for all events.  Must be called BEFORE analyze() in process
  void SaveEventVariables(const EventClass& evt);

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

  bool anDnLateResponseWithPIDAllTrks_Filled_;
  bool anDnLateResponseWithPID_Filled_;

  RooUnfoldResponse anDnLateResponseWithPID_;
  TH2D* hWithPIDTruthMomentum_;
  TH2D* hWithPIDTruthMomentumReco_;
  TH2D* hWithPIDMeasVsTruthMomentumTruProtons_;
  TH2D* hWithPIDMeasVsTruthMomentumTruDeuterons_;
  TH2D* hWithPIDTruthMomentumNotReco_;
  TH2D* hWithPIDMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponseWithPIDAllTrks_;
  TH2D* hWithPIDAllTrksTruthMomentum_;
  TH2D* hWithPIDAllTrksTruthMomentumReco_;
  TH2D* hWithPIDAllTrksMeasVsTruthMomentumTruProtons_;
  TH2D* hWithPIDAllTrksMeasVsTruthMomentumTruDeuterons_;
  TH2D* hWithPIDAllTrksTruthMomentumNotReco_;
  TH2D* hWithPIDAllTrksMeasuredMomentum_;

};

#endif/*HistMuCapRooUnfold_h*/
