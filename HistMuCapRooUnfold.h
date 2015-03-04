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

  bool anDnLateResponse_Filled_;

  bool anDnLateResponseWithPIDAllTrks_Filled_;
  bool anDnLateResponseWithPID_Filled_;
  bool anDnLateResponseContained_Filled_;
  bool anDnLateResponsePlnRngCutPln_Filled_;
  bool anDnLateResponsePlnVsPCutZone1_Filled_;
  bool anDnLateResponsePlnVsPCutZone4_Filled_;
  bool anDnLateResponsePlnVsPCutZone2_Filled_;
  bool anDnLateResponsePlnVsPCutZone3_Filled_;
  bool anDnLateResponsePlnVsPCutZone5_Filled_;
  bool anDnLateResponsePlnVsPCutZone6_Filled_;
  bool anDnLateResponsePlnVsPCutZone1AllTrks_Filled_;
  bool anDnLateResponsePlnVsPCutZone2AllTrks_Filled_;

  RooUnfoldResponse anDnLateResponse_;
  TH1D* hTruthMomentum_;
  TH1D* hTruthMomentumReco_;
  TH2D* hMeasVsTruthMomentum_;
  TH1D* hTruthMomentumNotReco_;
  TH1D* hMeasuredMomentum_;

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

  RooUnfoldResponse anDnLateResponseContained_;
  TH1D* hContainedTruthMomentum_;
  TH1D* hContainedTruthMomentumReco_;
  TH2D* hContainedMeasVsTruthMomentum_;
  TH1D* hContainedTruthMomentumNotReco_;
  TH1D* hContainedMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnRngCutPln_;
  TH1D* hPlnRngCutPlnTruthMomentum_;
  TH1D* hPlnRngCutPlnTruthMomentumReco_;
  TH2D* hPlnRngCutPlnMeasVsTruthMomentum_;
  TH1D* hPlnRngCutPlnTruthMomentumNotReco_;
  TH1D* hPlnRngCutPlnMeasuredMomentum_;

  // Contained tracks
  RooUnfoldResponse anDnLateResponsePlnVsPCutZone1_;
  TH1D* hPlnVsPCutZone1TruthMomentum_;
  TH1D* hPlnVsPCutZone1TruthMomentumReco_;
  TH2D* hPlnVsPCutZone1MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone1TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone1MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone2_;
  TH1D* hPlnVsPCutZone2TruthMomentum_;
  TH1D* hPlnVsPCutZone2TruthMomentumReco_;
  TH2D* hPlnVsPCutZone2MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone2TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone2MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone3_;
  TH1D* hPlnVsPCutZone3TruthMomentum_;
  TH1D* hPlnVsPCutZone3TruthMomentumReco_;
  TH2D* hPlnVsPCutZone3MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone3TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone3MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone4_;
  TH1D* hPlnVsPCutZone4TruthMomentum_;
  TH1D* hPlnVsPCutZone4TruthMomentumReco_;
  TH2D* hPlnVsPCutZone4MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone4TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone4MeasuredMomentum_;

  // Contained + uncontained tracks
  RooUnfoldResponse anDnLateResponsePlnVsPCutZone1AllTrks_;
  TH1D* hPlnVsPCutZone1AllTrksTruthMomentum_;
  TH1D* hPlnVsPCutZone1AllTrksTruthMomentumReco_;
  TH2D* hPlnVsPCutZone1AllTrksMeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone1AllTrksTruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone1AllTrksMeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone2AllTrks_;
  TH1D* hPlnVsPCutZone2AllTrksTruthMomentum_;
  TH1D* hPlnVsPCutZone2AllTrksTruthMomentumReco_;
  TH2D* hPlnVsPCutZone2AllTrksMeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone2AllTrksTruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone2AllTrksMeasuredMomentum_;

  // Uncontained tracks
  RooUnfoldResponse anDnLateResponsePlnVsPCutZone5_;
  TH1D* hPlnVsPCutZone5TruthMomentum_;
  TH1D* hPlnVsPCutZone5TruthMomentumReco_;
  TH2D* hPlnVsPCutZone5MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone5TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone5MeasuredMomentum_;

  RooUnfoldResponse anDnLateResponsePlnVsPCutZone6_;
  TH1D* hPlnVsPCutZone6TruthMomentum_;
  TH1D* hPlnVsPCutZone6TruthMomentumReco_;
  TH2D* hPlnVsPCutZone6MeasVsTruthMomentum_;
  TH1D* hPlnVsPCutZone6TruthMomentumNotReco_;
  TH1D* hPlnVsPCutZone6MeasuredMomentum_;

};

#endif/*HistMuCapRooUnfold_h*/
