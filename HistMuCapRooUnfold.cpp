// Andrei Gaponenko, 2014

#include "HistMuCapRooUnfold.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapRooUnfold::init(HistogramFactory& H,
                                     const std::string& hdir,
                                     const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  //----------------------------------------------------------------
  int NbBinP = 30;
  double MaxP = 300.;

  hMeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponse", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);
  hContainedMeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponseContained", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);
  hPlnRngCutPlnMeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnRngCutPln", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);
  hPlnVsPCutZone1MeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);
  hPlnVsPCutZone2MeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);
  hPlnVsPCutZone3MeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone3", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);
  hPlnVsPCutZone4MeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone4", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NbBinP, 0., MaxP);

  hWithPIDMeasuredMomentum_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]", 2,0,2.,NbBinP, 0., MaxP);


  //----------------------------------------------------------------
  if(doMCTruth_) {
    anDnLateResponse_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponse_, "anDnLateResponse", hdir);
    anDnLateResponseContained_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponseContained_, "anDnLateResponseContained", hdir);
    anDnLateResponsePlnRngCutPln_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnRngCutPln_, "anDnLateResponsePlnRngCutPln", hdir);

    anDnLateResponsePlnVsPCutZone1_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone1_, "anDnLateResponsePlnVsPCutZone1", hdir);
    anDnLateResponsePlnVsPCutZone2_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone2_, "anDnLateResponsePlnVsPCutZone2", hdir);
    anDnLateResponsePlnVsPCutZone3_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone3_, "anDnLateResponsePlnVsPCutZone3", hdir);
    anDnLateResponsePlnVsPCutZone4_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone4_, "anDnLateResponsePlnVsPCutZone4", hdir);

    // Contained + uncontained tracks
    anDnLateResponsePlnVsPCutZone1AllTrks_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone1AllTrks_, "anDnLateResponsePlnVsPCutZone1AllTrks", hdir);
    anDnLateResponsePlnVsPCutZone2AllTrks_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone2AllTrks_, "anDnLateResponsePlnVsPCutZone2AllTrks", hdir);

    // Contained + tracks
    anDnLateResponsePlnVsPCutZone5_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone5_, "anDnLateResponsePlnVsPCutZone5", hdir);
    anDnLateResponsePlnVsPCutZone6_.Setup(NbBinP, 0., MaxP, NbBinP, 0., MaxP);
    H.Store(&anDnLateResponsePlnVsPCutZone6_, "anDnLateResponsePlnVsPCutZone6", hdir);

    // Temporary histos just to define the 2D response functions
    TH2D *MeasuredTmp = new TH2D("MeasuredMomentumVsPID", "Measured momentum vs PID;PID;Momentum", 2,0,2., NbBinP, 0., MaxP);
    TH2D *TrueTmp = new TH2D("TrueMomentumVsPID", "True momentum vs PID;PID;Momentum", 2,0,2., NbBinP, 0., MaxP);
    anDnLateResponseWithPID_.Setup(MeasuredTmp, TrueTmp);
    H.Store(&anDnLateResponseWithPID_, "anDnLateResponseWithPID", hdir);
    anDnLateResponseWithPIDAllTrks_.Setup(MeasuredTmp, TrueTmp);
    H.Store(&anDnLateResponseWithPIDAllTrks_, "anDnLateResponseWithPIDAllTrks", hdir);

    hTruthMomentum_ = H.DefineTH1D(hdir+"/LateResponse", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hTruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponse", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hMeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponse", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hTruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponse", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hContainedTruthMomentum_ = H.DefineTH1D(hdir+"/LateResponseContained", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hContainedTruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponseContained", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hContainedMeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponseContained", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hContainedTruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponseContained", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnRngCutPlnTruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnRngCutPln", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnRngCutPlnTruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnRngCutPln", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnRngCutPlnMeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnRngCutPln", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnRngCutPlnTruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnRngCutPln", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone1TruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone1TruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone1MeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone1", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone1TruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone2TruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone2TruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone2MeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone2", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone2TruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone1AllTrksTruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1AllTrks", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone1AllTrksTruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1AllTrks", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone1AllTrksMeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone1AllTrks", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone1AllTrksTruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone1AllTrks", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone2AllTrksTruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2AllTrks", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone2AllTrksTruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2AllTrks", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone2AllTrksMeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone2AllTrks", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone2AllTrksTruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone2AllTrks", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone5TruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone5", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone5TruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone5", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone5MeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone5", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone5TruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone5", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone6TruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone6", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone6TruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone6", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone6MeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone6", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone6TruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone6", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone3TruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone3", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone3TruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone3", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone3MeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone3", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone3TruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone3", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hPlnVsPCutZone4TruthMomentum_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone4", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone4TruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone4", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NbBinP, 0., MaxP);
    hPlnVsPCutZone4MeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponsePlnVsPCutZone4", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hPlnVsPCutZone4TruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponsePlnVsPCutZone4", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NbBinP, 0., MaxP);

    hWithPIDTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MCTruthMomentum", "True momentum used in response function, contained trks;PID;Momentum [MeV/c]",2,0,2., NbBinP, 0., MaxP);
    hWithPIDTruthMomentumReco_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MCTruthMomentumReco", "True momentum of reconstructed tracks, contained trks;PID;Momentum [MeV/c]",2,0,2., NbBinP, 0., MaxP);
    hWithPIDMeasVsTruthMomentumTruProtons_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MCMeasVsTruthMomentumTruProtons", "Measured vs. true momentum used in response function for tru protons, contained trks;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hWithPIDMeasVsTruthMomentumTruDeuterons_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MCMeasVsTruthMomentumTruDeuterons", "Measured vs. true momentum used in response function for tru deuterons, contained trks;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hWithPIDTruthMomentumNotReco_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed, contained trks;PID;Momentum [MeV/c]",2,0,2., NbBinP, 0., MaxP);

    hWithPIDAllTrksTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponseWithPIDAllTrks", "MCTruthMomentum", "True momentum used in response function;PID;Momentum [MeV/c]",2,0,2., NbBinP, 0., MaxP);
    hWithPIDAllTrksTruthMomentumReco_ = H.DefineTH2D(hdir+"/LateResponseWithPIDAllTrks", "MCTruthMomentumReco", "True momentum of reconstructed tracks;PID;Momentum [MeV/c]",2,0,2., NbBinP, 0., MaxP);
    hWithPIDAllTrksMeasVsTruthMomentumTruProtons_ = H.DefineTH2D(hdir+"/LateResponseWithPIDAllTrks", "MCMeasVsTruthMomentumTruProtons", "Measured vs. true momentum used in response function for tru protons;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hWithPIDAllTrksMeasVsTruthMomentumTruDeuterons_ = H.DefineTH2D(hdir+"/LateResponseWithPIDAllTrks", "MCMeasVsTruthMomentumTruDeuterons", "Measured vs. true momentum used in response function for tru deuterons;True momentum [MeV/c];Measured momentum [MeV/c]",NbBinP, 0., MaxP,NbBinP, 0., MaxP);
    hWithPIDAllTrksTruthMomentumNotReco_ = H.DefineTH2D(hdir+"/LateResponseWithPIDAllTrks", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;PID;Momentum [MeV/c]",2,0,2., NbBinP, 0., MaxP);
  }

}


void HistMuCapRooUnfold::SaveEventVariables(const EventClass& evt) {
  truCaptEvt_ = false;
  if(doMCTruth_) {

    const int imctrk = evt.iCaptureMcTrk;
    p_true_ = evt.mcvertex_ptot[evt.iCaptureMcVtxStart];
    TruePIDProton_ = -1;
    if ( imctrk > -1 ){
      if ( evt.mctrack_pid[imctrk] == MuCapUtilities::PID_G3_PROTON){
        TruePIDProton_ = 0;
      } else if ( evt.mctrack_pid[imctrk] == MuCapUtilities::PID_G3_DEUTERON){
        TruePIDProton_ = 1;
      }
    }
    truCaptEvt_ = evt.iCaptureMcVtxStart != -1;
  }

  anDnLateResponse_Filled_ = false;
  anDnLateResponseWithPIDAllTrks_Filled_ = false;
  anDnLateResponseWithPID_Filled_ = false;
  anDnLateResponseContained_Filled_ = false;
  anDnLateResponsePlnRngCutPln_Filled_ = false;
  anDnLateResponsePlnVsPCutZone1_Filled_ = false;
  anDnLateResponsePlnVsPCutZone4_Filled_ = false;
  anDnLateResponsePlnVsPCutZone2_Filled_ = false;
  anDnLateResponsePlnVsPCutZone3_Filled_ = false;
  anDnLateResponsePlnVsPCutZone5_Filled_ = false;
  anDnLateResponsePlnVsPCutZone6_Filled_ = false;
  anDnLateResponsePlnVsPCutZone1AllTrks_Filled_ = false;
  anDnLateResponsePlnVsPCutZone2AllTrks_Filled_ = false;
}

//================================================================
void HistMuCapRooUnfold::Fill(const EventClass& evt, int iPosTrack, int iNegTrack, bool isPosTrackContained, double rangePIDVar) {

  if(iNegTrack == -1) { // Veto DIO events
    if(iPosTrack != -1) { // Got a reconstructed capture track
      hMeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
      if(isPosTrackContained) {
        // The "contained tracks" analysis channel
        hContainedMeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
        double trackEnd = double(evt.hefit_pstop[iPosTrack]);
        int RecoPIDProton = int(double(trackEnd-28) < (0.40 * evt.ptot[iPosTrack] - 22.));
        hWithPIDMeasuredMomentum_->Fill(RecoPIDProton, evt.ptot[iPosTrack]);
        if( (trackEnd-28) > 10){
          hPlnRngCutPlnMeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
        }
        if ( RecoPIDProton == 0 ){
          if( (trackEnd-28) > 10){
            // Zone 1
            hPlnVsPCutZone1MeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
          } else {
            // Zone 4
            hPlnVsPCutZone4MeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
          }
        } else {
          if( (trackEnd-28) > 14){
            // Zone 2
            hPlnVsPCutZone2MeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
          } else {
            // Zone 3
            hPlnVsPCutZone3MeasuredMomentum_->Fill(evt.ptot[iPosTrack]);
          }
        }
      }
    }
  }

  if(doMCTruth_ && truCaptEvt_) {
    if(evt.iCaptureMcVtxStart != -1) {  // ANT: Should be removed ???
      bool IsContained = false;
      bool IsInTrkRange = false;
      bool IsInTrkRangeProton = false;
      bool IsInTrkRangeDeuteron = false;

      if( iNegTrack == -1 && iPosTrack != -1) {
        anDnLateResponse_Filled_ = true;
        anDnLateResponse_.Fill(evt.ptot[iPosTrack], p_true_);
        hTruthMomentum_->Fill(p_true_);
        hTruthMomentumReco_->Fill(p_true_);
        hMeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);

        double trackEnd = double(evt.hefit_pstop[iPosTrack]);
        IsContained = isPosTrackContained;
        IsInTrkRange = (trackEnd-28) > 10;
        IsInTrkRangeProton = (trackEnd-28) > 10;
        IsInTrkRangeDeuteron = (trackEnd-28) > 14;
        // 0 for protons, 1 for deuterons
        int RecoPIDProton = int(double(trackEnd-28) < (0.40 * evt.ptot[iPosTrack] - 22.));
        anDnLateResponseWithPIDAllTrks_Filled_ = true;
        anDnLateResponseWithPIDAllTrks_.Fill(RecoPIDProton, evt.ptot[iPosTrack], TruePIDProton_, p_true_);
        hWithPIDAllTrksTruthMomentum_->Fill(TruePIDProton_, p_true_);
        hWithPIDAllTrksTruthMomentumReco_->Fill(TruePIDProton_, p_true_);
        if (TruePIDProton_ == 0){
          hWithPIDAllTrksMeasVsTruthMomentumTruProtons_->Fill(p_true_,evt.ptot[iPosTrack]);
        } else {                        
          hWithPIDAllTrksMeasVsTruthMomentumTruDeuterons_->Fill(p_true_,evt.ptot[iPosTrack]);
        }
        if (IsContained ) {
          anDnLateResponseWithPID_Filled_ = true;
          anDnLateResponseWithPID_.Fill(RecoPIDProton, evt.ptot[iPosTrack], TruePIDProton_, p_true_);
          hWithPIDTruthMomentum_->Fill(TruePIDProton_, p_true_);
          hWithPIDTruthMomentumReco_->Fill(TruePIDProton_, p_true_);
          if (TruePIDProton_ == 0){
            hWithPIDMeasVsTruthMomentumTruProtons_->Fill(p_true_,evt.ptot[iPosTrack]);
          } else {                        
            hWithPIDMeasVsTruthMomentumTruDeuterons_->Fill(p_true_,evt.ptot[iPosTrack]);
          }

          anDnLateResponseContained_Filled_ = true;
          anDnLateResponseContained_.Fill(evt.ptot[iPosTrack], p_true_);
          hContainedTruthMomentum_->Fill(p_true_);
          hContainedTruthMomentumReco_->Fill(p_true_);
          hContainedMeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
          if ( (trackEnd-28) > 10) {
            anDnLateResponsePlnRngCutPln_Filled_ = true;
            anDnLateResponsePlnRngCutPln_.Fill(evt.ptot[iPosTrack], p_true_);
            hPlnRngCutPlnTruthMomentum_->Fill(p_true_);
            hPlnRngCutPlnTruthMomentumReco_->Fill(p_true_);
            hPlnRngCutPlnMeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
          }
          if ( RecoPIDProton == 0 ){
            if ( IsInTrkRangeProton) {
              // Zone 1
              anDnLateResponsePlnVsPCutZone1_Filled_ = true;
              anDnLateResponsePlnVsPCutZone1_.Fill(evt.ptot[iPosTrack], p_true_);
              hPlnVsPCutZone1TruthMomentum_->Fill(p_true_);
              hPlnVsPCutZone1TruthMomentumReco_->Fill(p_true_);
              hPlnVsPCutZone1MeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
            } else {
              // Zone 4
              anDnLateResponsePlnVsPCutZone4_Filled_ = true;
              anDnLateResponsePlnVsPCutZone4_.Fill(evt.ptot[iPosTrack], p_true_);
              hPlnVsPCutZone4TruthMomentum_->Fill(p_true_);
              hPlnVsPCutZone4TruthMomentumReco_->Fill(p_true_);
              hPlnVsPCutZone4MeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
            }
          } else {
            if ( IsInTrkRangeDeuteron ){
              // Zone 2
              anDnLateResponsePlnVsPCutZone2_Filled_ = true;
              anDnLateResponsePlnVsPCutZone2_.Fill(evt.ptot[iPosTrack], p_true_);
              hPlnVsPCutZone2TruthMomentum_->Fill(p_true_);
              hPlnVsPCutZone2TruthMomentumReco_->Fill(p_true_);
              hPlnVsPCutZone2MeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
            } else {
              // Zone 3
              anDnLateResponsePlnVsPCutZone3_Filled_ = true;
              anDnLateResponsePlnVsPCutZone3_.Fill(evt.ptot[iPosTrack], p_true_);
              hPlnVsPCutZone3TruthMomentum_->Fill(p_true_);
              hPlnVsPCutZone3TruthMomentumReco_->Fill(p_true_);
              hPlnVsPCutZone3MeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
            }
          }
        } else {
          // uncontained tracks
          if ( RecoPIDProton == 0 ){
            if ( IsInTrkRangeProton) {
              // Zone 1
              anDnLateResponsePlnVsPCutZone5_Filled_ = true;
              anDnLateResponsePlnVsPCutZone5_.Fill(evt.ptot[iPosTrack], p_true_);
              hPlnVsPCutZone5TruthMomentum_->Fill(p_true_);
              hPlnVsPCutZone5TruthMomentumReco_->Fill(p_true_);
              hPlnVsPCutZone5MeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
            }
          } else {
            if ( IsInTrkRangeDeuteron ){
              // Zone 2
              anDnLateResponsePlnVsPCutZone6_Filled_ = true;
              anDnLateResponsePlnVsPCutZone6_.Fill(evt.ptot[iPosTrack], p_true_);
              hPlnVsPCutZone6TruthMomentum_->Fill(p_true_);
              hPlnVsPCutZone6TruthMomentumReco_->Fill(p_true_);
              hPlnVsPCutZone6MeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
            }
          }
        } 
        // contained + uncontained tracks
        if ( RecoPIDProton == 0 ){
          if ( IsInTrkRangeProton) {
            // Zone 1
            anDnLateResponsePlnVsPCutZone1AllTrks_Filled_ = true;
            anDnLateResponsePlnVsPCutZone1AllTrks_.Fill(evt.ptot[iPosTrack], p_true_);
            hPlnVsPCutZone1AllTrksTruthMomentum_->Fill(p_true_);
            hPlnVsPCutZone1AllTrksTruthMomentumReco_->Fill(p_true_);
            hPlnVsPCutZone1AllTrksMeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
          }
        } else {
          if ( IsInTrkRangeDeuteron ){
            // Zone 2
            anDnLateResponsePlnVsPCutZone2AllTrks_Filled_ = true;
            anDnLateResponsePlnVsPCutZone2AllTrks_.Fill(evt.ptot[iPosTrack], p_true_);
            hPlnVsPCutZone2AllTrksTruthMomentum_->Fill(p_true_);
            hPlnVsPCutZone2AllTrksTruthMomentumReco_->Fill(p_true_);
            hPlnVsPCutZone2AllTrksMeasVsTruthMomentum_->Fill(p_true_,evt.ptot[iPosTrack]);
          }
        }
      }
    }
  }

}

//================================================================
void HistMuCapRooUnfold::FillAndMiss(){
  if(doMCTruth_ && truCaptEvt_) {

    if( !anDnLateResponse_Filled_){
      anDnLateResponse_.Miss( p_true_);
      hTruthMomentum_->Fill(p_true_);
      hTruthMomentumReco_->Fill(p_true_);
    }


    // If contained in false, then we haven't saved
    // the event yet and we must save it now.
    if( !anDnLateResponseWithPIDAllTrks_Filled_){
      anDnLateResponseWithPIDAllTrks_.Miss(TruePIDProton_, p_true_);
      hWithPIDAllTrksTruthMomentum_->Fill(TruePIDProton_, p_true_);
      hWithPIDAllTrksTruthMomentumNotReco_->Fill(TruePIDProton_, p_true_);
    }
    if( !anDnLateResponseWithPID_Filled_){
      anDnLateResponseWithPID_.Miss(TruePIDProton_, p_true_);
      hWithPIDTruthMomentum_->Fill(TruePIDProton_, p_true_);
      hWithPIDTruthMomentumNotReco_->Fill(TruePIDProton_, p_true_);
    }
    if( !anDnLateResponseContained_Filled_){
      anDnLateResponseContained_.Miss(p_true_);
      hContainedTruthMomentum_->Fill(p_true_);
      hContainedTruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnRngCutPln_Filled_){
      anDnLateResponsePlnRngCutPln_.Miss(p_true_);
      hPlnRngCutPlnTruthMomentum_->Fill(p_true_);
      hPlnRngCutPlnTruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone1_Filled_){
      anDnLateResponsePlnVsPCutZone1_.Miss(p_true_);
      hPlnVsPCutZone1TruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone1TruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone4_Filled_){
      anDnLateResponsePlnVsPCutZone4_.Miss(p_true_);
      hPlnVsPCutZone4TruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone4TruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone2_Filled_){
      anDnLateResponsePlnVsPCutZone2_.Miss(p_true_);
      hPlnVsPCutZone2TruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone2TruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone3_Filled_){
      anDnLateResponsePlnVsPCutZone3_.Miss(p_true_);
      hPlnVsPCutZone3TruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone3TruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone5_Filled_){
      anDnLateResponsePlnVsPCutZone5_.Miss(p_true_);
      hPlnVsPCutZone5TruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone5TruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone6_Filled_){
      anDnLateResponsePlnVsPCutZone6_.Miss(p_true_);
      hPlnVsPCutZone6TruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone6TruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone1AllTrks_Filled_){
      anDnLateResponsePlnVsPCutZone1AllTrks_.Miss(p_true_);
      hPlnVsPCutZone1AllTrksTruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone1AllTrksTruthMomentumNotReco_->Fill(p_true_);
    }
    if( !anDnLateResponsePlnVsPCutZone2AllTrks_Filled_){
      anDnLateResponsePlnVsPCutZone2AllTrks_.Miss(p_true_);
      hPlnVsPCutZone2AllTrksTruthMomentum_->Fill(p_true_);
      hPlnVsPCutZone2AllTrksTruthMomentumNotReco_->Fill(p_true_);
    }
  }
}
