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
HistMuCapRooUnfold::HistUnfold1D::HistUnfold1D(HistogramFactory &H, std::string Dir, std::string Name, bool MCTruth, int NBinP, double MaxP){
  hMeasuredMomentum_ = H.DefineTH1D(Dir+"/LateResponse"+Name, "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NBinP, 0., MaxP);

  if(MCTruth){
    Response_.Setup(NBinP, 0., MaxP, NBinP, 0., MaxP);
    H.Store(&Response_, "anDnLateResponse", Dir);

    hTruthMomentum_        = H.DefineTH1D(Dir+"/LateResponse"+Name, "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NBinP, 0., MaxP);
    hTruthMomentumReco_    = H.DefineTH1D(Dir+"/LateResponse"+Name, "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NBinP, 0., MaxP);
    hMeasVsTruthMomentum_  = H.DefineTH2D(Dir+"/LateResponse"+Name, "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NBinP, 0., MaxP,NBinP, 0., MaxP);
    hTruthMomentumNotReco_ = H.DefineTH1D(Dir+"/LateResponse"+Name, "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NBinP, 0., MaxP);
  }
}
  

//================================================================
void HistMuCapRooUnfold::HistUnfold1D::Reset(){
  Selected = false;
}

//================================================================
void HistMuCapRooUnfold::HistUnfold1D::FillMeasured(double mom){
  hMeasuredMomentum_->Fill(mom);
}

//================================================================
void HistMuCapRooUnfold::HistUnfold1D::FillTruth(double reco, double tru){
  Selected = true;
  Response_.Fill(reco, tru);
  hTruthMomentum_->Fill(tru);
  hTruthMomentumReco_->Fill(reco);
  hMeasVsTruthMomentum_->Fill(tru, reco);
}

//================================================================
void HistMuCapRooUnfold::HistUnfold1D::MissTruth(double tru){
  if (Selected)
    return;
  Response_.Miss(tru);
  hTruthMomentum_->Fill(tru);
  hTruthMomentumReco_->Fill(tru);
}

//****************************************************************


//================================================================
void HistMuCapRooUnfold::init(HistogramFactory& H,
                                     const std::string& hdir,
                                     const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  //----------------------------------------------------------------
  int NbBinP = 30;
  double MaxP = 300.;

  FullSpectrum_ = new HistUnfold1D(H, hdir, "", doMCTruth_, NbBinP, MaxP);
  Contained_ = new HistUnfold1D(H, hdir, "Contained", doMCTruth_, NbBinP, MaxP);
  PlnRngCutPln_ = new HistUnfold1D(H, hdir, "PlnRngCutPln", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone1_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone1", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone2_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone2", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone3_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone3", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone4_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone4", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone5_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone5", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone6_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone6", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone1AllTrks_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone1AllTrks", doMCTruth_, NbBinP, MaxP);
  PlnVsPCutZone2AllTrks_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone2AllTrks", doMCTruth_, NbBinP, MaxP);

  hWithPIDMeasuredMomentum_ = H.DefineTH2D(hdir+"/LateResponseWithPID", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]", 2,0,2.,NbBinP, 0., MaxP);


  //----------------------------------------------------------------
  if(doMCTruth_) {

    // Temporary histos just to define the 2D response functions
    TH2D *MeasuredTmp = new TH2D("MeasuredMomentumVsPID", "Measured momentum vs PID;PID;Momentum", 2,0,2., NbBinP, 0., MaxP);
    TH2D *TrueTmp = new TH2D("TrueMomentumVsPID", "True momentum vs PID;PID;Momentum", 2,0,2., NbBinP, 0., MaxP);
    anDnLateResponseWithPID_.Setup(MeasuredTmp, TrueTmp);
    H.Store(&anDnLateResponseWithPID_, "anDnLateResponseWithPID", hdir);
    anDnLateResponseWithPIDAllTrks_.Setup(MeasuredTmp, TrueTmp);
    H.Store(&anDnLateResponseWithPIDAllTrks_, "anDnLateResponseWithPIDAllTrks", hdir);

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

  anDnLateResponseWithPIDAllTrks_Filled_ = false;
  anDnLateResponseWithPID_Filled_ = false;
}

//================================================================
void HistMuCapRooUnfold::Fill(const EventClass& evt, int iPosTrack, int iNegTrack, bool isPosTrackContained, double rangePIDVar) {

  if(iNegTrack == -1) { // Veto DIO events
    if(iPosTrack != -1) { // Got a reconstructed capture track
      FullSpectrum_->FillMeasured(evt.ptot[iPosTrack]);
      double trackEnd = double(evt.hefit_pstop[iPosTrack]);
      int RecoPIDProton = int(double(trackEnd-28) < (0.40 * evt.ptot[iPosTrack] - 22.));
      if(isPosTrackContained) {
        // The "contained tracks" analysis channel
        Contained_->FillMeasured(evt.ptot[iPosTrack]);
        hWithPIDMeasuredMomentum_->Fill(RecoPIDProton, evt.ptot[iPosTrack]);
        if( (trackEnd-28) > 10){
          PlnRngCutPln_->FillMeasured(evt.ptot[iPosTrack]);
        }
        if ( RecoPIDProton == 0 ){
          if( (trackEnd-28) > 10){
            // Zone 1
            PlnVsPCutZone1_->FillMeasured(evt.ptot[iPosTrack]);
          } else {
            // Zone 4
            PlnVsPCutZone4_->FillMeasured(evt.ptot[iPosTrack]);
          }
        } else {
          if( (trackEnd-28) > 14){
            // Zone 2
            PlnVsPCutZone2_->FillMeasured(evt.ptot[iPosTrack]);
          } else {
            // Zone 3
            PlnVsPCutZone3_->FillMeasured(evt.ptot[iPosTrack]);
          }
        }
      }
      if ( RecoPIDProton == 0 ){
        if( (trackEnd-28) > 10){
          PlnVsPCutZone1AllTrks_->FillMeasured(evt.ptot[iPosTrack]);
        }
      } else {
        if( (trackEnd-28) > 14){
          PlnVsPCutZone2AllTrks_->FillMeasured(evt.ptot[iPosTrack]);
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
        FullSpectrum_->FillTruth(evt.ptot[iPosTrack], p_true_);

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

          Contained_->FillTruth(evt.ptot[iPosTrack], p_true_);
          if ( (trackEnd-28) > 10) {
            PlnRngCutPln_->FillTruth(evt.ptot[iPosTrack], p_true_);
          }
          if ( RecoPIDProton == 0 ){
            if ( IsInTrkRangeProton) {
              // Zone 1
              PlnVsPCutZone1_->FillTruth(evt.ptot[iPosTrack], p_true_);
            } else {
              // Zone 4
              PlnVsPCutZone4_->FillTruth(evt.ptot[iPosTrack], p_true_);
            }
          } else {
            if ( IsInTrkRangeDeuteron ){
              // Zone 2
              PlnVsPCutZone2_->FillTruth(evt.ptot[iPosTrack], p_true_);
            } else {
              // Zone 3
              PlnVsPCutZone3_->FillTruth(evt.ptot[iPosTrack], p_true_);
            }
          }
        } else {
          // uncontained tracks
          if ( RecoPIDProton == 0 ){
            if ( IsInTrkRangeProton) {
              // Zone 5
              PlnVsPCutZone5_->FillTruth(evt.ptot[iPosTrack], p_true_);
            }
          } else {
            if ( IsInTrkRangeDeuteron ){
              // Zone 6
              PlnVsPCutZone6_->FillTruth(evt.ptot[iPosTrack], p_true_);
            }
          }
        } 
        // contained + uncontained tracks
        if ( RecoPIDProton == 0 ){
          if ( IsInTrkRangeProton) {
            // Zone 1
            PlnVsPCutZone1AllTrks_->FillTruth(evt.ptot[iPosTrack], p_true_);
          }
        } else {
          if ( IsInTrkRangeDeuteron ){
            // Zone 2
            PlnVsPCutZone2AllTrks_->FillTruth(evt.ptot[iPosTrack], p_true_);
          }
        }
      }
    }
  }

}

//================================================================
void HistMuCapRooUnfold::FillAndMiss(){
  if(doMCTruth_ && truCaptEvt_) {

    FullSpectrum_->MissTruth(p_true_);


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
    Contained_->MissTruth(p_true_);
    PlnRngCutPln_->MissTruth(p_true_);
    PlnVsPCutZone1_->MissTruth(p_true_);
    PlnVsPCutZone4_->MissTruth(p_true_);
    PlnVsPCutZone2_->MissTruth(p_true_);
    PlnVsPCutZone3_->MissTruth(p_true_);
    PlnVsPCutZone5_->MissTruth(p_true_);
    PlnVsPCutZone6_->MissTruth(p_true_);
    PlnVsPCutZone1AllTrks_->MissTruth(p_true_);
    PlnVsPCutZone2AllTrks_->MissTruth(p_true_);

  }
}
