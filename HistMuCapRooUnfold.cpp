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
  hMeasuredMomentum_ = H.DefineTH1D(Dir+"/rooUnfold/"+Name, "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",NBinP, 0., MaxP);

  if(MCTruth){
    Response_.Setup(NBinP, 0., MaxP, NBinP, 0., MaxP);
    H.Store(&Response_, "unfoldMatrix"+Name, Dir+"/rooUnfold");

    hTruthMomentum_        = H.DefineTH1D(Dir+"/rooUnfold/"+Name, "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",NBinP, 0., MaxP);
    hTruthMomentumReco_    = H.DefineTH1D(Dir+"/rooUnfold/"+Name, "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",NBinP, 0., MaxP);
    hMeasVsTruthMomentum_  = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",NBinP, 0., MaxP,NBinP, 0., MaxP);
    hTruthMomentumNotReco_ = H.DefineTH1D(Dir+"/rooUnfold/"+Name, "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",NBinP, 0., MaxP);
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
  hTruthMomentumNotReco_->Fill(tru);
}

//****************************************************************


//================================================================
HistMuCapRooUnfold::HistUnfold2D::HistUnfold2D(HistogramFactory &H, std::string Dir, std::string Name, bool MCTruth, int NBinP, double MaxP){
  hMeasuredMomentum_  = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MeasuredMomentum", "Measured momentum spectrum;PID;Momentum [MeV/c]",2,0.,2.,NBinP, 0., MaxP);

  if(MCTruth){
    std::string tmpStr = "MeasuredMomentumVsPID"+Name;
    TH2D *MeasuredTmp = new TH2D(tmpStr.c_str(), "Measured momentum vs PID;PID;Momentum", 2,0,2., NBinP, 0., MaxP);
    tmpStr = "TrueMomentumVsPID"+Name;
    TH2D *TrueTmp     = new TH2D(tmpStr.c_str(), "True momentum vs PID;PID;Momentum", 2,0,2., NBinP, 0., MaxP);
    Response_.Setup(MeasuredTmp, TrueTmp);
    H.Store(&Response_, "unfoldMatrix"+Name, Dir+"/rooUnfold");

    hTruthMomentum_        = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MCTruthMomentum", "True momentum used in response function;PID;Momentum [MeV/c]",2,0,2., NBinP, 0., MaxP);
    hTruthMomentumReco_    = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MCTruthMomentumReco", "True momentum of reconstructed tracks;PID;Momentum [MeV/c]",2,0,2., NBinP, 0., MaxP);
    hMeasVsTruthMomentumTruProtons_   = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MCMeasVsTruthMomentumTruProtons", "Measured vs. true momentum used in response function for tru protons;True momentum [MeV/c];Measured momentum [MeV/c]",NBinP, 0., MaxP,NBinP, 0., MaxP);
    hMeasVsTruthMomentumTruDeuterons_ = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MCMeasVsTruthMomentumTruDeuterons", "Measured vs. true momentum used in response function for tru deuterons;True momentum [MeV/c];Measured momentum [MeV/c]",NBinP, 0., MaxP,NBinP, 0., MaxP);
    hTruthMomentumNotReco_ = H.DefineTH2D(Dir+"/rooUnfold/"+Name, "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;PID;Momentum [MeV/c]",2,0,2., NBinP, 0., MaxP);
  }
}

//================================================================
void HistMuCapRooUnfold::HistUnfold2D::Reset(){
  Selected = false;
}

//================================================================
void HistMuCapRooUnfold::HistUnfold2D::FillMeasured(int PID, double mom){
  hMeasuredMomentum_->Fill(PID, mom);
}

//================================================================
void HistMuCapRooUnfold::HistUnfold2D::FillTruth(int recoPID, int truPID, double reco, double tru){
  Selected = true;
  Response_.Fill(recoPID, reco, truPID, tru);
  hTruthMomentum_->Fill(truPID, tru);
  hTruthMomentumReco_->Fill(recoPID, reco);
  if (truPID == 0){
    hMeasVsTruthMomentumTruProtons_->Fill(tru,reco);
  } else {                        
    hMeasVsTruthMomentumTruDeuterons_->Fill(tru,reco);
  }
}

//================================================================
void HistMuCapRooUnfold::HistUnfold2D::MissTruth(int truPID, double tru){
  if (Selected)
    return;
  Response_.Miss(truPID, tru);
  hTruthMomentum_->Fill(truPID, tru);
  hTruthMomentumNotReco_->Fill(truPID, tru);
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
  int N1D = 0;

  FullSpectrum_   = new HistUnfold1D(H, hdir, "", doMCTruth_, NbBinP, MaxP);               AllUnfold1D_[N1D++] = FullSpectrum_;
  Contained_      = new HistUnfold1D(H, hdir, "Contained", doMCTruth_, NbBinP, MaxP);      AllUnfold1D_[N1D++] = Contained_;
  PlnRngCutPln_   = new HistUnfold1D(H, hdir, "PlnRngCutPln", doMCTruth_, NbBinP, MaxP);   AllUnfold1D_[N1D++] = PlnRngCutPln_;
  PlnVsPCutZone1_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone1", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone1_;
  PlnVsPCutZone2_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone2", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone2_;
  PlnVsPCutZone3_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone3", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone3_;
  PlnVsPCutZone4_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone4", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone4_;
  PlnVsPCutZone5_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone5", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone5_;
  PlnVsPCutZone6_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone6", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone6_;
  PlnVsPCutZone1AllTrks_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone1AllTrks", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone1AllTrks_;
  PlnVsPCutZone2AllTrks_ = new HistUnfold1D(H, hdir, "PlnVsPCutZone2AllTrks", doMCTruth_, NbBinP, MaxP); AllUnfold1D_[N1D++] = PlnVsPCutZone2AllTrks_;

  WithPID_ = new HistUnfold2D(H, hdir, "WithPID", doMCTruth_, NbBinP, MaxP);
  WithPIDAllTrks_ = new HistUnfold2D(H, hdir, "WithPIDAllTrks", doMCTruth_, NbBinP, MaxP);


}


void HistMuCapRooUnfold::SaveEventVariables(const EventClass& evt) {
  truCaptEvt_ = false;
  if(doMCTruth_) {
    p_true_ = 0.0;
    const int imctrk = evt.iCaptureMcTrk;
    TruePIDProton_ = -1;
    if ( imctrk > -1 ){
      p_true_ = evt.mcvertex_ptot[evt.iCaptureMcVtxStart];
      if ( evt.mctrack_pid[imctrk] == MuCapUtilities::PID_G3_PROTON){
        TruePIDProton_ = 0;
      } else if ( evt.mctrack_pid[imctrk] == MuCapUtilities::PID_G3_DEUTERON){
        TruePIDProton_ = 1;
      }
    }
    truCaptEvt_ = evt.iCaptureMcVtxStart != -1;
  }

  WithPID_->Reset();
  WithPIDAllTrks_->Reset();

  for (int n = 0; n < NBUNFOLD1D; n++)
    AllUnfold1D_[n]->Reset();

}

//================================================================
void HistMuCapRooUnfold::FillMeasured(const EventClass& evt, int iPosTrack, int iNegTrack, bool isPosTrackContained, double rangePIDVar) {

  if(iNegTrack == -1) { // Veto DIO events
    if(iPosTrack != -1) { // Got a reconstructed capture track
      FullSpectrum_->FillMeasured(evt.ptot[iPosTrack]);
      double trackEnd = double(evt.hefit_pstop[iPosTrack]);
      int RecoPIDProton = int(double(trackEnd-28) < (0.40 * evt.ptot[iPosTrack] - 22.));
      WithPIDAllTrks_->FillMeasured(RecoPIDProton, evt.ptot[iPosTrack]);
      if(isPosTrackContained) {
        // The "contained tracks" analysis channel
        Contained_->FillMeasured(evt.ptot[iPosTrack]);
        WithPID_->FillMeasured(RecoPIDProton, evt.ptot[iPosTrack]);
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
}

void HistMuCapRooUnfold::Fill(const EventClass& evt, int iPosTrack, int iNegTrack, bool isPosTrackContained, double rangePIDVar) {
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

        WithPIDAllTrks_->FillTruth(RecoPIDProton, TruePIDProton_, evt.ptot[iPosTrack], p_true_);
        if (IsContained ) {
          WithPID_->FillTruth(RecoPIDProton, TruePIDProton_, evt.ptot[iPosTrack], p_true_);

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

    WithPID_->MissTruth(TruePIDProton_, p_true_);
    WithPIDAllTrks_->MissTruth(TruePIDProton_, p_true_);

    for (int n = 0; n < NBUNFOLD1D; n++)
      AllUnfold1D_[n]->MissTruth(p_true_);

  }
}
