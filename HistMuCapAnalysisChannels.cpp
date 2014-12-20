// Andrei Gaponenko, 2014

#include "HistMuCapAnalysisChannels.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapAnalysisChannels::init(HistogramFactory& hf,
                                     const std::string& hdir,
                                     const DetectorGeo& geom,
                                     const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  //----------------------------------------------------------------
  // "channel" analysis histograms
  const int recoContPNbins = 88;
  const double recoContPMin = 30.;
  const double recoContPMax = 250.;
  const int recoContRNbins = 6;
  const double recoContRMin = 8.;
  const double recoContRMax = 32.;

  const int recoUncPNbins = 88;
  const double recoUncPMin = 30.;
  const double recoUncPMax = 250.;

  //----------------------------------------------------------------
  contained_prange_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP", "Last plane hit/|cos(theta)| vs p, contained",
                                    recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
  contained_prange_->SetOption("colz");
  contained_prange_->GetXaxis()->SetTitle("p [MeV/c]");
  contained_prange_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

  if(doMCTruth_) {
    contained_prange_mcproton_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP_mcproton", "Last plane hit/|cos(theta)| vs p, contained mcproton",
                                               recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
    contained_prange_mcproton_->SetOption("colz");
    contained_prange_mcproton_->GetXaxis()->SetTitle("p [MeV/c]");
    contained_prange_mcproton_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

    contained_prange_mcdeuteron_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP_mcdeuteron", "Last plane hit/|cos(theta)| vs p, contained mcdeuteron",
                                                 recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
    contained_prange_mcdeuteron_->SetOption("colz");
    contained_prange_mcdeuteron_->GetXaxis()->SetTitle("p [MeV/c]");
    contained_prange_mcdeuteron_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

    contained_prange_mcdio_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP_mcdio", "Last plane hit/|cos(theta)| vs p, contained mcdio",
                                            recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
    contained_prange_mcdio_->SetOption("colz");
    contained_prange_mcdio_->GetXaxis()->SetTitle("p [MeV/c]");
    contained_prange_mcdio_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");
  }

  //----------------------------------------------------------------
  uncontained_p_ = hf.DefineTH1D(hdir+"/uncontained", "ptot", "ptot, non contained", recoUncPNbins, recoUncPMin, recoUncPMax);
  uncontained_p_->GetXaxis()->SetTitle("p, MeV/c");

  if(doMCTruth_) {
    uncontained_p_mcproton_ = hf.DefineTH1D(hdir+"/uncontained", "ptot_mcproton", "ptot, non contained mcproton", recoUncPNbins, recoUncPMin, recoUncPMax);
    uncontained_p_mcproton_->GetXaxis()->SetTitle("p, MeV/c");

    uncontained_p_mcdeuteron_ = hf.DefineTH1D(hdir+"/uncontained", "ptot_mcdeuteron", "ptot, non contained mcdeuteron", recoUncPNbins, recoUncPMin, recoUncPMax);
    uncontained_p_mcdeuteron_->GetXaxis()->SetTitle("p, MeV/c");

    uncontained_p_mcdio_ = hf.DefineTH1D(hdir+"/uncontained", "ptot_mcdio", "ptot, non contained mcdio", recoUncPNbins, recoUncPMin, recoUncPMax);
    uncontained_p_mcdio_->GetXaxis()->SetTitle("p, MeV/c");
  }

  //----------------------------------------------------------------
  hitbased_.init(hf, hdir+"/hitbased", geom, conf);

  hTDCWidthContained_.init(hf, hdir+"/contained/tdcwidth", geom, conf);
  hTDCWidthUncontained_.init(hf, hdir+"/uncontained/tdcwidth", geom, conf);

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning
    const int gen1nbins = 400; // 2.5 MeV/c bins
    const double gen1pmin = 0.;
    const double gen1pmax = 400.;

    // True distribution of lost events
    mclost2_ptot_ = hf.DefineTH1D(hdir, "mclost2_ptot", "mcptot, of lost events for 2 channel analysis", gen1nbins, gen1pmin, gen1pmax);
    mclost2_count_ = hf.DefineTH1D(hdir, "mclost2_count", "count of lost events for 2 channel analysis ", 1, -0.5, 0.5);

    mclost3_ptot_ = hf.DefineTH1D(hdir, "mclost3_ptot", "mcptot, of lost events for 2 channel analysis", gen1nbins, gen1pmin, gen1pmax);
    mclost3_count_ = hf.DefineTH1D(hdir, "mclost3_count", "count of lost events for 2 channel analysis ", 1, -0.5, 0.5);

    // True distribution of all events used for the unfolding
    mcin_proton_ptot_ = hf.DefineTH1D(hdir, "mcin_proton_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcin_deuteron_ptot_ = hf.DefineTH1D(hdir, "mcin_deuteron_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcin_dio_count_ = hf.DefineTH1D(hdir, "mcin_dio_count", "noncapture count, input", 1, -0.5, 0.5);

    // Migration matrices for the contained channel, two generator binnings
    containedMigration_ = hf.DefineTH3D(hdir+"/contained", "migration",
                                        "Contained channel migration",
                                        gen1nbins, gen1pmin, gen1pmax,
                                        recoContPNbins, recoContPMin, recoContPMax,
                                        recoContRNbins, recoContRMin, recoContRMax);

    containedMigration_->GetXaxis()->SetTitle("p true, MeV/c");
    containedMigration_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedMigration_->GetZaxis()->SetTitle("range");

    // Migration matrices for the non contained channel, two generator binnings
    uncontainedMigration_ = hf.DefineTH2D(hdir+"/uncontained","migration",
                                          "Non-contained channel migration",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedMigration_->SetOption("colz");
    uncontainedMigration_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedMigration_->GetYaxis()->SetTitle("p reco, MeV/c");

    hTruthTrkContained_.init(hf, hdir+"/contained/MCTruthTrk", conf);
    hTruthTrkUncontained_.init(hf, hdir+"/uncontained/MCTruthTrk", conf);

    hResolutionContained_.init(hf, hdir+"/contained/resolution", conf);
    hResolutionUncontained_.init(hf, hdir+"/uncontained/resolution", conf);
  }

}

//================================================================
void HistMuCapAnalysisChannels::fill(const EventClass& evt,
                                     int iPosTrack,
                                     int iNegTrack,
                                     bool isPosTrackContained,
                                     double rangePIDVar,
                                     const ClustersByPlane& globalPlaneClusters
                                     )
{
  //----------------------------------------------------------------
  // Figure out an exclusive analysis channel for this event

  bool eventUsedInAChannel2 = false;
  bool eventUsedInHitBased = false;

  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  // Simulated DIO have no easily accessible MC truth.  We'll tread PID=zero as DIO down in this code.
  const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

  if(iNegTrack == -1) { // Veto DIO events
    if(iPosTrack != -1) { // Got a reconstructed capture track

      const double prec = evt.ptot[iPosTrack];

      if(isPosTrackContained) {
        // The "contained tracks" analysis channel
        eventUsedInAChannel2 = true;
        contained_prange_->Fill(prec, rangePIDVar);
        hTDCWidthContained_.fill(evt, globalPlaneClusters);

        if(doMCTruth_) {

          if(imcvtxStart  != -1) {
            containedMigration_->Fill(evt.mcvertex_ptot[imcvtxStart], prec, rangePIDVar);
          }

          switch(mcParticle) {
          case MuCapUtilities::PID_G3_PROTON:
            contained_prange_mcproton_->Fill(prec, rangePIDVar);
            break;
          case MuCapUtilities::PID_G3_DEUTERON:
            contained_prange_mcdeuteron_->Fill(prec, rangePIDVar);
            break;
          case 0:
            contained_prange_mcdio_->Fill(prec, rangePIDVar);
            break;
          }

          hTruthTrkContained_.fill(evt);
          hResolutionContained_.fill(evt, iPosTrack);
        }
      }
      else { // The non-contained tracks channel
        eventUsedInAChannel2 = true;
        uncontained_p_->Fill(prec);

        hTDCWidthUncontained_.fill(evt, globalPlaneClusters);

        if(doMCTruth_) {

          if(imcvtxStart  != -1) {
            uncontainedMigration_->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
          }

          switch(mcParticle) {
          case MuCapUtilities::PID_G3_PROTON:
            uncontained_p_mcproton_->Fill(prec);
            break;
          case MuCapUtilities::PID_G3_DEUTERON:
            uncontained_p_mcdeuteron_->Fill(prec);
            break;
          case 0:
            uncontained_p_mcdio_->Fill(prec);
            break;
          }

          hTruthTrkUncontained_.fill(evt);
          hResolutionUncontained_.fill(evt, iPosTrack);
        }
      }
    }
  }

  if(!eventUsedInAChannel2) {
    // No good positive tracks. Do a hit-based analysis here.
    eventUsedInHitBased = hitbased_.accepted(evt, globalPlaneClusters, iNegTrack);
  }

  if(doMCTruth_) {
    // Truth momentum with the binning used in the unfolding
    // Keep protons and deuterons separately to compare hadd-ed
    // pseudodata truth to unfolding results.
    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      mcin_proton_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      mcin_deuteron_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    case 0:
      mcin_dio_count_->Fill(0.);
      break;
    }

    if(!eventUsedInAChannel2) {
      mclost2_count_->Fill(0.);
      if(imcvtxStart  != -1) {
        mclost2_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      }

      if(!eventUsedInHitBased) {
        mclost3_count_->Fill(0.);
        if(imcvtxStart  != -1) {
          mclost3_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
        }
      }
    }
  }
}

//================================================================
