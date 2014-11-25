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
                                     const DetectorGeo&,
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

  contained_prange_ = hf.DefineTH2D(hdir, "containedRangecosVsP", "Last plane hit/|cos(theta)| vs p, contained",
                                    recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
  contained_prange_->SetOption("colz");
  contained_prange_->GetXaxis()->SetTitle("p [MeV/c]");
  contained_prange_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

  uncontained_p_ = hf.DefineTH1D(hdir, "uncontained_ptot", "ptot, non contained", recoUncPNbins, recoUncPMin, recoUncPMax);
  uncontained_p_->GetXaxis()->SetTitle("p, MeV/c");

  if(doMCTruth_) {
    // truth level binning
    const int gen1nbins = 400; // 2.5 MeV/c bins
    const double gen1pmin = 0.;
    const double gen1pmax = 400.;

    // True distribution of lost events
    lost1_ptot_ = hf.DefineTH1D(hdir, "lost1_ptot", "mcptot, lost", gen1nbins, gen1pmin, gen1pmax);
    noncapture_lost_ = hf.DefineTH1D(hdir, "noncapture_lost", "count of lost non-capture events", 1, -0.5, 0.5);

    // True distribution of all events used for the unfolding
    mcproton_ptot_ = hf.DefineTH1D(hdir, "mcproton_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcdeuteron_ptot_ = hf.DefineTH1D(hdir, "mcdeuteron_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    noncapture_count_ = hf.DefineTH1D(hdir, "noncapture_count", "noncapture count, input", 1, -0.5, 0.5);

    // Migration matrices for the contained channel, two generator binnings
    containedMigration1_ = hf.DefineTH3D(hdir, "containedMigration1",
                                         "Contained channel migration, generator binning 1",
                                         gen1nbins, gen1pmin, gen1pmax,
                                         recoContPNbins, recoContPMin, recoContPMax,
                                         recoContRNbins, recoContRMin, recoContRMax);

    containedMigration1_->GetXaxis()->SetTitle("p true, MeV/c");
    containedMigration1_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedMigration1_->GetZaxis()->SetTitle("range");

    // Migration matrices for the non contained channel, two generator binnings
    uncontainedMigration1_ = hf.DefineTH2D(hdir, "uncontainedMigration1",
                                           "Non-contained channel migration, generator binning 1",
                                           gen1nbins, gen1pmin, gen1pmax,
                                           recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedMigration1_->SetOption("colz");
    uncontainedMigration1_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedMigration1_->GetYaxis()->SetTitle("p reco, MeV/c");

    hTruthTrkContained_.init(hf, hdir+"/MCTruthTrkContained", conf);
    hTruthTrkUncontained_.init(hf, hdir+"/MCTruthTrkUncontained", conf);

    hResolutionContained_.init(hf, hdir+"/ResolutionContained", conf);
    hResolutionUncontained_.init(hf, hdir+"/ResolutionUncontained", conf);
  }

}

//================================================================
void HistMuCapAnalysisChannels::fill(const EventClass& evt,
                                     int iPosTrack,
                                     int iNegTrack,
                                     bool isPosTrackContained,
                                     double rangePIDVar)
{
  //----------------------------------------------------------------
  // Figure out an exclusive analysis channel for this event

  bool eventUsedInAChannel = false;

  if(iNegTrack == -1) { // Veto DIO events
    if(iPosTrack != -1) { // Got a reconstructed capture track

      const double prec = evt.ptot[iPosTrack];

      if(isPosTrackContained) {
        // The "contained tracks" analysis channel
        eventUsedInAChannel = true;
        contained_prange_->Fill(prec, rangePIDVar);
        if(doMCTruth_) {
          hTruthTrkContained_.fill(evt);
          hResolutionContained_.fill(evt, iPosTrack);

          const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
          if(imcvtxStart  != -1) {
            containedMigration1_->Fill(evt.mcvertex_ptot[imcvtxStart], prec, rangePIDVar);
          }
        }
      }
      else { // The non-contained tracks channel
        eventUsedInAChannel = true;
        uncontained_p_->Fill(prec);
        if(doMCTruth_) {
          hTruthTrkUncontained_.fill(evt);
          hResolutionUncontained_.fill(evt, iPosTrack);

          const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
          if(imcvtxStart  != -1) {
            uncontainedMigration1_->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
          }
        }
      }
    }
    else {
      // No good positive tracks. Can do a hit-based analysis here.
    }
  }
  if(doMCTruth_ && !eventUsedInAChannel) {
    const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
    if(imcvtxStart  != -1) {
      lost1_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
    }
    else {
      noncapture_lost_->Fill(0.);
    }
  }

  // Truth momentum with the binning used in the unfolding
  // Keep protons and deuterons separately to compare hadd-ed
  // pseudodata truth to unfolding results.
  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  if(imcvtxStart  != -1) {
    const int mcParticle = evt.mctrack_pid[evt.iCaptureMcTrk];
    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      mcproton_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      mcdeuteron_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
      break;
    }
  }
  else {
    noncapture_count_->Fill(0.);
  }
}

//================================================================
