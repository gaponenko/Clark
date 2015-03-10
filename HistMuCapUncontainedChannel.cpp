// Andrei Gaponenko, 2014

#include "HistMuCapUncontainedChannel.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapUncontainedChannel::init(HistogramFactory& hf,
                                     const std::string& hgrandtopdir,
                                     const std::string& channelsetname,
                                     const DetectorGeo& geom,
                                     const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  const std::string hdir = hgrandtopdir+"/"+channelsetname+"/contained";

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());
  h_cuts_r->SetOption("hist text");

  h_cuts_p = hf.DefineTH1D(hdir, "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);
  h_cuts_p->SetOption("hist text");

  //----------------------------------------------------------------
  // "channel" analysis histograms
  const int recoUncPNbins = 88;
  const double recoUncPMin = 30.;
  const double recoUncPMax = 250.;

  //----------------------------------------------------------------
  reco_ = hf.DefineTH1D(hdir, "ptot", "ptot, non contained", recoUncPNbins, recoUncPMin, recoUncPMax);
  reco_->GetXaxis()->SetTitle("p, MeV/c");

  if(doMCTruth_) {
    reco_mcproton_ = hf.DefineTH1D(hdir, "ptot_mcproton", "ptot, non contained mcproton", recoUncPNbins, recoUncPMin, recoUncPMax);
    reco_mcproton_->GetXaxis()->SetTitle("p, MeV/c");

    reco_mcdeuteron_ = hf.DefineTH1D(hdir, "ptot_mcdeuteron", "ptot, non contained mcdeuteron", recoUncPNbins, recoUncPMin, recoUncPMax);
    reco_mcdeuteron_->GetXaxis()->SetTitle("p, MeV/c");

    reco_mcdio_ = hf.DefineTH1D(hdir, "ptot_mcdio", "ptot, non contained mcdio", recoUncPNbins, recoUncPMin, recoUncPMax);
    reco_mcdio_->GetXaxis()->SetTitle("p, MeV/c");
  }

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning must be consistent for all channels
    const int gen1nbins = conf.read<int>("MuCapture/"+channelsetname+"/numGeneratorBins");
    const double gen1pmin = conf.read<double>("MuCapture/"+channelsetname+"/genpmin");
    const double gen1pmax = conf.read<double>("MuCapture/"+channelsetname+"/genpmax");

    // Migration matrices for the non contained channel
    migration_ = hf.DefineTH2D(hdir,"migration",
                                          "Non-contained channel migration",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          recoUncPNbins, recoUncPMin, recoUncPMax);

    migration_->SetOption("colz");
    migration_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_->GetYaxis()->SetTitle("p reco, MeV/c");

    migration_mcproton_ = hf.DefineTH2D(hdir,"migration_mcproton",
                                                   "Non-contained channel migration, proton",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    migration_mcproton_->SetOption("colz");
    migration_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_mcproton_->GetYaxis()->SetTitle("p reco, MeV/c");

    migration_mcdeuteron_ = hf.DefineTH2D(hdir,"migration_mcdeuteron",
                                                   "Non-contained channel migration, deuteron",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    migration_mcdeuteron_->SetOption("colz");
    migration_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_mcdeuteron_->GetYaxis()->SetTitle("p reco, MeV/c");

    // Contamination matrices for the non contained channel
    contamination_ = hf.DefineTH2D(hdir,"contamination",
                                          "Non-contained channel contamination",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          recoUncPNbins, recoUncPMin, recoUncPMax);

    contamination_->SetOption("colz");
    contamination_->GetXaxis()->SetTitle("p true, MeV/c");
    contamination_->GetYaxis()->SetTitle("p reco, MeV/c");

    contamination_mcproton_ = hf.DefineTH2D(hdir,"contamination_mcproton",
                                                   "Non-contained channel contamination, proton",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    contamination_mcproton_->SetOption("colz");
    contamination_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    contamination_mcproton_->GetYaxis()->SetTitle("p reco, MeV/c");

    contamination_mcdeuteron_ = hf.DefineTH2D(hdir,"contamination_mcdeuteron",
                                                   "Non-contained channel contamination, deuteron",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    contamination_mcdeuteron_->SetOption("colz");
    contamination_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    contamination_mcdeuteron_->GetYaxis()->SetTitle("p reco, MeV/c");
  }

}

//================================================================
bool HistMuCapUncontainedChannel::accepted(const EventClass& evt,
                                         bool referenceSampleAccepted,
                                         int iPosTrack,
                                         int iNegTrack)

{
  CutNumber c = analyzeEvent(evt, referenceSampleAccepted, iPosTrack, iNegTrack);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c == CUTS_ACCEPTED;
}

//================================================================
HistMuCapUncontainedChannel::CutNumber
HistMuCapUncontainedChannel::analyzeEvent(const EventClass& evt,
                                        bool referenceSampleAccepted,
                                        int iPosTrack,
                                        int iNegTrack)
{
  if(iPosTrack == -1) {
    return CUT_POSTRK;
  }

  if(iNegTrack != -1) {
    return CUT_DIO;
  }

  const double prec = evt.ptot[iPosTrack];
  reco_->Fill(prec);

  if(doMCTruth_) {

    const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
    // Simulated DIO have no easily accessible MC truth.  We'll tread PID=zero as DIO down in this code.
    const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

    if(imcvtxStart  != -1) {
      (referenceSampleAccepted ? migration_ : contamination_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
    }

    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      reco_mcproton_->Fill(prec);
      (referenceSampleAccepted ? migration_mcproton_ : contamination_mcproton_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      reco_mcdeuteron_->Fill(prec);
      (referenceSampleAccepted ? migration_mcdeuteron_ : contamination_mcdeuteron_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
      break;
    case 0:
      reco_mcdio_->Fill(prec);
      break;
    }
  }

  return CUTS_ACCEPTED;
}

//================================================================
