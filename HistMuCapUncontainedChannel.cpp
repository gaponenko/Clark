// Andrei Gaponenko, 2014

#include "HistMuCapUncontainedChannel.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

#define SETUP_RECOMC_HIST_UNCONTAINED(hname)                            \
  do {                                                                  \
    reco##hname##_ = hf.DefineTH1D(hdir, "ptot"#hname,                  \
                                   "ptot, non contained"#hname,         \
                                   recopnbins, recopmin, recopmax);     \
                                                                        \
    reco##hname##_->GetXaxis()->SetTitle("p, MeV/c");                   \
  } while(0)

#define SETUP_RESPONSE_HIST_UNCONTAINTED(kind, hname)                   \
  do {                                                                  \
    kind##hname##_ = hf.DefineTH2D(hdir, #kind #hname,                  \
                                   "Non-contained channel "#kind#hname, \
                                   gen1nbins, gen1pmin, gen1pmax,       \
                                   recopnbins, recopmin, recopmax);     \
                                                                        \
                                                                        \
    kind##hname##_->SetOption("colz");                                  \
    kind##hname##_->GetXaxis()->SetTitle("p true, MeV/c");              \
    kind##hname##_->GetYaxis()->SetTitle("p reco, MeV/c");              \
  } while(0)




//================================================================
void HistMuCapUncontainedChannel::init(HistogramFactory& hf,
                                       const std::string& hgrandtopdir,
                                       const std::string& channelsetname,
                                       const DetectorGeo& geom,
                                       const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  const std::string htopdir = hgrandtopdir+"/"+channelsetname;
  const std::string hdir = htopdir+"/uncontained";

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
  const int recopnbins = conf.read<int>(hdir+"/recopnbins");
  const double recopmin = conf.read<double>(hdir+"/recopmin");
  const double recopmax = conf.read<double>(hdir+"/recopmax");

  //----------------------------------------------------------------
  SETUP_RECOMC_HIST_UNCONTAINED();

  if(doMCTruth_) {
    SETUP_RECOMC_HIST_UNCONTAINED(_mcproton);
    SETUP_RECOMC_HIST_UNCONTAINED(_mcdeuteron);
    SETUP_RECOMC_HIST_UNCONTAINED(_mctriton);
    SETUP_RECOMC_HIST_UNCONTAINED(_mcalpha);
    SETUP_RECOMC_HIST_UNCONTAINED(_mcdio);
  }

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning must be consistent for all channels
    const int gen1nbins = conf.read<int>(htopdir+"/numGeneratorBins");
    const double gen1pmin = conf.read<double>(htopdir+"/genpmin");
    const double gen1pmax = conf.read<double>(htopdir+"/genpmax");

    // Migration matrices for the non contained channel
    SETUP_RESPONSE_HIST_UNCONTAINTED(migration,);
    SETUP_RESPONSE_HIST_UNCONTAINTED(migration, _mcproton);
    SETUP_RESPONSE_HIST_UNCONTAINTED(migration, _mcdeuteron);
    SETUP_RESPONSE_HIST_UNCONTAINTED(migration, _mctriton);
    SETUP_RESPONSE_HIST_UNCONTAINTED(migration, _mcalpha);

    // Contamination matrices for the non contained channel
    SETUP_RESPONSE_HIST_UNCONTAINTED(contamination,);
    SETUP_RESPONSE_HIST_UNCONTAINTED(contamination, _mcproton);
    SETUP_RESPONSE_HIST_UNCONTAINTED(contamination, _mcdeuteron);
    SETUP_RESPONSE_HIST_UNCONTAINTED(contamination, _mctriton);
    SETUP_RESPONSE_HIST_UNCONTAINTED(contamination, _mcalpha);
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
    case MuCapUtilities::PID_G3_TRITON:
      reco_mctriton_->Fill(prec);
      (referenceSampleAccepted ? migration_mctriton_ : contamination_mctriton_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
      break;
    case MuCapUtilities::PID_G3_ALPHA:
      reco_mcalpha_->Fill(prec);
      (referenceSampleAccepted ? migration_mcalpha_ : contamination_mcalpha_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
      break;
    case 0:
      reco_mcdio_->Fill(prec);
      break;
    default:
      std::ostringstream os;
      os<<"HistMuCapUncontainedChannel: unknown capture PID="<<mcParticle;
      throw std::runtime_error(os.str());
    }
  }

  return CUTS_ACCEPTED;
}

//================================================================
