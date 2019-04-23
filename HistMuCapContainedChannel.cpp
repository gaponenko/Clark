// Andrei Gaponenko, 2014

#include "HistMuCapContainedChannel.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"


#define SETUP_RECOMC_HIST_CONTAINED(hname)                              \
  do {                                                                  \
    reco##hname##_ = hf.DefineTH2D(hdir, "reco"#hname, "reco"#hname,    \
                                   xvarnbins, xvarmin, xvarmax,         \
                                   yvarnbins, yvarmin, yvarmax);        \
    reco##hname##_->SetOption("colz");                                  \
    reco##hname##_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());       \
    reco##hname##_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());       \
  } while(0)

#define SETUP_RESPONSE_HIST_CONTAINTED(kind, hname)                     \
  do {                                                                  \
    kind##hname##_ = hf.DefineTH3D(hdir, #kind #hname,                  \
                                   "Contained channel "#kind#hname,     \
                                   gen1nbins, gen1pmin, gen1pmax,       \
                                   xvarnbins, xvarmin, xvarmax,         \
                                   yvarnbins, yvarmin, yvarmax);        \
                                                                        \
    kind##hname##_->GetXaxis()->SetTitle("p true [MeV/c]");             \
    kind##hname##_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());       \
    kind##hname##_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());       \
  } while(0)

//================================================================
void HistMuCapContainedChannel::init(HistogramFactory& hf,
                                     const std::string& hgrandtopdir,
                                     const std::string& channelsetname,
                                     const DetectorGeo& geom,
                                     const ConfigFile& conf,
                                     MuCapContainedVars::IVarProcessor& cvp)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  const std::string htopdir = hgrandtopdir+"/"+channelsetname;
  const std::string hdir = htopdir+"/contained";

  cvp_ = &cvp;
  cvp_->init(hdir+"/ccut", hf, geom, conf);

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
  const int xvarnbins = conf.read<int>(hdir+"/xvarnbins");
  const double xvarmin = conf.read<double>(hdir+"/xvarmin");
  const double xvarmax = conf.read<double>(hdir+"/xvarmax");

  const int yvarnbins = conf.read<int>(hdir+"/yvarnbins");
  const double yvarmin = conf.read<double>(hdir+"/yvarmin");
  const double yvarmax = conf.read<double>(hdir+"/yvarmax");

  //----------------------------------------------------------------
  SETUP_RECOMC_HIST_CONTAINED();
  if(doMCTruth_) {
    SETUP_RECOMC_HIST_CONTAINED(_mcproton);
    SETUP_RECOMC_HIST_CONTAINED(_mcdeuteron);
    SETUP_RECOMC_HIST_CONTAINED(_mctriton);
    SETUP_RECOMC_HIST_CONTAINED(_mcalpha);
    SETUP_RECOMC_HIST_CONTAINED(_mcdio);
  }

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning must be consistent for all channels
    const int gen1nbins = conf.read<int>(htopdir+"/numGeneratorBins");
    const double gen1pmin = conf.read<double>(htopdir+"/genpmin");
    const double gen1pmax = conf.read<double>(htopdir+"/genpmax");

    // Migration matrices for the contained channel
    SETUP_RESPONSE_HIST_CONTAINTED(migration,);
    SETUP_RESPONSE_HIST_CONTAINTED(migration, _mcproton);
    SETUP_RESPONSE_HIST_CONTAINTED(migration, _mcdeuteron);
    SETUP_RESPONSE_HIST_CONTAINTED(migration, _mctriton);
    SETUP_RESPONSE_HIST_CONTAINTED(migration, _mcalpha);

    // Contamination matrices for the contained channel
    SETUP_RESPONSE_HIST_CONTAINTED(contamination,);
    SETUP_RESPONSE_HIST_CONTAINTED(contamination, _mcproton);
    SETUP_RESPONSE_HIST_CONTAINTED(contamination, _mcdeuteron);
    SETUP_RESPONSE_HIST_CONTAINTED(contamination, _mctriton);
    SETUP_RESPONSE_HIST_CONTAINTED(contamination, _mcalpha);
  }

  //----------------------------------------------------------------
  if(conf.read<bool>("MuCapture/channels/"+channelsetname+"/doPCosthSlices", false)) {
    pcosrange_.resize(yvarnbins);
    for(unsigned i=0; i<yvarnbins; ++i) {
      std::ostringstream os;
      os<<"ptotVsCosth_rangebin"<<1+i;
      pcosrange_[i] = hf.DefineTH2D(hdir, os.str(), os.str(), 650, 0., 650., 100, -1., +1.);
      pcosrange_[i]->SetOption("colz");
      pcosrange_[i]->GetXaxis()->SetTitle("ptot, MeV/c");
      pcosrange_[i]->GetYaxis()->SetTitle("cos(theta)");
    }
  }
  //----------------------------------------------------------------
  if(conf.read<bool>("MuCapture/channels/"+channelsetname+"/doExtraHistos", false)) {
    // hardcoded binning for the one-off plot for the paper
    pidplot_ = hf.DefineTH2D(hdir, "pidplot", "pidplot", 300, 0., 300., 30, 0.5, 30.5);
    pidplot_->SetOption("colz");
    pidplot_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());
    pidplot_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());

    // fine binned ptot to extract peak position number for the paper
    pidplotFine_ = hf.DefineTH2D(hdir, "pidplotfine", "pidplotfine", 6000, 0., 300., 25, 5., 30.);
    pidplotFine_->SetOption("colz");
    pidplotFine_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());
    pidplotFine_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());

    pplaneFine_ = hf.DefineTH2D(hdir, "pplanefine", "pplanefine", 6000, 0., 300., 24, 32.5, 56.5);
    pplaneFine_->SetOption("colz");
    pplaneFine_->GetXaxis()->SetTitle("ptot");
    pplaneFine_->GetYaxis()->SetTitle("lastplane");
  }
  //----------------------------------------------------------------

}

//================================================================
bool HistMuCapContainedChannel::accepted(const EventClass& evt,
                                         bool referenceSampleAccepted,
                                         int iPosTrack,
                                         int iNegTrack,
                                         const ClustersByPlane& protonGlobalClusters)

{
  CutNumber c = analyzeEvent(evt, referenceSampleAccepted, iPosTrack, iNegTrack, protonGlobalClusters);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c == CUTS_ACCEPTED;
}

//================================================================
HistMuCapContainedChannel::CutNumber
HistMuCapContainedChannel::analyzeEvent(const EventClass& evt,
                                        bool referenceSampleAccepted,
                                        int iPosTrack,
                                        int iNegTrack,
                                        const ClustersByPlane& protonGlobalClusters)
{
  if(iPosTrack == -1) {
    return CUT_POSTRK;
  }

  if(iNegTrack != -1) {
    return CUT_DIO;
  }

  MuCapContainedVars::Result cv = cvp_->compute(evt,iPosTrack,protonGlobalClusters);
  if(!cv.contained) {
    return CUT_CONTAINED;
  }

  reco_->Fill(cv.xvar, cv.yvar);
  if(!pcosrange_.empty()) {
    int bin = reco_->GetYaxis()->FindFixBin(cv.yvar);
    if((0 < bin) && (bin <= pcosrange_.size() )) {
      pcosrange_[bin-1]->Fill(evt.ptot[iPosTrack], evt.costh[iPosTrack]);
    }
  }
  if(pidplot_) {
    pidplot_->Fill(cv.xvar, cv.yvar);
  }
  if(pidplotFine_) {
    pidplotFine_->Fill(cv.xvar, cv.yvar);
  }
  if(pplaneFine_) {
    int extendedLastPlane = MuCapUtilities::findExtendedLastPlane(evt, iPosTrack, protonGlobalClusters);
    pplaneFine_->Fill(evt.ptot[iPosTrack], extendedLastPlane);
  }

  if(doMCTruth_) {

    const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
    // Simulated DIO have no easily accessible MC truth.  We'll tread PID=zero as DIO down in this code.
    const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

    if(imcvtxStart  != -1) {
      (referenceSampleAccepted ? migration_ : contamination_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], cv.xvar, cv.yvar);
    }

    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      reco_mcproton_->Fill(cv.xvar, cv.yvar);
      (referenceSampleAccepted ? migration_mcproton_ : contamination_mcproton_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], cv.xvar, cv.yvar);
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      reco_mcdeuteron_->Fill(cv.xvar, cv.yvar);
      (referenceSampleAccepted ? migration_mcdeuteron_ : contamination_mcdeuteron_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], cv.xvar, cv.yvar);
      break;

    case MuCapUtilities::PID_G3_TRITON:
      reco_mctriton_->Fill(cv.xvar, cv.yvar);
      (referenceSampleAccepted ? migration_mctriton_ : contamination_mctriton_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], cv.xvar, cv.yvar);
      break;

    case MuCapUtilities::PID_G3_ALPHA:
      reco_mcalpha_->Fill(cv.xvar, cv.yvar);
      (referenceSampleAccepted ? migration_mcalpha_ : contamination_mcalpha_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], cv.xvar, cv.yvar);
      break;

    case 0:
      reco_mcdio_->Fill(cv.xvar, cv.yvar);
      break;
    default:
      std::ostringstream os;
      os<<"HistMuCapContainedChannel: unknown capture PID="<<mcParticle;
      throw std::runtime_error(os.str());
    }
  }


  return CUTS_ACCEPTED;
}

//================================================================
