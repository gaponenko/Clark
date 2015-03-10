// Andrei Gaponenko, 2014

#include "HistMuCapContainedChannel.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

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
  dnPosTrkContainment_.init(hdir+"/dnPosTrkContainment", hf, geom, conf);

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
  reco_ = hf.DefineTH2D(hdir, "reco", "reco",
                        xvarnbins, xvarmin, xvarmax,
                        yvarnbins, yvarmin, yvarmax);

  reco_->SetOption("colz");
  reco_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());
  reco_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());

  if(doMCTruth_) {
    reco_mcproton_ = hf.DefineTH2D(hdir, "reco_mcproton", "reco, mcproton",
                                   xvarnbins, xvarmin, xvarmax, yvarnbins, yvarmin, yvarmax);
    reco_mcproton_->SetOption("colz");
    reco_mcproton_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());
    reco_mcproton_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());

    reco_mcdeuteron_ = hf.DefineTH2D(hdir, "reco_mcdeuteron", "reco, mcdeuteron",
                                     xvarnbins, xvarmin, xvarmax, yvarnbins, yvarmin, yvarmax);
    reco_mcdeuteron_->SetOption("colz");
    reco_mcdeuteron_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());
    reco_mcdeuteron_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());

    reco_mcdio_ = hf.DefineTH2D(hdir, "reco_mcdio", "reco, mcdio",
                                xvarnbins, xvarmin, xvarmax, yvarnbins, yvarmin, yvarmax);
    reco_mcdio_->SetOption("colz");
    reco_mcdio_->GetXaxis()->SetTitle(cvp_->xtitle().c_str());
    reco_mcdio_->GetYaxis()->SetTitle(cvp_->ytitle().c_str());
  }

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning must be consistent for all channels
    const int gen1nbins = conf.read<int>(htopdir+"/numGeneratorBins");
    const double gen1pmin = conf.read<double>(htopdir+"/genpmin");
    const double gen1pmax = conf.read<double>(htopdir+"/genpmax");

    // Migration matrices for the contained channel
    migration_ = hf.DefineTH3D(hdir, "migration",
                               "Contained channel migration",
                               gen1nbins, gen1pmin, gen1pmax,
                               xvarnbins, xvarmin, xvarmax,
                               yvarnbins, yvarmin, yvarmax);

    migration_->GetXaxis()->SetTitle("p true [MeV/c]");
    migration_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());
    migration_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());

    migration_mcproton_ = hf.DefineTH3D(hdir, "migration_mcproton",
                                        "Contained channel migration, proton",
                                        gen1nbins, gen1pmin, gen1pmax,
                                        xvarnbins, xvarmin, xvarmax,
                                        yvarnbins, yvarmin, yvarmax);

    migration_mcproton_->GetXaxis()->SetTitle("p true [MeV/c]");
    migration_mcproton_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());
    migration_mcproton_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());

    migration_mcdeuteron_ = hf.DefineTH3D(hdir, "migration_mcdeuteron",
                                          "Contained channel migration, deuteron",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          xvarnbins, xvarmin, xvarmax,
                                          yvarnbins, yvarmin, yvarmax);

    migration_mcdeuteron_->GetXaxis()->SetTitle("p true [MeV/c]");
    migration_mcdeuteron_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());
    migration_mcdeuteron_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());

    // Contamination matrices for the contained channel
    contamination_ = hf.DefineTH3D(hdir, "contamination",
                                   "Contained channel contamination",
                                   gen1nbins, gen1pmin, gen1pmax,
                                   xvarnbins, xvarmin, xvarmax,
                                   yvarnbins, yvarmin, yvarmax);

    contamination_->GetXaxis()->SetTitle("p true [MeV/c]");
    contamination_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());
    contamination_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());

    contamination_mcproton_ = hf.DefineTH3D(hdir, "contamination_mcproton",
                                            "Contained channel contamination, proton",
                                            gen1nbins, gen1pmin, gen1pmax,
                                            xvarnbins, xvarmin, xvarmax,
                                            yvarnbins, yvarmin, yvarmax);

    contamination_mcproton_->GetXaxis()->SetTitle("p true [MeV/c]");
    contamination_mcproton_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());
    contamination_mcproton_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());

    contamination_mcdeuteron_ = hf.DefineTH3D(hdir, "contamination_mcdeuteron",
                                              "Contained channel contamination, deuteron",
                                              gen1nbins, gen1pmin, gen1pmax,
                                              xvarnbins, xvarmin, xvarmax,
                                              yvarnbins, yvarmin, yvarmax);

    contamination_mcdeuteron_->GetXaxis()->SetTitle("p true [MeV/c]");
    contamination_mcdeuteron_->GetYaxis()->SetTitle(cvp_->xtitle().c_str());
    contamination_mcdeuteron_->GetZaxis()->SetTitle(cvp_->ytitle().c_str());
  }

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

  if(!dnPosTrkContainment_.contained(evt,iPosTrack,protonGlobalClusters)) {
    return CUT_CONTAINED;
  }

  MuCapContainedVars::Result cv = cvp_->compute(evt,iPosTrack,protonGlobalClusters);
  reco_->Fill(cv.xvar, cv.yvar);

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
    case 0:
      reco_mcdio_->Fill(cv.xvar, cv.yvar);
      break;
    }
  }


  return CUTS_ACCEPTED;
}

//================================================================
