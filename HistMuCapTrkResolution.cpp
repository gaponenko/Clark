// Andrei Gaponenko, 2013

#include "HistMuCapTrkResolution.h"

#include "TH1.h"
#include "TProfile.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapTrkResolution::init(HistogramFactory& hf,
                                  const std::string& hdir,
                                  const ConfigFile& conf)
{

  hGlobalResolution_ = hf.DefineTH1D(hdir, "globalRes", "Momentum resolution in the fiducial region", 300, -100., 50.);
  hGlobalResolution_->GetXaxis()->SetTitle("pmc [MeV/c]");

  hMomResVsMom_ = hf.DefineTProfile(hdir, "momResVsMom", "Momentum resolution vs MC momentum", 450, 0., 450., "s");
  hMomResVsMom_->GetXaxis()->SetTitle("pmc [MeV/c]");
  hMomResVsMom_->GetYaxis()->SetTitle("prec-pmc [MeV/c]");

  hMomResVsCosth_ = hf.DefineTProfile(hdir, "momResVsCosth", "Momentum resolution vs MC cos(theta)", 100, -1., 1., "s");
  hMomResVsCosth_->GetXaxis()->SetTitle("cos(theta) MC");
  hMomResVsCosth_->GetYaxis()->SetTitle("prec-pmc [MeV/c]");

  hMomBias2DMC_ = hf.DefineTProfile2D(hdir, "momBias2DMC", "Momentum bias in MC (p,costh)", 450, 0., 450., 100, -1., 1., "s");
  hMomBias2DMC_->GetXaxis()->SetTitle("pmc [MeV/c]");
  hMomBias2DMC_->GetYaxis()->SetTitle("cos(theta) MC");
  hMomBias2DMC_->SetOption("colz");

  hMomBias2DReco_ = hf.DefineTProfile2D(hdir, "momBias2DReco", "Momentum bias in reco (p,costh)", 450, 0., 450., 100, -1., 1., "s");
  hMomBias2DReco_->GetXaxis()->SetTitle("prec [MeV/c]");
  hMomBias2DReco_->GetYaxis()->SetTitle("cos(theta) rec");
  hMomBias2DReco_->SetOption("colz");
}

//================================================================
void HistMuCapTrkResolution::fill(const EventClass& evt, int iRecoTrack)
{
  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  if(imcvtxStart != -1) {
    const double ptrue = evt.mcvertex_ptot[imcvtxStart];
    const double costrue = evt.mcvertex_costh[imcvtxStart];

    const double prec = evt.ptot[iRecoTrack];
    const double cosrec = evt.costh[iRecoTrack];

    const double dp = prec - ptrue;

    hGlobalResolution_->Fill(dp);
    hMomResVsMom_->Fill(ptrue, dp);
    hMomResVsCosth_->Fill(costrue, dp);

    hMomBias2DMC_->Fill(ptrue, costrue, dp);
    hMomBias2DReco_->Fill(prec, cosrec, dp);
  }
}

//================================================================
