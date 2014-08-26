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

  hMomResVsMom_ = hf.DefineTProfile(hdir, "momResVsMom", "Momentum resolution vs MC momentum", 300, 0., 300., "s");

  hMomResVsCosth_ = hf.DefineTProfile(hdir, "momResVsCosth", "Momentum resolution vs MC cos(theta)", 100, -1., 1., "s");
}

//================================================================
void HistMuCapTrkResolution::fill(const EventClass& evt, int iRecoTrack)
{
  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  if(imcvtxStart != -1) {
    const double ptrue = evt.mcvertex_ptot[imcvtxStart];
    const double costrue = evt.mcvertex_costh[imcvtxStart];
    const double prec = evt.ptot[iRecoTrack];

    hGlobalResolution_->Fill(prec - ptrue);
    hMomResVsMom_->Fill(ptrue, prec-ptrue);
    hMomResVsCosth_->Fill(costrue, prec-ptrue);
  }
}

//================================================================
