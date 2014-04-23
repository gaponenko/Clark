// Andrei Gaponenko, 2013

#include "HistTrkQuality.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "Math/SpecFunc.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"
#include "EventClass.h"

//================================================================
void HistTrkQuality::init(HistogramFactory& hf,
                          const std::string& hdir,
                          const ConfigFile& conf)
{
  hchi2_ = hf.DefineTH1D(hdir, "chi2", "chi2", 500, 0., 250);

  hNDF_vs_chi2_ = hf.DefineTH2D(hdir, "ndfVsChi2", "NDF vs chi2", 500, 0., 250., 100, -0.5, 99.5);
  hNDF_vs_chi2_->SetOption("colz");;

  hprob_ = hf.DefineTH1D(hdir, "prob", "prob", 500, 0., 1.);

  hkineprob_ = hf.DefineTProfile2D(hdir, "kineprob", "prob profile2d", 300, 0., 300., 100, -1., +1.);
  hkineprob_->SetOption("colz");
}

//================================================================
void HistTrkQuality::fill(const EventClass& evt, int i)
{
  hchi2_->Fill(evt.hefit_chi2[i]);
  hNDF_vs_chi2_->Fill(evt.hefit_chi2[i], evt.hefit_ndof[i]);

  const double prob = 1. - ROOT::Math::inc_gamma(evt.hefit_ndof[i]/2., evt.hefit_chi2[i]);
  hprob_->Fill(prob);
  hkineprob_->Fill(evt.ptot[i], evt.costh[i], prob);
}

//================================================================
