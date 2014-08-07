// Andrei Gaponenko, 2014

#include "MuCapTrkContainmentCut.h"

#include <cmath>

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "EventClass.h"
#include "MuCapUtilities.h"

//================================================================
void MuCapTrkContainmentCut::init(const std::string& hdir,
                                  HistogramFactory &hf,
                                  const DetectorGeo&,
                                  const ConfigFile& conf)
{
  cutMaxPlane_ = conf.read<int>("MuCapture/TrkContainmentCut/cutMaxPlane");
  cutMaxRout_ = conf.read<double>("MuCapture/TrkContainmentCut/cutMaxRout");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "cuts_r", "Tracks rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());
  h_cuts_r->SetOption("hist text");

  h_cuts_p = hf.DefineTH1D(hdir, "cuts_p", "Tracks before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);
  h_cuts_p->SetOption("hist text");

  //----------------------------------------------------------------
  hExtendedLastPlane_ = hf.DefineTH1D(hdir, "extendedLastPlane", "Extended last plane", 56, 0.5, 56.5);
  hRout_ = hf.DefineTH1D(hdir, "rout", "Rout", 400, 0., 40.);
}

//================================================================
bool MuCapTrkContainmentCut::contained(const EventClass& evt, int itrack,  const ClustersByPlane& protonGlobalClusters)
{
  CutNumber c = analyzeTrack(evt, itrack, protonGlobalClusters);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c==CUTS_ACCEPTED;
}

//================================================================
MuCapTrkContainmentCut::CutNumber MuCapTrkContainmentCut::
analyzeTrack(const EventClass& evt, int i, const ClustersByPlane& protonGlobalClusters)
{

  int extendedLastPlane = MuCapUtilities::findExtendedLastPlane(evt, i, protonGlobalClusters);
  hExtendedLastPlane_->Fill(extendedLastPlane);
  if(extendedLastPlane > cutMaxPlane_) {
    return CUT_ZRANGE;
  }

  //----------------------------------------------------------------
  const double rout =
    sqrt(std::pow(evt.hefit_ucenter[i], 2) + std::pow(evt.hefit_vcenter[i], 2))
    + evt.radius[i];

  hRout_->Fill(rout);
  if(rout > cutMaxRout_) {
    return CUT_ROUT;
  }

  //----------------------------------------------------------------
  return CUTS_ACCEPTED;
}
