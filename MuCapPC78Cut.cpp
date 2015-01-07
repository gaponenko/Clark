// Andrei Gaponenko, 2014

#include "MuCapPC78Cut.h"

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
void MuCapPC78Cut::init(HistogramFactory& hf,
                        const std::string& hdir,
                        const DetectorGeo& geom,
                        const ConfigFile& conf)
{
  cutPC7width_ = conf.read<double>("MuCapture/CutPC78/cutTDCwidth");
  cutPChalfsumWidth_ = cutPC7width_;

  pc7width_before_ = hf.DefineTH1D(hdir, "pc7widthBeforeCut", "PC7 TDC width, before cut", 1000, 0., 1000.);
  pc7width_after_ = hf.DefineTH1D(hdir, "pc7widthAfterCut", "PC7 TDC width, after cut", 1000, 0., 1000.);

  pchalfsumwidth_before_ = hf.DefineTH1D(hdir, "pchalfsumwidthBeforeCut", "(PC7+PC8 TDC width)/2, before cut", 1000, 0., 1000.);
  pchalfsumwidth_after_ = hf.DefineTH1D(hdir, "pchalfsumwidthAfterCut", "(PC7+PC8 TDC width)/2, after cut", 1000, 0., 1000.);

  pc78width_before_ = hf.DefineTH2D(hdir, "pc78widthBeforeCut", "TDC width 8 vs 7, widest hits, before cut", 50, 0., 1000., 50, 0., 1000.);
  pc78width_before_->SetOption("colz");
  pc78width_before_->GetXaxis()->SetTitle("PC7 hit width");
  pc78width_before_->GetYaxis()->SetTitle("PC8 hit width");

  pc78width_after_ = hf.DefineTH2D(hdir, "pc78widthAfterCut", "TDC width 8 vs 7, widest hits, after cut", 50, 0., 1000., 50, 0., 1000.);
  pc78width_after_->SetOption("colz");
  pc78width_after_->GetXaxis()->SetTitle("PC7 hit width");
  pc78width_after_->GetYaxis()->SetTitle("PC8 hit width");
}

//================================================================
bool MuCapPC78Cut::accepted(const TDCHitWPPtr& pc7WidestHit) {
  pc7width_before_->Fill(pc7WidestHit->width());
  const bool res = (cutPC7width_ <= pc7WidestHit->width());
  if(res) {
    pc7width_after_->Fill(pc7WidestHit->width());
  }
  return res;
}

//================================================================
bool MuCapPC78Cut::accepted(const TDCHitWPPtr& hit7, const TDCHitWPPtr& hit8) {
  const double w7 = hit7->width();
  const double w8 = hit8->width();
  const double whalfsum = (w7+w8)/2;
  pc78width_before_->Fill(w7, w8);
  pchalfsumwidth_before_->Fill(whalfsum);
  const bool res = (cutPChalfsumWidth_ <= whalfsum);
  if(res) {
    pc78width_after_->Fill(w7, w8);
    pchalfsumwidth_after_->Fill(whalfsum);
  }
  return res;
}

//================================================================
bool MuCapPC78Cut::accepted(const EventClass& evt, const ClustersByPlane& protonGlobalClusters) {

  // require a PC7 hit
  if(protonGlobalClusters[29].empty()) {
      return false;
  }

  if(protonGlobalClusters[30].empty()) { // No PC8 info
    TDCHitWPPtr hit7 = maxTDCWidthHit(protonGlobalClusters.at(29));
    return accepted(hit7);
  }
  else {
    // Use both PC7 and PC8
    TDCHitWPPtr hit7 = maxTDCWidthHit(protonGlobalClusters.at(29));
    TDCHitWPPtr hit8 = maxTDCWidthHit(protonGlobalClusters.at(30));
    return accepted(hit7, hit8);
  }
}

//================================================================
