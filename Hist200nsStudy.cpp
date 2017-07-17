// Andrei Gaponenko, 2013

#include "Hist200nsStudy.h"

#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"
#include "TDCHitWP.h"
#include "WireCluster.h"

//================================================================
void Hist200nsStudy::init(const std::string& hdir,
                          HistogramFactory &hf,
                          const DetectorGeo& geom,
                          const ConfigFile& cconf)
{
  // settings for the rcp1 channel, copied from DefaultConfig.cpp
  ConfigFile& conf = const_cast<ConfigFile&>(cconf);
  conf.add(hdir+"/ccut/cutMaxPlane", 51);
  conf.add(hdir+"/ccut/cutMaxRout", 15.);
  conf.add(hdir+"/ccut/useExtendedRange", true);

  cvp_.init(hdir+"/ccut", hf, geom, conf);


  hitTrackdtDC_ =  hf.DefineTH1D(hdir, "hitTrackdtDC", "t(hit)-t(track)", 240, -100., 1100.);
  hitTrackdtPC_ =  hf.DefineTH1D(hdir, "hitTrackdtPC", "t(hit)-t(track)", 240, -100., 1100.);

  //----------------------------------------------------------------
  // Proton/Deuteron PID from range vs momentum.

  rangemom_contained_ = hf.DefineTH2D(hdir, "rangemom_contained", "track range vs momentum, contained", 250, 0., 250., 12, 8., 32.);;
  rangemom_contained_->SetOption("colz");
  rangemom_contained_->GetXaxis()->SetTitle(cvp_.xtitle().c_str());
  rangemom_contained_->GetYaxis()->SetTitle(cvp_.ytitle().c_str());


  rangemom_proton_ = hf.DefineTH2D(hdir, "rangemom_proton", "track range vs momentum, contained", 250, 0., 250., 12, 8., 32.);;
  rangemom_proton_->SetOption("colz");
  rangemom_proton_->GetXaxis()->SetTitle(cvp_.xtitle().c_str());
  rangemom_proton_->GetYaxis()->SetTitle(cvp_.ytitle().c_str());

  rangemom_deuteron_ = hf.DefineTH2D(hdir, "rangemom_deuteron", "track range vs momentum, contained", 250, 0., 250., 12, 8., 32.);;
  rangemom_deuteron_->SetOption("colz");
  rangemom_deuteron_->GetXaxis()->SetTitle(cvp_.xtitle().c_str());
  rangemom_deuteron_->GetYaxis()->SetTitle(cvp_.ytitle().c_str());

  //----------------
  hits23MaxWidth_contained_ = hf.DefineTH1D(hdir, "dc23maxhitwidth_contained", "max hit width in DC23, contained", 150, 0., 1500.);
  hits23MaxWidth_proton_ = hf.DefineTH1D(hdir, "dc23maxhitwidth_proton", "max hit width in DC23, proton", 150, 0., 1500.);
  hits23MaxWidth_deuteron_ = hf.DefineTH1D(hdir, "dc23maxhitwidth_deuteron", "max hit width in DC23, deuteron", 150, 0., 1500.);

  //----------------------------------------------------------------
  hitsMaxWidth24vs23_all_ = hf.DefineTH2D(hdir, "maxwidth24vs23_all", "max hit width in DC24 vs 23, all", 150, 0., 1500., 150, 0., 1500.);
  hitsMaxWidth24vs23_all_->SetOption("colz");
  hitsMaxWidth24vs23_all_->GetXaxis()->SetTitle("DC23 max width");
  hitsMaxWidth24vs23_all_->GetYaxis()->SetTitle("DC24 max width");

  hitsMaxWidth24vs23_contained_ = hf.DefineTH2D(hdir, "maxwidth24vs23_contained", "max hit width in DC24 vs 23, contained", 150, 0., 1500., 150, 0., 1500.);
  hitsMaxWidth24vs23_contained_->SetOption("colz");
  hitsMaxWidth24vs23_contained_->GetXaxis()->SetTitle("DC23 max width");
  hitsMaxWidth24vs23_contained_->GetYaxis()->SetTitle("DC24 max width");

  hitsMaxWidthPC8vsDC23_all_ = hf.DefineTH2D(hdir, "maxwidthPC8vsDC23_all", "max hit width in PC8 vs DC23, all", 150, 0., 1500., 150, 0., 1500.);
  hitsMaxWidthPC8vsDC23_all_->SetOption("colz");
  hitsMaxWidthPC8vsDC23_all_->GetXaxis()->SetTitle("DC23 max width");
  hitsMaxWidthPC8vsDC23_all_->GetYaxis()->SetTitle("PC8 max width");

  hitsMaxWidthPC8vsDC23_contained_ = hf.DefineTH2D(hdir, "maxwidthPC8vsDC23_contained", "max hit width in PC8 vs DC23, contained", 150, 0., 1500., 150, 0., 1500.);
  hitsMaxWidthPC8vsDC23_contained_->SetOption("colz");
  hitsMaxWidthPC8vsDC23_contained_->GetXaxis()->SetTitle("DC23 max width");
  hitsMaxWidthPC8vsDC23_contained_->GetYaxis()->SetTitle("PC8 max width");

  //----------------------------------------------------------------
}

//================================================================
void Hist200nsStudy::fill(const EventClass& evt, int iTrack, const ClustersByPlane& protonGlobalClusters) {
  // use same values as in windowing
  // NB: here these are w.r.t. the track instead of w.r.t. PC time
  const double winDCStart = -100.;
  const double winDCEnd = 1050.;

  // ???
  const double winPCStart = -100.;
  const double winPCEnd = +100.;

  const double trackTime = evt.hefit_time[iTrack];

  // Parameters for proton vs deuteron PID from Anthony Hillairet's
  // posting on 2016-06-01
  // https://twist.triumf.ca/~e614/forum/view.php?site=twist&bn=twist_physics&key=1464800416
  //
  // "diagonal cut: y = 0.4*x - 22"

  static const double cutPID_slope = 0.4;
  static const double cutPID_offset = -22.;
  // an additional cut to get higher purity samples
  static const double cutPID_ymin = 12.;

  std::vector<double> maxWidthDC(45, 0.);
  for(unsigned i=0; i< evt.dc_hits().size(); ++i) {
    const double dt = evt.dc_hits()[i].time() - trackTime;
    hitTrackdtDC_->Fill(dt);
    if((winDCStart < dt) && (dt < winDCEnd)) {
      const int plane = evt.dc_hits()[i].plane();
      double cw = evt.dc_hits()[i].width();
      if(maxWidthDC[plane] < cw) {
        maxWidthDC[plane] = cw;
      }
    }
  }

  std::vector<double> maxWidthPC(13, 0.);
  for(unsigned i=0; i< evt.pc_hits().size(); ++i) {
    const double dt = evt.pc_hits()[i].time() - trackTime;
    hitTrackdtPC_->Fill(dt);
    if((winPCStart < dt) && (dt < winPCEnd)) {
      const int plane = evt.pc_hits()[i].plane();
      double cw = evt.pc_hits()[i].width();
      if(maxWidthPC[plane] < cw) {
        maxWidthPC[plane] = cw;
      }
    }
  }

  hitsMaxWidth24vs23_all_->Fill(maxWidthDC[23], maxWidthDC[24]);
  hitsMaxWidthPC8vsDC23_all_->Fill(maxWidthDC[23], maxWidthPC[8]);

  MuCapContainedVars::Result cv = cvp_.compute(evt,iTrack,protonGlobalClusters);
  if(cv.contained) {

    hitsMaxWidth24vs23_contained_->Fill(maxWidthDC[23], maxWidthDC[24]);
    hitsMaxWidthPC8vsDC23_contained_->Fill(maxWidthDC[23], maxWidthPC[8]);

    hits23MaxWidth_contained_->Fill(maxWidthDC[23]);
    rangemom_contained_->Fill(cv.xvar, cv.yvar);

    if(cutPID_ymin < cv.yvar) { // Range is not too smal - fill the PID-separated histos
      if(cutPID_slope * cv.xvar + cutPID_offset < cv.yvar) {
        // proton
        rangemom_proton_->Fill(cv.xvar, cv.yvar);
        hits23MaxWidth_proton_->Fill(maxWidthDC[23]);
      }
      else {
        // deuteron
        rangemom_deuteron_->Fill(cv.xvar, cv.yvar);
        hits23MaxWidth_deuteron_->Fill(maxWidthDC[23]);
      }
    }
  }

}

//================================================================
