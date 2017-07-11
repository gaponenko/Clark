// Andrei Gaponenko, 2013

#include "HistXT4.h"

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
void HistXT4WireGap::init(const std::string& hdir,
                             unsigned cutWireGap,
                             const std::string& suffix,
                             HistogramFactory& hf,
                             const ConfigFile& conf)
{
  cutWireGap_ = cutWireGap;

  //----------------
  outwide_dt_all_ =  hf.DefineTH1D(hdir, "outwide_dt_all", "t(out)- t(wide), all hits"+suffix, 400, -1000., 1000.);

  outwide_ww_all_ =  hf.DefineTH2D(hdir, "outwide_ww_all", "w(out) vs w(wide), all hits"+suffix, 100, 0., 1000., 100, 0., 1000.);
  outwide_ww_all_->SetOption("colz");

  outwide_wt_all_ =  hf.DefineTH2D(hdir, "outwide_wt_all", "w(out) vs dt, all hits"+suffix, 400, -1000., 1000., 100, 0., 1000.);
  outwide_wt_all_->SetOption("colz");

  //----------------
  outwide_dt_xt_ =  hf.DefineTH1D(hdir, "outwide_dt_xt", "t(out)- t(wide), xt hits"+suffix, 400, -1000., 1000.);

  outwide_ww_xt_ =  hf.DefineTH2D(hdir, "outwide_ww_xt", "w(out) vs w(wide), xt hits"+suffix, 100, 0., 1000., 100, 0., 1000.);
  outwide_ww_xt_->SetOption("colz");

  outwide_wt_xt_ =  hf.DefineTH2D(hdir, "outwide_wt_xt", "w(out) vs dt, xt hits"+suffix, 400, -1000., 1000., 100, 0., 1000.);
  outwide_wt_xt_->SetOption("colz");

  //----------------
  outwide_dt_nxt_ =  hf.DefineTH1D(hdir, "outwide_dt_nxt", "t(out)- t(wide), nxt hits"+suffix, 400, -1000., 1000.);

  outwide_ww_nxt_ =  hf.DefineTH2D(hdir, "outwide_ww_nxt", "w(out) vs w(wide), nxt hits"+suffix, 100, 0., 1000., 100, 0., 1000.);
  outwide_ww_nxt_->SetOption("colz");

  outwide_wt_nxt_ =  hf.DefineTH2D(hdir, "outwide_wt_nxt", "w(out) vs dt, nxt hits"+suffix, 400, -1000., 1000., 100, 0., 1000.);
  outwide_wt_nxt_->SetOption("colz");

}

//================================================================
void HistXT4WireGap::fill(const TDCHitWPPtrCollection& planeHits) {
  if(!planeHits.empty()) {
    // Find the widest hits
    double iwide = 0;
    double maxw = planeHits[iwide]->width();
    for(unsigned i=1; i<planeHits.size(); ++i) {
      if(maxw < planeHits[i]->width()) {
        iwide = i;
        maxw = planeHits[iwide]->width();
      }
    }

    const TDCHitWPPtr& wideHit = planeHits[iwide];

    // Fill distributions for hits that are "far enough" from the widest hit
    for(unsigned i=0; i<planeHits.size(); ++i) {
      if(std::abs(wideHit->cell() - planeHits[i]->cell()) > cutWireGap_) {

        const double dt = planeHits[i]->time() - wideHit->time();

        outwide_dt_all_->Fill(dt);
        outwide_ww_all_->Fill(wideHit->width(), planeHits[i]->width());
        outwide_wt_all_->Fill(dt, planeHits[i]->width());

        if(planeHits[i]->xtalk()) {
          outwide_dt_xt_->Fill(dt);
          outwide_ww_xt_->Fill(wideHit->width(), planeHits[i]->width());
          outwide_wt_xt_->Fill(dt, planeHits[i]->width());
        }
        else {
          outwide_dt_nxt_->Fill(dt);
          outwide_ww_nxt_->Fill(wideHit->width(), planeHits[i]->width());
          outwide_wt_nxt_->Fill(dt, planeHits[i]->width());
        }
      }
    }
  }
}

//================================================================
void HistXT4::init(const std::string& hdir,
                   HistogramFactory &hf,
                   const ConfigFile& conf)
{
  trackCosth_ =  hf.DefineTH1D(hdir, "trackCosth", "track cos(theta)", 100, 0., 1.);
  trackAngle_ =  hf.DefineTH1D(hdir, "trackAngle", "track theta (degrees)", 100, 0., 90.);

  wg1_.init(hdir+"/wgap1", 1, ", wg1", hf, conf);
  wg2_.init(hdir+"/wgap2", 2, ", wg2", hf, conf);
}

//================================================================
void HistXT4::fill(const EventClass& evt, int iTrack) {
  // use same values as in windowing
  // NB: here these are w.r.t. the track instead of w.r.t. PC time
  const double winDCStart = -100.;
  const double winDCEnd = 1050.;

  const double trackTime = evt.hefit_time[iTrack];

  TDCHitWPPtrCollection dc23hitsAll;

  for(unsigned i=0; i< evt.dc_hits().size(); ++i) {
    if(evt.dc_hits()[i].plane() == 23) {
      const double dt = evt.dc_hits()[i].time() - trackTime;
      if((winDCStart < dt) && (dt < winDCEnd)) {
        dc23hitsAll.push_back(TDCHitWPPtr(evt.dc_hits(), i));
      }
    }
  }

  const double costh= evt.costh[iTrack];
  const double angle = atan2(evt.sinth[iTrack], evt.costh[iTrack]) * 180/(4.*atan(1.));

  trackCosth_->Fill(costh);
  trackAngle_->Fill(angle);

  if(angle < 40.) {
    wg1_.fill(dc23hitsAll);
  }
  else {
    wg2_.fill(dc23hitsAll);
  }

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //test: if(dc23hitsAll.empty()) {
  //test:   std::cout<<"Got an event with emtpy dc23hitsAll collection.  Dumping all DC hits:\n"
  //test:            <<evt.dc_hits()
  //test:            <<std::endl;
  //test:
  //test:   std::cout<<"Track time = "<<trackTime<<std::endl;
  //test: }
}

//================================================================
