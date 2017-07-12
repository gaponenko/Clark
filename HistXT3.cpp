// Andrei Gaponenko, 2013

#include "HistXT3.h"

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
void HistXT3TwoCluster::init(const std::string& hdir,
                             const std::string& suffix,
                             HistogramFactory& hf,
                             const ConfigFile& conf)
{
  twocluster_dt_ =  hf.DefineTH1D(hdir, "twocluster_dt", "cluster t(late)-t(early)  for 2-cluster events"+suffix, 200, 0., 1000.);

  twocluster_ww_ =  hf.DefineTH2D(hdir, "twocluster_ww", "cluster maxw(late) vs maxw(early)  for 2-cluster events"+suffix, 100, 0., 1000., 100, 0., 1000.);
  twocluster_ww_->SetOption("colz");

  twocluster_wt_ =  hf.DefineTH2D(hdir, "twocluster_wt", "cluster maxw(late) vs dt  for 2-cluster events"+suffix, 200, 0., 1000., 100, 0., 1000.);
  twocluster_wt_->SetOption("colz");

  twocluster_ss_ =  hf.DefineTH2D(hdir, "twocluster_ss", "cluster size late vs early for 2-cluster events"+suffix,
                                  10, 0.5, 10.5, 10, 0.5, 10.5);
  twocluster_ss_->SetOption("colz");

  //----------------
  hitcluster_dt_all_ =  hf.DefineTH1D(hdir, "hitcluster_dt_all", "hit t(late)- cluster t(early)  for 2-cluster events, all hits"+suffix, 200, 0., 1000.);

  hitcluster_ww_all_ =  hf.DefineTH2D(hdir, "hitcluster_ww_all", "hit w(late) vs cluster maxw(early)  for 2-cluster events, all hits"+suffix, 100, 0., 1000., 100, 0., 1000.);
  hitcluster_ww_all_->SetOption("colz");

  hitcluster_wt_all_ =  hf.DefineTH2D(hdir, "hitcluster_wt_all", "hit w(late) vs hit-cluster dt  for 2-cluster events, all hits"+suffix, 200, 0., 1000., 100, 0., 1000.);
  hitcluster_wt_all_->SetOption("colz");

  //----------------
  hitcluster_dt_xt_ =  hf.DefineTH1D(hdir, "hitcluster_dt_xt", "hit t(late)- cluster t(early)  for 2-cluster events, xt hits"+suffix, 200, 0., 1000.);

  hitcluster_ww_xt_ =  hf.DefineTH2D(hdir, "hitcluster_ww_xt", "hit w(late) vs cluster maxw(early)  for 2-cluster events, xt hits"+suffix, 100, 0., 1000., 100, 0., 1000.);
  hitcluster_ww_xt_->SetOption("colz");

  hitcluster_wt_xt_ =  hf.DefineTH2D(hdir, "hitcluster_wt_xt", "hit w(late) vs hit-cluster dt  for 2-cluster events, xt hits"+suffix, 200, 0., 1000., 100, 0., 1000.);
  hitcluster_wt_xt_->SetOption("colz");

  //----------------
  hitcluster_dt_nxt_ =  hf.DefineTH1D(hdir, "hitcluster_dt_nxt", "hit t(late)- cluster t(early)  for 2-cluster events, nxt hits"+suffix, 200, 0., 1000.);

  hitcluster_ww_nxt_ =  hf.DefineTH2D(hdir, "hitcluster_ww_nxt", "hit w(late) vs cluster maxw(early)  for 2-cluster events, nxt hits"+suffix, 100, 0., 1000., 100, 0., 1000.);
  hitcluster_ww_nxt_->SetOption("colz");

  hitcluster_wt_nxt_ =  hf.DefineTH2D(hdir, "hitcluster_wt_nxt", "hit w(late) vs hit-cluster dt  for 2-cluster events, nxt hits"+suffix, 200, 0., 1000., 100, 0., 1000.);
  hitcluster_wt_nxt_->SetOption("colz");

}

//================================================================
void HistXT3TwoCluster::fill(const WireClusterCollection& planeClusters) {
  if(planeClusters.size() == 2) {
    const unsigned iFirst = iearliestTimeCluster(planeClusters);

    const WireCluster& c1 = planeClusters[iFirst];
    const WireCluster& c2 = planeClusters[1-iFirst];

    double t1 = earliestTimeHit(c1)->time();
    double t2 = earliestTimeHit(c2)->time();

    twocluster_dt_->Fill(t2 - t1);
    twocluster_ww_->Fill(c1.maxTDCWidth(), c2.maxTDCWidth());
    twocluster_wt_->Fill(t2 - t1, c2.maxTDCWidth());

    twocluster_ss_->Fill(c1.numCells(), c2.numCells());

    for(unsigned i=0; i<c2.hits().size(); ++i) {
      const TDCHitWP& h2 = *c2.hits()[i];

      hitcluster_dt_all_->Fill(h2.time() - t1);
      hitcluster_ww_all_->Fill(c1.maxTDCWidth(), h2.width());
      hitcluster_wt_all_->Fill(h2.time() - t1, h2.width());

      if(h2.xtalk()) {
        hitcluster_dt_xt_->Fill(h2.time() - t1);
        hitcluster_ww_xt_->Fill(c1.maxTDCWidth(), h2.width());
        hitcluster_wt_xt_->Fill(h2.time() - t1, h2.width());
      }
      else {
        hitcluster_dt_nxt_->Fill(h2.time() - t1);
        hitcluster_ww_nxt_->Fill(c1.maxTDCWidth(), h2.width());
        hitcluster_wt_nxt_->Fill(h2.time() - t1, h2.width());
      }
    }
  }
}

//================================================================
void HistXT3::init(const std::string& hdir,
                   HistogramFactory &hf,
                   const ConfigFile& conf)
{
  trackCosth_ =  hf.DefineTH1D(hdir, "trackCosth", "track cos(theta)", 100, 0., 1.);
  trackAngle_ =  hf.DefineTH1D(hdir, "trackAngle", "track theta (degrees)", 100, 0., 90.);
  trackMomentum_ =  hf.DefineTH1D(hdir, "trackMomentum", "track momentum", 400, 0., 400.);

  inputHitMultiplicity23_ =  hf.DefineTH1D(hdir, "inputHitMultiplicity23", "hit multiplicity in DC23", 20, -0.5, 19.5);
  inputHitMultiplicity23xt_ =  hf.DefineTH1D(hdir, "inputHitMultiplicity23xt", "xt hit multiplicity in DC23", 20, -0.5, 19.5);
  inputHitMultiplicity23nxt_ =  hf.DefineTH1D(hdir, "inputHitMultiplicity23nxt", "nxt hit multiplicity in DC23", 20, -0.5, 19.5);

  hitWidthAll_ = hf.DefineTH1D(hdir, "hitWidthAll", "hit width for, all", 150, 0., 1500.);
  hitWidthXT_  = hf.DefineTH1D(hdir, "hitWidthXT",  "hit width, xt", 150, 0., 1500.);
  hitWidthNXT_ = hf.DefineTH1D(hdir, "hitWidthNXT", "hit width, nxt", 150, 0., 1500.);

  hitTrackTimeAll_ = hf.DefineTH1D(hdir, "hitTrackTimeAll", "hit time - track time, all time window hits", 240, -100., 1100.);

  inputWireMultiplicity23_ =  hf.DefineTH1D(hdir, "inputWireMultiplicity23", "wire multiplicity in DC23", 20, -0.5, 19.5);

  clusterMultiplicity23All_ =  hf.DefineTH1D(hdir, "clusterMultiplicity23All", "cluster multiplicity in DC23, all", 10, -0.5, 9.5);
  clusterMultiplicity23NXT_ =  hf.DefineTH1D(hdir, "clusterMultiplicity23NXT", "cluster multiplicity in DC23, nxt", 10, -0.5, 9.5);
  clusterMultiplicity23Clark_ =  hf.DefineTH1D(hdir, "clusterMultiplicity23Clark", "cluster multiplicity in DC23, clark", 10, -0.5, 9.5);

  singleHit23Dt_ =  hf.DefineTH1D(hdir, "dc23onehitdt", "t(hit)-t(track) for single hits in DC23", 240, -100., 1100.);
  singleHit23Width_ =  hf.DefineTH1D(hdir, "dc23onehitwidth", "hit width for single hits in DC23", 150, 0., 1500.);

  hits23MaxWidth_ =  hf.DefineTH1D(hdir, "dc23maxhitwidth", "max hit width in DC23", 150, 0., 1500.);
  hits23FirstWidth_ =  hf.DefineTH1D(hdir, "dc23firsthitwidth", "hit width for the earliest hit in DC23", 150, 0., 1500.);

  doubleClusterSizeAll_ =  hf.DefineTH2D(hdir, "doubleClusterSizeAll", "min vs max cluster size, all", 20, 0.5, 19.5, 5, -0.5, 4.5);
  doubleClusterSizeAll_->SetOption("colz");

  doubleClusterSizeNXT_ =  hf.DefineTH2D(hdir, "doubleClusterSizeNXT", "min vs max cluster size for 2-cluster events, nxt", 20, 0.5, 19.5, 5, -0.5, 4.5);
  doubleClusterSizeNXT_->SetOption("colz");

  doubleClusterSizeClark_ =  hf.DefineTH2D(hdir, "doubleClusterSizeClark", "min vs max cluster size for 2-cluster events, Clark", 20, 0.5, 19.5, 5, -0.5, 4.5);
  doubleClusterSizeClark_->SetOption("colz");

  //----------------
  twoClusterAll_.init(hdir+"/twoClusterAll", ", all", hf, conf);
  twoClusterNXT_.init(hdir+"/twoClusterNXT", ", nxt", hf, conf);
}

//================================================================
void HistXT3::fill(const EventClass& evt, int iTrack, const ClustersByPlane& globalProtonClusters) {
  // use same values as in windowing
  // NB: here these are w.r.t. the track instead of w.r.t. PC time
  const double winDCStart = -100.;
  const double winDCEnd = 1050.;

  const double costh= evt.costh[iTrack];
  const double angle = atan2(evt.sinth[iTrack], evt.costh[iTrack]) * 180/(4.*atan(1.));

  trackCosth_->Fill(costh);
  trackAngle_->Fill(angle);
  trackMomentum_->Fill(evt.ptot[iTrack]);

  const double trackTime = evt.hefit_time[iTrack];

  TDCHitWPPtrCollection dc23hitsAll, dc23hitsNXT;

  for(unsigned i=0; i< evt.dc_hits().size(); ++i) {
    if(evt.dc_hits()[i].plane() == 23) {
      const double dt = evt.dc_hits()[i].time() - trackTime;
      if((winDCStart < dt) && (dt < winDCEnd)) {
        dc23hitsAll.push_back(TDCHitWPPtr(evt.dc_hits(), i));
        hitWidthAll_->Fill(evt.dc_hits()[i].width());
        hitTrackTimeAll_->Fill(evt.dc_hits()[i].time() - trackTime);

        if(!evt.dc_hits()[i].xtalk()) {
          dc23hitsNXT.push_back(TDCHitWPPtr(evt.dc_hits(), i));

          hitWidthNXT_->Fill(evt.dc_hits()[i].width());
        }
        else {
          hitWidthXT_->Fill(evt.dc_hits()[i].width());
        }
      }
    }
  }

  inputHitMultiplicity23_->Fill(dc23hitsAll.size());
  int xtmult=0;
  for(unsigned i=0; i<dc23hitsAll.size(); ++i) {
    if(dc23hitsAll[i]->xtalk()) {
      ++xtmult;
    }
  }
  inputHitMultiplicity23xt_->Fill(xtmult);
  inputHitMultiplicity23nxt_->Fill(dc23hitsAll.size() - xtmult);


  if(dc23hitsAll.size() == 1) {
    const double dt = dc23hitsAll[0]->time() - trackTime;
    singleHit23Dt_->Fill(dt);
    singleHit23Width_->Fill(dc23hitsAll[0]->width());
  }

  if(!dc23hitsAll.empty()) {
    unsigned iwide = 0;
    double maxw = dc23hitsAll[iwide]->width();

    unsigned ifirst = 0;
    double tfirst = dc23hitsAll[ifirst]->time();

    for(unsigned i=1; i<dc23hitsAll.size(); ++i) {
      if(maxw < dc23hitsAll[i]->width()) {
        iwide = i;
        maxw = dc23hitsAll[iwide]->width();
      }
      if(dc23hitsAll[i]->time() < tfirst) {
        ifirst = i;
        tfirst = dc23hitsAll[ifirst]->time();
      }
    }

    hits23MaxWidth_->Fill(maxw);
    hits23FirstWidth_->Fill(dc23hitsAll[ifirst]->width());
  }

  //----------------------------------------------------------------
  ClustersByPlane clustersAll = constructPlaneClusters(44, dc23hitsAll);
  const WireClusterCollection& clusters23all = clustersAll[23];
  clusterMultiplicity23All_->Fill(clusters23all.size());

  int nwires = 0;
  for(unsigned i=0; i<clusters23all.size(); ++i) {
    nwires += clusters23all[i].numCells();
  }
  inputWireMultiplicity23_->Fill(nwires);

  if(clusters23all.size()==1) {
    doubleClusterSizeAll_->Fill(clusters23all[0].numCells(), 0.);
  }
  if(clusters23all.size()==2) {
    int s1 = clusters23all[0].numCells();
    int s2 = clusters23all[1].numCells();
    doubleClusterSizeAll_->Fill(std::max(s1, s2), std::min(s1, s2));
  }

  //----------------------------------------------------------------
  ClustersByPlane clustersNXT = constructPlaneClusters(44, dc23hitsNXT);
  const WireClusterCollection& clusters23nxt = clustersNXT[23];
  clusterMultiplicity23NXT_->Fill(clusters23nxt.size());

  if(clusters23nxt.size()==1) {
    doubleClusterSizeNXT_->Fill(clusters23nxt[0].numCells(), 0.);
  }
  if(clusters23nxt.size()==2) {
    int s1 = clusters23nxt[0].numCells();
    int s2 = clusters23nxt[1].numCells();
    doubleClusterSizeNXT_->Fill(std::max(s1, s2), std::min(s1, s2));
  }

  //----------------------------------------------------------------
  const WireClusterCollection& clusters23clark = globalProtonClusters[23 + 8];
  clusterMultiplicity23Clark_->Fill(clusters23clark.size());

  if(clusters23clark.size()==1) {
    doubleClusterSizeClark_->Fill(clusters23clark[0].numCells(), 0.);
  }
  if(clusters23clark.size()==2) {
    int s1 = clusters23clark[0].numCells();
    int s2 = clusters23clark[1].numCells();
    doubleClusterSizeClark_->Fill(std::max(s1, s2), std::min(s1, s2));
  }

  //----------------------------------------------------------------
  twoClusterAll_.fill(clusters23all);
  twoClusterNXT_.fill(clusters23nxt);

  //----------------------------------------------------------------

  // //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // static int count = 0;
  // if(count++ < 5) {
  //   std::cout<<"DC23 all hits = "<<clusters23all<<std::endl;
  //   std::cout<<"DC23 clark hits = "<<clusters23clark<<std::endl;
  //   //std::cout<<"globalProtonClusters = "<<clustersClark<<std::endl;
  // }
  // else {
  //   throw std::runtime_error("stop here");
  // }


}

//================================================================
