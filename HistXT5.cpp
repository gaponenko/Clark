// Andrei Gaponenko, 2013

#include "HistXT5.h"

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
void HistXT5Ana::init(const std::string& hdir,
                      const std::string& suffix,
                      HistogramFactory& hf,
                      const ConfigFile& conf)
{
  clustermult_ =  hf.DefineTH2D(hdir, "cluster_multiplicity", "DC24 vs DC23 cluster multiplicity, "+suffix, 10, -0.5, 9.5, 10, -0.5, 9.5);
  clustermult_->SetOption("colz");
  clustermult_->GetXaxis()->SetTitle("DC23 cluster multiplicity");
  clustermult_->GetYaxis()->SetTitle("DC24 cluster multiplicity");

  clustersize_ =  hf.DefineTH2D(hdir, "cluster_size", "DC24 vs DC23 max cluster size, "+suffix, 10, -0.5, 9.5, 10, -0.5, 9.5);
  clustersize_->SetOption("colz");
  clustersize_->GetXaxis()->SetTitle("DC23 max cluster size");
  clustersize_->GetYaxis()->SetTitle("DC24 max cluster size");

  logic_cm_ = hf.DefineTH2D(hdir, "logic_cm", "DC24 vs DC23 cluster mult. logic, "+suffix, 2, -0.5, 1.5, 2, -0.5, 1.5);
  logic_cm_->SetOption("colz,text");
  logic_cm_->GetXaxis()->SetTitle("DC23");
  logic_cm_->GetYaxis()->SetTitle("DC24");

  logic_cs_ = hf.DefineTH2D(hdir, "logic_cs", "DC24 vs DC23 cluster size logic, "+suffix, 2, -0.5, 1.5, 2, -0.5, 1.5);
  logic_cs_->SetOption("colz,text");
  logic_cs_->GetXaxis()->SetTitle("DC23");
  logic_cs_->GetYaxis()->SetTitle("DC24");

  logic_cmcs_ = hf.DefineTH2D(hdir, "logic_cmcs", "DC24 vs DC23 cm and cs logic, "+suffix, 2, -0.5, 1.5, 2, -0.5, 1.5);
  logic_cmcs_->SetOption("colz,text");
  logic_cmcs_->GetXaxis()->SetTitle("DC23");
  logic_cmcs_->GetYaxis()->SetTitle("DC24");
}

//================================================================
namespace {
  int maxClusterSize(const WireClusterCollection& clusters) {
    int res = 0;
    for(unsigned i=0; i<clusters.size(); ++i) {
      if(res < clusters[i].numCells()) {
        res = clusters[i].numCells();
      }
    }
    return res;
  }
}

void HistXT5Ana::fill(const TDCHitWPPtrCollection& selectHits) {
  ClustersByPlane planes = constructPlaneClusters(44, selectHits);
  const WireClusterCollection& clusters23 = planes[23];
  const WireClusterCollection& clusters24 = planes[24];

  clustermult_->Fill(clusters23.size(), clusters24.size());

  const int ms23 = maxClusterSize(clusters23);
  const int ms24 = maxClusterSize(clusters24);

  clustersize_->Fill(ms23, ms24);

  bool lcm23 = (clusters23.size() < 2);
  bool lcm24 = (clusters24.size() < 2);
  logic_cm_->Fill(!lcm23, !lcm24);

  bool lcs23 = (ms23 < 3);
  bool lcs24 = (ms24 < 3);
  logic_cs_->Fill(!lcs23, !lcs24);

  logic_cmcs_->Fill(!(lcm23&&lcs23), !(lcm24&&lcs24));
}

//================================================================
void HistXT5::init(const std::string& hdir,
                          HistogramFactory &hf,
                          const DetectorGeo& geom,
                          const ConfigFile& conf)
{
  hpcos_all_  = hf.DefineTH2D(hdir, "pcos_all", "pcos, all", 500, 0., 500., 100, -1., +1.);
  hpcos_all_->SetOption("colz");
  hpcos_passed_  = hf.DefineTH2D(hdir, "pcos_passed", "pcos, all", 500, 0., 500., 100, -1., +1.);
  hpcos_passed_->SetOption("colz");

  trackAngle_all_ =  hf.DefineTH1D(hdir, "trackAngle_all", "track theta (degrees), before cut", 100, 0., 90.);
  trackAngle_passed_ =  hf.DefineTH1D(hdir, "trackAngle_passed", "track theta (degrees), after cut", 100, 0., 90.);

  xtanaAll_.init(hdir+"/anaAll", ", all", hf, conf);
  xtanaNXT_.init(hdir+"/anaNXT", ", nxt", hf, conf);
}

//================================================================
void HistXT5::fill(const EventClass& evt, int iTrack) {
  // use same values as in windowing
  // NB: here these are w.r.t. the track instead of w.r.t. PC time
  const double winDCStart = -100.;
  const double winDCEnd = 1050.;

  const double costh= evt.costh[iTrack];
  const double angle = atan2(evt.sinth[iTrack], evt.costh[iTrack]) * 180/(4.*atan(1.));
  const double ptot = evt.ptot[iTrack];

  hpcos_all_->Fill(ptot,costh);
  trackAngle_all_->Fill(angle);

  if(angle < 40.) {
    hpcos_passed_->Fill(ptot,costh);
    trackAngle_passed_->Fill(angle);

    // Select DC hits in the track's time window
    const double trackTime = evt.hefit_time[iTrack];

    TDCHitWPPtrCollection selectHitsAll, selectHitsNXT;

    for(unsigned i=0; i< evt.dc_hits().size(); ++i) {
      const double dt = evt.dc_hits()[i].time() - trackTime;
      if((winDCStart < dt) && (dt < winDCEnd)) {

        // We are interested only in DC23 and DC24 hits
        // skip others for efficiency
        if((23 <= evt.dc_hits()[i].plane())&&(evt.dc_hits()[i].plane() <= 24)) {
          selectHitsAll.push_back(TDCHitWPPtr(evt.dc_hits(), i));
          if(!evt.dc_hits()[i].xtalk()) {
            selectHitsNXT.push_back(TDCHitWPPtr(evt.dc_hits(), i));
          }
        }
      }
    } // select hits

    xtanaAll_.fill(selectHitsAll);
    xtanaNXT_.fill(selectHitsNXT);

  } // angle < 40
}

//================================================================
