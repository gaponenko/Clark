// Andrei Gaponenko, 2013

#include "MuCapProtonWindow.h"

#include <algorithm>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

#include "PlaneRange.h"


//================================================================
namespace { // local helpers
  double getMinTime(const WireClusterCollection& coll) {
    assert(!coll.empty());
    double t = coll.front().hits().front()->time;
    for(WireClusterCollection::const_iterator c = coll.begin(); c!=coll.end(); ++c) {
      for(TDCHitWPPtrCollection::const_iterator h = c->hits().begin(); h != c->hits().end(); ++h) {
        if(t > (*h)->time) {
          t = (*h)->time;
        }
      }
    }
    return t;
  }
}

//================================================================
void MuCapProtonWindow::init(HistogramFactory &hf, const ConfigFile& conf) {
  maxPlane_ = conf.read<double>("MuCapture/ProtonWindow/maxPlane");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D("MuCapture/ProtonWindow", "protonCuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D("MuCapture/ProtonWindow", "protonCuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hNumClusters_ = hf.DefineTH2D("MuCapture/ProtonWindow", "numClustersVsPlane", "Num cluster vs plane", 56, 0.5, 56.5,  11, -0.5, 10.5);

  hStartPos_ = hf.DefineTH2D("MuCapture/ProtonWindow", "startPos", "Proton start V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);

  hStartOffset_ = hf.DefineTH2D("MuCapture/ProtonWindow", "startOffset", "Proton start offset V vs U (cell units)", 41, -5.125, 5.125,  41, -5.125, 5.125);

  hStartClusterSize_ = hf.DefineTH2D("MuCapture/ProtonWindow", "startClusterSize", "Proton start cluster V vs U width (cell units)", 25, -0.5, 24.5,  25, -0.5, 24.5);

  hLastPlane_ = hf.DefineTH1D("MuCapture/ProtonWindow", "lastPlane", "Proton last plane", 56, 0.5, 56.5);
  hProtonTime_ = hf.DefineTH1D("MuCapture/ProtonWindow", "protonTime", "Proton time", 1000, 0., 10000.);

  hwidthPCProtonWin_.init("MuCapture/pcWidthProton", "pcpwidth", 12, hf, conf);
  hwidthDCProtonWin_.init("MuCapture/dcWidthProton", "dcpwidth", 44, hf, conf);

  uvan_.init("MuCapture/ProtonWindow/UVAnalysis", hf, conf);
}

//================================================================
void MuCapProtonWindow::process(const ROOT::Math::XYPoint& muStopUV,
                                const TimeWindow& protonWindowPC,
                                const TimeWindow& protonWindowDC,
                                const TDCHitWPPtrCollection& unassignedDCHits,
                                const EventClass& evt
                                )
{
  EventCutNumber c = analyze(muStopUV, protonWindowPC, protonWindowDC, unassignedDCHits, evt);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
}

//================================================================
MuCapProtonWindow::EventCutNumber MuCapProtonWindow::
analyze(const ROOT::Math::XYPoint& muStopUV,
        const TimeWindow& protonWindowPC,
        const TimeWindow& protonWindowDC,
        const TDCHitWPPtrCollection& unassignedDCHits,
        const EventClass& evt
        )
{
  const TDCHitWPPtrCollection& protonPCHits = protonWindowPC.hits;
  const TDCHitWPPtrCollection& protonDCHits = protonWindowDC.hits;

  const ClustersByPlane clustersPC = constructPlaneClusters(12, protonPCHits);
  const ClustersByPlane clustersDC = constructPlaneClusters(44, protonDCHits);
  const ClustersByPlane global = globalPlaneClusters(clustersPC, clustersDC);

  const PlaneRange gr = findPlaneRange(global);
  if(gr.min < 29) {
    return CUT_UPSTREAM;
  }

  if(clustersPC[7].empty()) {
    return CUT_NOPC7;
  }

  for(int plane=1; plane <= 56; ++plane) {
    hNumClusters_->Fill(plane, global[plane].size());
  }

  if(!gr.noGaps) {
    return CUT_RANGE_GAPS;
  }

  uvan_.process(evt,  protonWindowPC.tstart, global, muStopUV);

  hLastPlane_->Fill(gr.max);
  if(gr.max > maxPlane_) {
    return CUT_MAX_PLANE;
  }

//  if(clustersPC[7].size() > 1) {
//    std::cout<<"pc 7 clusters:\n"<<clustersPC[7]<<std::endl;
//  }

  // If we use window start time here, longer proton tracks will get sharper timing
  // because of more hits.  Restrict the amount of data we use for proton time to
  // just PC7, so that all tracks get equally bad timing.
  hProtonTime_->Fill(getMinTime(clustersPC[7]));

  //std::cout<<"accepted: "<<evt.nrun<<"\t"<<evt.nevt<<"\t"<<gr.max<<std::endl;

  hwidthPCProtonWin_.fill(protonPCHits);
  hwidthDCProtonWin_.fill(protonDCHits);

  return CUTS_ACCEPTED;
}
