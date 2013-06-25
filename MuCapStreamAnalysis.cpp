// Andrei Gaponenko, 2013

#include "MuCapStreamAnalysis.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "PlaneRange.h"
#include "EventList.h"
#include "MuCapUtilities.h"

//================================================================
namespace {
  template<class Container> bool nonEmpty(const Container& c) {
    return !c.empty();
  }
}

//================================================================
void MuCapStreamAnalysis::init(HistogramFactory &hf, const std::string& hdir,
                               const DetectorGeo& geom, const ConfigFile& conf,
                               TimeWindow::StreamType cutWinStream, double cutWinTimeMin)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");
  cutStream_ = cutWinStream;
  cutWinTimeMin_ = cutWinTimeMin;
  cutWinTimeMax_ = conf.read<double>("MuCapture/cutWinTimeMax");
  cutZContainedNumPlanes_ = conf.read<unsigned>("MuCapture/StreamAnalysis/cutZContainedNumPlanes");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "protonCuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D(hdir, "protonCuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hMultiWindowPCHits_ = hf.DefineTH2D(hdir, "multiWindowPCHits", "after trig PC hits2 vs hits1",
                                      20, -0.5, 19.5, 20, -0.5, 19.5);
  hMultiWindowPCHits_->SetOption("colz");

  hMultiWindowHits_ = hf.DefineTH2D(hdir, "multiWindowHits", "after trig hits2 vs hits1",
                                    140, -0.5, 139.5, 140, -0.5, 139.5);
  hMultiWindowHits_->SetOption("colz");

  hMultiWindowPCPlanes_ = hf.DefineTH2D(hdir, "multiWindowPCPlanes", "after trig PC planes2 vs planes1",
                                        13, -0.5, 12.5, 13, -0.5, 12.5);
  hMultiWindowPCPlanes_->SetOption("colz");

  hMultiWindowPCPlanesmm_ = hf.DefineTH2D(hdir, "multiWindowPCPlanesmm", "min vs max PC planes in after trig win1, 2",
                                          13, -0.5, 12.5, 13, -0.5, 12.5);
  hMultiWindowPCPlanesmm_->SetOption("colz");

  hMultiWindowPlanesmm_ = hf.DefineTH2D(hdir, "multiWindowPlanesmm", "min vs max num planes in after-trig windows 1, 2",
                                    57, -0.5, 56.5, 57, -0.5, 56.5);
  hMultiWindowPlanesmm_->SetOption("colz");

  hMultiWindowRanges_ = hf.DefineTH2D(hdir, "multiWindowRanges", "min vs max num ranges in after-trig windows 1, 2",
                                    10, -0.5, 9.5, 10, -0.5, 9.5);
  hMultiWindowRanges_->SetOption("colz");

  hOccupancyPCwin1_.init(hdir, "multiWindowPCMap1", 12, 160, hf, conf);
  hOccupancyPCwin2_.init(hdir, "multiWindowPCMap2", 12, 160, hf, conf);

  hNumAfterTrigWindows_ = hf.DefineTH1D(hdir, "numAfterTrigTimeWindows", "numAfterTrigTimeWindows", 10, -0.5, 9.5);
  hWindowTimeBefore_ = hf.DefineTH1D(hdir, "windowTimeBeforeCut", "Proton window start time", 1000, 0., 10000.);
  hWindowTimeAfter_ = hf.DefineTH1D(hdir, "windowTimeAfterCut", "Proton window start time after the cut", 1000, 0., 10000.);

  hhsAfterTimeCuts_.init(hdir+"/hsAfterTimeCuts", hf, geom, conf);

  hNumVetoHits_ = hf.DefineTH1D(hdir, "numZContainVetoHits", "Num hits in veto planes", 10, -0.5, 9.5);

  hhsZContained_.init(hdir+"/hsZContained", hf, geom, conf);

  //----------------------------------------------------------------
  uvan_.init(hdir+"/UVAnalysis", hf, conf, TimeWindow::MIXED);
  hwidthPCTightProtons_.init(hdir+"/pcWidthTightProtons", "pcpwidth", 12, hf, conf);
  hwidthDCTightProtons_.init(hdir+"/dcWidthTightProtons", "dcpwidth", 44, hf, conf);
  hwidthPCTightDIO_.init(hdir+"/pcWidthTightDIO", "pcpwidth", 12, hf, conf);
  hwidthDCTightDIO_.init(hdir+"/dcWidthTightDIO", "dcpwidth", 44, hf, conf);
  hcLooseProtons_.init(hf, hdir+"/TDCLooseProtons", geom, conf, cutStream_);
  hcTightProtons_.init(hf, hdir+"/TDCTightProtons", geom, conf, cutStream_);
  hcdio_.init(hf, hdir+"/TDCDIO", geom, conf, cutStream_);
}

//================================================================
void MuCapStreamAnalysis::process(const EventClass& evt,
                                  const TimeWindowingResults& wres,
                                  const ROOT::Math::XYPoint& muStopUV,
                                  const std::vector<ClustersByPlane>& afterTrigGlobalClusters)
{
  EventCutNumber c = analyze(evt, wres, muStopUV, afterTrigGlobalClusters);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }

  if(gEventList.requested(evt)) {
    const int cbin = h_cuts_p->GetXaxis()->FindFixBin(c);
    std::cout<<__func__<<": run "<<evt.nrun<<" event "<<evt.nevt
             <<" status "<<c<<": "<<h_cuts_p->GetXaxis()->GetBinLabel(cbin)
             <<std::endl;
  }
}

//================================================================
MuCapStreamAnalysis::EventCutNumber MuCapStreamAnalysis::
analyze(const EventClass& evt,
        const TimeWindowingResults& wres,
        const ROOT::Math::XYPoint& muStopUV,
        const std::vector<ClustersByPlane>& afterTrigGlobalClusters)
{
  //----------------------------------------------------------------
  hNumAfterTrigWindows_->Fill(afterTrigGlobalClusters.size());
  if(afterTrigGlobalClusters.empty()) {
    return CUT_NOHITS;
  }

  if(afterTrigGlobalClusters.size() > 1) {
    const TimeWindow& win1 = wres.windows[wres.iTrigWin + 1];
    const TimeWindow& win2 = wres.windows[wres.iTrigWin + 2];

    hMultiWindowPCHits_->Fill(win1.pcHits.size(), win2.pcHits.size());

    hMultiWindowHits_->Fill(win1.pcHits.size() + win1.dcHits.size(),
                            win2.pcHits.size() + win2.dcHits.size());

    unsigned const planes1 = MuCapUtilities::numPlanes(win1.pcHits);
    unsigned const planes2 = MuCapUtilities::numPlanes(win2.pcHits);
    hMultiWindowPCPlanes_->Fill(planes1, planes2);
    hMultiWindowPCPlanesmm_->Fill(std::max(planes1,planes2), std::min(planes1, planes2));

    hOccupancyPCwin1_.fill(win1.pcHits);
    hOccupancyPCwin2_.fill(win2.pcHits);

    if((cutWinTimeMin_ < win1.tstart) && (win2.tstart <= cutWinTimeMax_ )) {
      // Avoid DC overlaps for these plots
      if(win2.tstart - win1.tstart > 1050.) {

        const unsigned n1 = std::count_if(afterTrigGlobalClusters[0].begin(), afterTrigGlobalClusters[0].end(), nonEmpty<WireClusterCollection>);
        const unsigned n2 = std::count_if(afterTrigGlobalClusters[1].begin(), afterTrigGlobalClusters[1].end(), nonEmpty<WireClusterCollection>);
        hMultiWindowPlanesmm_->Fill(std::max(n1,n2), std::min(n1, n2));

        PlaneRange r1 = findPlaneRange(afterTrigGlobalClusters[0]);
        PlaneRange r2 = findPlaneRange(afterTrigGlobalClusters[1]);
        hMultiWindowRanges_->Fill(std::max(r1.segments().size(), r2.segments().size()),
                                  std::min(r1.segments().size(), r2.segments().size()));

        //std::cout<<"t1 = "<<win1.tstart<<", r1 = "<<r1<<",\tt2 = "<<win2.tstart<<", r2="<<r2<<std::endl;
      }
    }
  }

  //----------------------------------------------------------------
  if(afterTrigGlobalClusters.size() > 1) {
    return CUT_MULTIWIN;
  }

  const unsigned iProtonWin = 1 + wres.iTrigWin;
  const TimeWindow& protonWindow = wres.windows[iProtonWin];
  const ClustersByPlane& protonGlobalClusters = afterTrigGlobalClusters[0];

  //----------------------------------------------------------------
  // Trig time is 0, dt from that rather than from less precise trigWin time
  hWindowTimeBefore_->Fill(protonWindow.tstart);
  if( (protonWindow.tstart < cutWinTimeMin_) || (cutWinTimeMax_ < protonWindow.tstart)) {
    return CUT_WINTIME;
  }
  hWindowTimeAfter_->Fill(protonWindow.tstart);

  // Count the number of planes hit up and dn of the target
  int nHitPlanesUp(0), nHitPlanesDn(0);
  for(unsigned i=1; i<protonGlobalClusters.size(); ++i) {
    if(!protonGlobalClusters[i].empty()) {
      ++(i <= protonGlobalClusters.size()/2 ? nHitPlanesUp : nHitPlanesDn);
    }
  }
  //std::cout<<"AG: nHitPlanesUp = "<<nHitPlanesUp<<", nHitPlanesDn = "<<nHitPlanesDn<<std::endl;
  hhsAfterTimeCuts_.fill(protonGlobalClusters);

  const unsigned numDIO = uvan_.process(evt, protonWindow.tstart, protonGlobalClusters, muStopUV);
  if(numDIO) {
    hwidthPCTightDIO_.fill(protonWindow.pcHits);
    hwidthDCTightDIO_.fill(protonWindow.dcHits);
    hcdio_.fill(protonGlobalClusters);
  }

  //----------------------------------------------------------------
  // Count the number of hit planes in the veto region
  int numVetoHits(0);
  for(unsigned i=1; i<=cutZContainedNumPlanes_; ++i) {
    if(!protonGlobalClusters[i].empty()) {
      ++numVetoHits;
    }
    if(!protonGlobalClusters[protonGlobalClusters.size() - i].empty()) {
      ++numVetoHits;
    }
  }

  hNumVetoHits_->Fill(numVetoHits);
  if(numVetoHits) {
    return CUT_Z_CONTAINED;
  }


  // Here we have a sample of events with charged particle emission from mu capture
  // (plus some DIOs hitting the glass)
  // Normalize to number of captures in accepted deltaT (corrected for accidental effects)
  // and get a measurement of charged particles per capture.
  hhsZContained_.fill(protonGlobalClusters);


  // FIXME: just do the downstream analysis, but at earlier times
  // Can use the already-studies TDC width cut to select protons
  // Can allow some upstream hits...  Note that upstream hits
  // would be eaten by the muon window for early times anyway.


//  if(protonWindow.stream != cutStream_) {
//    return CUT_STREAM;
//  }

//  if(clustersPC[7].empty()) {
//    return CUT_NOPC7;
//  }
//
//  for(int plane=1; plane <= 56; ++plane) {
//    hNumClusters_->Fill(plane, global[plane].size());
//  }
//
//  if(!gr.noGaps) {
//    return CUT_RANGE_GAPS;
//  }
//
//  const unsigned numDIO = uvan_.process(evt,  protonWindow.tstart, global, muStopUV);
//  if(numDIO) {
//    hwidthPCTightDIO_.fill(protonPCHits);
//    hwidthDCTightDIO_.fill(protonDCHits);
//    hcdio_.fill(global);
//  }
//
//  hLastPlane_->Fill(gr.max);
//  if(doMCTruth_ && (evt.nmcvtx == 2)) {
//    hLastPlaneVsMCPstart_->Fill(evt.mcvertex_ptot[0], gr.max);
//  }
//
//  if(gr.max > cutMaxPlane_) {
//    return CUT_MAX_RANGE;
//  }
//
//  //----------------------------------------------------------------
//  // CUTS_LOOSE_PROTONS are passed at this point
//
//  hcLooseProtons_.fill(global);
//
//  if(doMCTruth_) {
//    htruthLoose_.fill(evt);
//  }
//
//  //hpw_.fill(global, evt);
//
//  //----------------------------------------------------------------
//  if(gr.max < 30) {
//    return CUT_NOPC8;
//  }
//
//  if(doMCTruth_) {
//    htruthPC8_.fill(evt);
//  }
//
//  //----------------------------------------------------------------
//  // the containment check requires at least 2U and 2V planes
//  if(gr.max < 32) {
//    return CUT_MIN_RANGE;
//  }
//
//  if(doMCTruth_) {
//    htruthMinRange_.fill(evt);
//  }
//
//  const double rext = rcheckProtonCandidates_.rmax(gr.max, global);
//  hCCRvsPlaneProtons_->Fill(gr.max, rext);
//  if(doMCTruth_) {
//    hrtruth_.fill(evt, gr.max, rext);
//  }
//  if(rext > cutRextMax_) {
//    return CUT_REXT;
//  }
//
//  if(doMCTruth_) {
//    htruthTight_.fill(evt);
//  }
//
//  hwidthPCTightProtons_.fill(protonPCHits);
//  hwidthDCTightProtons_.fill(protonDCHits);
//  hcTightProtons_.fill(global);

  return CUTS_TIGHT_PROTONS;
}
