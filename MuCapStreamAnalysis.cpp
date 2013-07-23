// Andrei Gaponenko, 2013

#include "MuCapStreamAnalysis.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "PlaneRange.h"
#include "EventList.h"
#include "MuCapUtilities.h"

//================================================================
void MuCapStreamAnalysis::init(HistogramFactory &hf, const std::string& hdir,
                               const DetectorGeo& geom, const ConfigFile& conf,
                               TimeWindow::StreamType cutWinStream, double cutWinTimeMin)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");
  cutStream_ = cutWinStream;

  cutBeamVetoMaxPCplanes_ = conf.read<double>("MuCapture/StreamAnalysis/cutBeamVetoMaxPCplanes");
  cutWinTimeMin_ = cutWinTimeMin;
  cutWinTimeMax_ = conf.read<double>("MuCapture/cutWinTimeMax");
  cutMultiwinNextdt_ = conf.read<double>("MuCapture/StreamAnalysis/cutMultiwinNextdt");

  cutZContainedNumToCheck_ = conf.read<unsigned>("MuCapture/StreamAnalysis/cutZContainedNumToCheck");
  cutZContainedMaxHitPlanes_ = conf.read<unsigned>("MuCapture/StreamAnalysis/cutZContainedMaxHitPlanes");
  cutRextMax_ = conf.read<double>("MuCapture/StreamAnalysis/cutRextMax");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "protonCuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D(hdir, "protonCuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hBeamVetoNumHitPlanes_ = hf.DefineTH1D(hdir, "beamVetoNumHitPlanes", "beamVetoNumHitPlanes", 6, -0.5, 5.5);
  hHitPCsAterBeamVeto_ = hf.DefineTH1D(hdir, "hitUpsteamPCsAfterBeamVeto", "hitUpsteamPCsAfterBeamVeto", 6, -0.5, 5.5);

  hWindowTimeBefore_ = hf.DefineTH1D(hdir, "windowTimeBeforeCut", "Proton window start time", 1000, 0., 10000.);
  hWindowTimeAfter_ = hf.DefineTH1D(hdir, "windowTimeAfterCut", "Proton window start time after the cut", 1000, 0., 10000.);

  hNumAfterTrigWindows_ = hf.DefineTH1D(hdir, "numAfterTrigTimeWindows", "numAfterTrigTimeWindows", 10, -0.5, 9.5);

  hWindow2Time_ = hf.DefineTH1D(hdir, "window2Time", "Second window start time", 1000, 0., 10000.);
  hWindow2dt_ = hf.DefineTH1D(hdir, "window2dt", "Second window time - proton win time", 1000, 0., 10000.);

  hZContaintedNumHitPlanesUp_ = hf.DefineTH1D(hdir, "hZContaintedNumHitPlanesUp", "hZContaintedNumHitPlanesUp", 10, -0.5, 9.5);
  hZContaintedNumHitPlanesDn_ = hf.DefineTH1D(hdir, "hZContaintedNumHitPlanesDn", "hZContaintedNumHitPlanesDn", 10, -0.5, 9.5);

  hNumPC7Clusters_ = hf.DefineTH1D(hdir, "numPC7Clusters", "Num PC7 clusters before cut", 10, -0.5, 9.5);
  hNumPC7WiresVsClusters_ = hf.DefineTH2D(hdir, "numPC7WiresVsClusters",
                                          "Num PC7 wires vs clusters before cut",
                                          10, -0.5, 9.5, 40, -0.5, 39.5);
  hNumPC7WiresVsClusters_->SetOption("colz");

  hLastPlaneLoose_ = hf.DefineTH1D(hdir, "lastPlaneLoose", "Loose proton last plane", 56, 0.5, 56.5);

  //----------------------------------------------------------------
  uvan_.init(hdir+"/UVAnalysis", hf, conf, TimeWindow::DOWNSTREAM, cutWinTimeMin - 100./*FIXME*/);
  hRangeDIO_.init(hdir+"/rangeDIO", hf, geom, conf);
  hdriftPCFiltered_.init(hdir+"/driftTimePCFiltered", hf, geom.numPCs(), 1000./*ns*/,
                         conf.read<double>("MuCapture/HistDriftTime/cutEffTrackHitDtPC"),
                         conf);

  //----------------------------------------------------------------
  hhsZContained_.init(hdir+"/hsZContained", hf, geom, conf);
  hRangeAfterPC7Cuts_.init(hdir+"/rangeAfterPC7Cuts", hf, geom, conf);

  //----------------------------------------------------------------
  hwidthPCTightDIO_.init(hdir+"/pcWidthTightDIO", "pcpwidth", 12, hf, conf);
  hwidthDCTightDIO_.init(hdir+"/dcWidthTightDIO", "dcpwidth", 44, hf, conf);
  hwidthPCLooseProtons_.init(hdir+"/pcWidthLooseProtons", "pcpwidth", 12, hf, conf);
  hwidthDCLooseProtons_.init(hdir+"/dcWidthLooseProtons", "dcpwidth", 44, hf, conf);
  hwidthPCTightProtons_.init(hdir+"/pcWidthTightProtons", "pcpwidth", 12, hf, conf);
  hwidthDCTightProtons_.init(hdir+"/dcWidthTightProtons", "dcpwidth", 44, hf, conf);

  hcdio_.init(hf, hdir+"/clDIO", geom, conf, cutStream_);
  hcLooseProtons_.init(hf, hdir+"/clLooseProtons", geom, conf, cutStream_);
  hcTightProtons_.init(hf, hdir+"/clTightProtons", geom, conf, cutStream_);
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
  if(afterTrigGlobalClusters.empty()) {
    return CUT_NOHITS;
  }

  const unsigned iProtonWin = 1 + wres.iTrigWin;
  const TimeWindow& protonWindow = wres.windows[iProtonWin];
  const ClustersByPlane& protonGlobalClusters = afterTrigGlobalClusters[0];
  const PlaneRange protonDnClusters = findDownstreamPlaneRange(protonGlobalClusters);

  //----------------------------------------------------------------
  // Veto accidental beam particles (also upstream DIOs and some protons)

  int numBeamVetoHitPlanes(0);
  for(int i=1; i<=4; ++i) {
    if(!protonGlobalClusters[i].empty()) {
      ++numBeamVetoHitPlanes;
    }
  }

  hBeamVetoNumHitPlanes_->Fill(numBeamVetoHitPlanes);
  if(numBeamVetoHitPlanes > cutBeamVetoMaxPCplanes_) {
    return CUT_BEAM_VETO;
  }

  for(int i=1; i<=4; ++i) {
    if(!protonGlobalClusters[i].empty()) {
      hHitPCsAterBeamVeto_->Fill(i);
    }
  }

  //----------------------------------------------------------------
  // Trig time is 0, dt from that rather than from less precise trigWin time
  hWindowTimeBefore_->Fill(protonWindow.tstart);
  if( (protonWindow.tstart < cutWinTimeMin_) || (cutWinTimeMax_ < protonWindow.tstart)) {
    return CUT_WINTIME;
  }
  hWindowTimeAfter_->Fill(protonWindow.tstart);

  //----------------------------------------------------------------
  hNumAfterTrigWindows_->Fill(afterTrigGlobalClusters.size());
  if(afterTrigGlobalClusters.size() > 1) {
    const TimeWindow& win1 = wres.windows[wres.iTrigWin + 1];
    const TimeWindow& win2 = wres.windows[wres.iTrigWin + 2];
    const double dt2 = win2.tstart - win1.tstart;
    hWindow2Time_->Fill(win2.tstart);
    hWindow2dt_->Fill(dt2);
    if(dt2 < cutMultiwinNextdt_) {
      return CUT_MULTIWIN_NEXTDT;
    }
  }

  //----------------------------------------------------------------
  // UVAnalysis: DIOs for normalization
  const int iDIO = uvan_.process(evt, muStopUV);
  if(iDIO >= 0) {
    hwidthPCTightDIO_.fill(protonWindow.pcHits);
    hwidthDCTightDIO_.fill(protonWindow.dcHits);
    hcdio_.fill(protonGlobalClusters);

    hRangeDIO_.fill(protonDnClusters);
    hdriftPCFiltered_.fill(evt, iDIO, protonWindow.pcHits);
  }

  //----------------------------------------------------------------
  // Z containment check
  int numZVetoHitsUp(0), numZVetoHitsDn(0);
  for(unsigned i=1; i<=cutZContainedNumToCheck_; ++i) {
    if(!protonGlobalClusters[i].empty()) {
      ++numZVetoHitsUp;
    }
    if(!protonGlobalClusters[protonGlobalClusters.size() - i].empty()) {
      ++numZVetoHitsDn;
    }
  }

  hZContaintedNumHitPlanesUp_->Fill(numZVetoHitsUp);
  hZContaintedNumHitPlanesDn_->Fill(numZVetoHitsDn);
  if((numZVetoHitsUp > cutZContainedMaxHitPlanes_) || (numZVetoHitsDn > cutZContainedMaxHitPlanes_)) {
    return CUT_Z_CONTAINED;
  }

  hhsZContained_.fill(protonGlobalClusters);

  //----------------------------------------------------------------
  const int PC7GlobalIndex = 29;
  hNumPC7Clusters_->Fill(protonGlobalClusters[PC7GlobalIndex].size());
  hNumPC7WiresVsClusters_->Fill(protonGlobalClusters[PC7GlobalIndex].size(),
                                MuCapUtilities::numWires(protonGlobalClusters[PC7GlobalIndex]));

  if(protonGlobalClusters[PC7GlobalIndex].empty()) {
    return CUT_PC7_HIT;
  }

  //----------------------------------------------------------------
//  if() {
//    return CUT_PC7_COORDINATE;
//  }

  hRangeAfterPC7Cuts_.fill(protonDnClusters);

  //  //----------------------------------------------------------------
  //  // CUTS_LOOSE_PROTONS are passed at this point
  //
  const PlaneRange gr = findPlaneRange(protonGlobalClusters);
  hLastPlaneLoose_->Fill(gr.max());

  hwidthPCLooseProtons_.fill(protonWindow.pcHits);
  hwidthDCLooseProtons_.fill(protonWindow.dcHits);
  hcLooseProtons_.fill(protonGlobalClusters);


  // FIXME: just do the downstream analysis, but at earlier times
  // Can use the already-studies TDC width cut to select protons
  // Can allow some upstream hits...  Note that upstream DC hits
  // would be eaten by the muon window for early times anyway.


  //  if(!gr.noGaps) {
  //    return CUT_RANGE_GAPS;
  //  }
  //
  //  if(doMCTruth_ && (evt.nmcvtx == 2)) {
  //    hLastPlaneVsMCPstart_->Fill(evt.mcvertex_ptot[0], gr.max);
  //  }
  //
  //
  //  if(doMCTruth_) {
  //    htruthLoose_.fill(evt);
  //  }
  //
  //  //hpw_.fill(global, evt);
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
