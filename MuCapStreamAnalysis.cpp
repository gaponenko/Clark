// Andrei Gaponenko, 2013

#include "MuCapStreamAnalysis.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "PlaneRange.h"
#include "EventList.h"
#include "MuCapUtilities.h"

//================================================================
void MuCapStreamAnalysis::init(HistogramFactory &hf, const std::string& hdir,
                               const DetectorGeo& geom, const ConfigFile& conf,
                               RecoResMuCapTrk *resTrk,
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
  cutPC7MaxDistanceToMuStop_ = conf.read<double>("MuCapture/StreamAnalysis/cutPC7MaxDistanceToMuStop");
  cutRextMax_ = conf.read<double>("MuCapture/StreamAnalysis/cutRextMax");

  commonSkimOutFileName_ = conf.read<std::string>("MuCapture/StreamAnalysis/commonSkimOutFileName", "");
  if(!commonSkimOutFileName_.empty()) {
    commonSkimOutFile_.open(commonSkimOutFileName_.c_str());
    if(!commonSkimOutFile_) {
      throw std::runtime_error("Error opening output file "+commonSkimOutFileName_);
    }
  }

  //----------------------------------------------------------------
  geom_ = &geom;

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

  hPC7DistanceToMuStop_ = hf.DefineTH1D(hdir, "pc7DistanceToMuStop", "pc7DistanceToMuStop",
                                        //101, -0.0125, 2.5125); //cm
                                        200, 0, 5.); //cm

  hLastPlaneLoose_ = hf.DefineTH1D(hdir, "lastPlaneLoose", "Loose proton last plane", 56, 0.5, 56.5);

  //----------------------------------------------------------------
  uvan_.init(hdir+"/UVAnalysis", hf, conf, TimeWindow::DOWNSTREAM, cutWinTimeMin - 100./*FIXME*/);
  muCapTrkHF_.init(hdir+"/TrkHF", hf, conf, TimeWindow::DOWNSTREAM,
                   conf.read<double>("MuCapture/TrkAnalysisHF/cutTimeMin") /*FIXME*/,
                   resTrk);

  hRangeDIO_.init(hdir+"/rangeDIO", hf, geom, conf);
  hdriftPCFiltered_.init(hdir+"/driftTimePCFiltered", hf, geom.numPCs(), 1000./*ns*/,
                         conf.read<double>("MuCapture/HistDriftTime/cutEffTrackHitDtPC"),
                         conf);

  //----------------------------------------------------------------
  hhsZContained_.init(hdir+"/hsZContained", hf, geom, conf);
  hhsLooseProtons_.init(hdir+"/hsLooseProtons", hf, geom, conf);
  hRangeAfterPC7Cuts_.init(hdir+"/rangeAfterPC7Cuts", hf, geom, conf);

  //----------------------------------------------------------------
  hwidthPCDIO_.init(hdir+"/pcWidthDIO", "pcpwidth", 12, hf, conf);
  hwidthDCDIO_.init(hdir+"/dcWidthDIO", "dcpwidth", 44, hf, conf);
  hwidthPCLooseProtons_.init(hdir+"/pcWidthLooseProtons", "pcpwidth", 12, hf, conf);
  hwidthDCLooseProtons_.init(hdir+"/dcWidthLooseProtons", "dcpwidth", 44, hf, conf);

  hcdio_.init(hf, hdir+"/clDIO", geom, conf, cutStream_);
  hcLooseProtons_.init(hf, hdir+"/clLooseProtons", geom, conf, cutStream_);

  hfLoose_.init(hf, hdir+"/finalLoose", geom, conf);

  //----------------------------------------------------------------
  if(doMCTruth_) {
    htruthAfterNoHits_.init(hf, hdir + "/TruthAfterNoHits", conf);
    htruthAfterBeamVeto_.init(hf, hdir + "/TruthAfterBeamVeto", conf);
    htruthAfterMutliwinTime_.init(hf, hdir + "/TruthAfterMultiwinTime", conf);
    htruthAfterWinTime_.init(hf, hdir + "/TruthAfterWinTime", conf);
    htruthContained_.init(hf, hdir + "/TruthContained", conf);
    htruthLoose_.init(hf, hdir + "/TruthLoose", conf);
    htruth2planes_.init(hf, hdir + "/Truth2planes", conf);
    htruth3planes_.init(hf, hdir + "/Truth3planes", conf);
    htruth4planes_.init(hf, hdir + "/Truth4planes", conf);
  }
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
  if(doMCTruth_) {
    htruthAfterNoHits_.fill(evt);
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
  if(doMCTruth_) {
    htruthAfterBeamVeto_.fill(evt);
  }

  for(int i=1; i<=4; ++i) {
    if(!protonGlobalClusters[i].empty()) {
      hHitPCsAterBeamVeto_->Fill(i);
    }
  }

  //----------------------------------------------------------------
  // Deal with multiple after-trigger time windows
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
  if(doMCTruth_) {
    htruthAfterMutliwinTime_.fill(evt);
  }

  //----------------------------------------------------------------
  // Run the track based analysis, and write out the "common" sample
  // before Z containment and window time cuts

  muCapTrkHF_.process(evt, muStopUV, protonWindow);
  if(commonSkimOutFile_) {
    commonSkimOutFile_<<evt.nrun<<" "<<evt.nevt<<std::endl;
  }

  //----------------------------------------------------------------
  // Trig time is 0, dt from that rather than from less precise trigWin time
  hWindowTimeBefore_->Fill(protonWindow.tstart);
  if( (protonWindow.tstart < cutWinTimeMin_) || (cutWinTimeMax_ < protonWindow.tstart)) {
    return CUT_WINTIME;
  }
  hWindowTimeAfter_->Fill(protonWindow.tstart);

  if(doMCTruth_) {
    htruthAfterWinTime_.fill(evt);
  }

  //----------------------------------------------------------------
  // UVAnalysis: DIOs for normalization
  const int iDIO = uvan_.process(evt, muStopUV);
  if(iDIO >= 0) {
    hwidthPCDIO_.fill(protonWindow.pcHits);
    hwidthDCDIO_.fill(protonWindow.dcHits);
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
  if(doMCTruth_) {
    htruthContained_.fill(evt);
  }

  //----------------------------------------------------------------
  const WireClusterCollection& pc7clusters = protonGlobalClusters[29];
  hNumPC7Clusters_->Fill(pc7clusters.size());
  hNumPC7WiresVsClusters_->Fill(pc7clusters.size(), MuCapUtilities::numWires(pc7clusters));
  if(pc7clusters.empty()) {
    return CUT_PC7_HIT;
  }

  //----------------------------------------------------------------
  double distanceToMuStop = std::numeric_limits<double>::max();
  for(WireClusterCollection::const_iterator ic = pc7clusters.begin(); ic != pc7clusters.end(); ++ic) {
    WirePlane::Measurement m7 = geom_->pc(7).measurement(ic->centralCell());

    const double offset = m7.coordinate - (m7.dir == WirePlane::U ? muStopUV.x() : muStopUV.y());
    //std::cout<<"m7.dir = "<<m7.dir<<", coord = "<<m7.coordinate<<", muStopUV = "<<muStopUV<<", offset = "<<offset<<std::endl;
    const double current = std::max(0., std::abs(offset) - geom_->pc(7).wireSpacing() * ic->numCells()/2.);
    distanceToMuStop = std::min(distanceToMuStop, current);
  }
  hPC7DistanceToMuStop_->Fill(distanceToMuStop);
  if(distanceToMuStop > cutPC7MaxDistanceToMuStop_) {
    return CUT_PC7_COORDINATE;
  }

  hRangeAfterPC7Cuts_.fill(protonDnClusters);

  //----------------------------------------------------------------
  // CUTS_LOOSE_PROTONS are passed at this point

  const PlaneRange gr = findPlaneRange(protonGlobalClusters);
  hLastPlaneLoose_->Fill(gr.max());

  hwidthPCLooseProtons_.fill(protonWindow.pcHits);
  hwidthDCLooseProtons_.fill(protonWindow.dcHits);
  hcLooseProtons_.fill(protonGlobalClusters);
  hhsLooseProtons_.fill(protonGlobalClusters);

  hfLoose_.fill(protonGlobalClusters, evt);

  if(doMCTruth_) {
    htruthLoose_.fill(evt);
    if(gr.max() >= 30) {
      htruth2planes_.fill(evt);
    }
    if(gr.max() >= 31) {
      htruth3planes_.fill(evt);
    }
    if(gr.max() >= 32) {
      htruth4planes_.fill(evt);
    }
  }

  return CUTS_LOOSE_PROTONS;
}
