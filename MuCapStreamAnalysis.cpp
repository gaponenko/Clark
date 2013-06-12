// Andrei Gaponenko, 2013

#include "MuCapStreamAnalysis.h"

#include <algorithm>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "PlaneRange.h"
#include "EventList.h"

//================================================================
void MuCapStreamAnalysis::init(HistogramFactory &hf, const std::string& hdir,
                               const DetectorGeo& geom, const ConfigFile& conf,
                               TimeWindow::StreamType cutWinStream, double cutWinTimeMin)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");
  cutStream_ = cutWinStream;
  cutWinTimeMin_ = cutWinTimeMin;
  cutWinTimeMax_ = conf.read<double>("MuCapture/cutWinTimeMax");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "protonCuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D(hdir, "protonCuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hWindowTime_ = hf.DefineTH1D(hdir, "windowTime", "Proton window start time", 1000, 0., 10000.);

  //hNumClusters_ = hf.DefineTH2D(hdir, "numClustersVsPlane", "Num cluster vs plane", 56, 0.5, 56.5,  11, -0.5, 10.5);
  //hNumClusters_->SetOption("colz");

  //----------------------------------------------------------------
  uvan_.init(hdir+"/UVAnalysis", hf, conf, cutStream_);

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
                                  unsigned iProtonWin,
                                  const ROOT::Math::XYPoint& muStopUV,
                                  const ClustersByPlane& protonGlobalClusters)
{
  EventCutNumber c = analyze(evt, wres, iProtonWin, muStopUV, protonGlobalClusters);
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
        unsigned iProtonWin,
        const ROOT::Math::XYPoint& muStopUV,
        const ClustersByPlane& protonGlobalClusters)
{
  const TimeWindow& protonWindow = wres.windows[iProtonWin];

  // Trig time is 0, dt from that rather than from less precise trigWin time
  if( (protonWindow.tstart < cutWinTimeMin_) || (cutWinTimeMax_ < protonWindow.tstart)) {
    return CUT_WINTIME;
  }

  const unsigned numDIO = uvan_.process(evt, protonWindow.tstart, protonGlobalClusters, muStopUV);
  if(numDIO) {
    hwidthPCTightDIO_.fill(protonWindow.pcHits);
    hwidthDCTightDIO_.fill(protonWindow.dcHits);
    hcdio_.fill(protonGlobalClusters);
  }


//  if(protonWindow.stream != cutStream_) {
//    return CUT_STREAM;
//  }
//
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
