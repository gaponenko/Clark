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
    double t = coll.front().hits().front()->time();
    for(WireClusterCollection::const_iterator c = coll.begin(); c!=coll.end(); ++c) {
      for(TDCHitWPPtrCollection::const_iterator h = c->hits().begin(); h != c->hits().end(); ++h) {
        if(t > (*h)->time()) {
          t = (*h)->time();
        }
      }
    }
    return t;
  }
}

//================================================================
void MuCapProtonWindow::init(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile& conf,
                             TimeWindow::StreamType cutWinStream, double cutAfterTrigTimeSep)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  cutStream_ = cutWinStream;
  if(cutStream_ != TimeWindow::DOWNSTREAM) {
    throw std::runtime_error("MuCapProtonWindow: requested cutWinStream value is not implemented");
  }

  cutAfterTrigTimeSep_ = cutAfterTrigTimeSep;

  cutMaxPlane_ = conf.read<double>("MuCapture/ProtonWindow/maxPlane");
  cutRextMax_ = conf.read<double>("MuCapture/ProtonWindow/RextMax");
  tightProtonsOutFileName_ = conf.read<std::string>("MuCapture/ProtonWindow/TightProtonsFileName", "");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D("MuCapture/ProtonWindow", "protonCuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = hf.DefineTH1D("MuCapture/ProtonWindow", "protonCuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hNumClusters_ = hf.DefineTH2D("MuCapture/ProtonWindow", "numClustersVsPlane", "Num cluster vs plane", 56, 0.5, 56.5,  11, -0.5, 10.5);
  hNumClusters_->SetOption("colz");

  hStartPos_ = hf.DefineTH2D("MuCapture/ProtonWindow", "startPos", "Proton start V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);
  hStartPos_->SetOption("colz");

  hStartOffset_ = hf.DefineTH2D("MuCapture/ProtonWindow", "startOffset", "Proton start offset V vs U (cell units)", 41, -5.125, 5.125,  41, -5.125, 5.125);
  hStartOffset_->SetOption("colz");

  hStartClusterSize_ = hf.DefineTH2D("MuCapture/ProtonWindow", "startClusterSize", "Proton start cluster V vs U width (cell units)", 25, -0.5, 24.5,  25, -0.5, 24.5);
  hStartClusterSize_->SetOption("colz");

  hLastPlane_ = hf.DefineTH1D("MuCapture/ProtonWindow", "lastPlane", "Proton last plane", 56, 0.5, 56.5);
  hProtonTime_ = hf.DefineTH1D("MuCapture/ProtonWindow", "protonTime", "Proton time", 1000, 0., 10000.);

  hwidthPCTightProtons_.init("MuCapture/ProtonWindow/pcWidthTightProtons", "pcpwidth", 12, hf, conf);
  hwidthDCTightProtons_.init("MuCapture/ProtonWindow/dcWidthTightProtons", "dcpwidth", 44, hf, conf);
  hwidthPCTightDIO_.init("MuCapture/ProtonWindow/pcWidthTightDIO", "pcpwidth", 12, hf, conf);
  hwidthDCTightDIO_.init("MuCapture/ProtonWindow/dcWidthTightDIO", "dcpwidth", 44, hf, conf);

  hcLooseProtons_.init(hf, "MuCapture/ProtonWindow/clLooseProtons", geom, conf);
  hcTightProtons_.init(hf, "MuCapture/ProtonWindow/clTightProtons", geom, conf);
  hcdio_.init(hf, "MuCapture/ProtonWindow/cldio", geom, conf);

  hCCRvsPlaneDIO_ = hf.DefineTH2D("MuCapture/ProtonWindow", "RvsPlaneDIO",
                                  "Extrpolated Rmax vs plane for DIO",
                                  56, 0.5, 56.5, 250, 0., 25.);

  hCCRvsPlaneDIO_->SetOption("colz");

  hCCRvsPlaneProtons_ = hf.DefineTH2D("MuCapture/ProtonWindow", "RvsPlaneProtons",
                                      "Extrpolated Rmax vs plane for Protons",
                                      56, 0.5, 56.5, 250, 0., 25.);

  hCCRvsPlaneProtons_->SetOption("colz");

  uvan_.init("MuCapture/ProtonWindow/UVAnalysis", hf, conf);
  hpw_.init(hf, "MuCapture/ProtonWindow/Final", geom, conf);
  rcheckDIO_.init(hf, "MuCapture/ProtonWindow/DIORCheck", geom, conf);
  rcheckProtonCandidates_.init(hf, "MuCapture/ProtonWindow/PCRCheck", geom, conf);

  if(doMCTruth_) {
    hrtruth_.init(hf, "MuCapture/ProtonWindow/RTruth", conf);
    htruthLoose_.init(hf, "MuCapture/ProtonWindow/TruthLoose", conf);
    htruthPC8_.init(hf, "MuCapture/ProtonWindow/TruthPC8", conf);
    htruthMinRange_.init(hf, "MuCapture/ProtonWindow/TruthMinRange", conf);
    htruthTight_.init(hf, "MuCapture/ProtonWindow/TruthTight", conf);

    hLastPlaneVsMCPstart_ = hf.DefineTH2D("MuCapture/ProtonWindow", "lastPlaneVsMCPstart",
                                          "Last plane vs MC pstart",
                                          200, 0., 200., 56, 0.5, 56.5);

    hLastPlaneVsMCPstart_->SetOption("colz");
  }

  if(!tightProtonsOutFileName_.empty()) {
    tightProtonsOutFile_.open(tightProtonsOutFileName_.c_str());
    if(!tightProtonsOutFile_) {
      throw std::runtime_error("Error opening output file "+tightProtonsOutFileName_);
    }
  }
}

//================================================================
void MuCapProtonWindow::process(const ROOT::Math::XYPoint& muStopUV,
                                const TimeWindowingResults& wres,
                                const EventClass& evt
                                )
{
  EventCutNumber c = analyze(muStopUV, wres, evt);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
}

//================================================================
MuCapProtonWindow::EventCutNumber MuCapProtonWindow::
analyze(const ROOT::Math::XYPoint& muStopUV,
        const TimeWindowingResults& wres,
        const EventClass& evt
        )
{
  const TimeWindow& protonWindow = wres.windows[1+wres.iTrigWin];

  // Trig time is 0, dt from that rather than from less precise trigWin time
  if(protonWindow.tstart < cutAfterTrigTimeSep_) {
    return CUT_TRIGSEP;
  }

  const TDCHitWPPtrCollection& protonPCHits = protonWindow.pcHits;
  const TDCHitWPPtrCollection& protonDCHits = protonWindow.dcHits;
  const ClustersByPlane clustersPC = constructPlaneClusters(12, protonPCHits);
  const ClustersByPlane clustersDC = constructPlaneClusters(44, protonDCHits);
  const ClustersByPlane global = globalPlaneClusters(clustersPC, clustersDC);
  const PlaneRange gr = findPlaneRange(global);

  if(protonWindow.stream != cutStream_) {
    return CUT_STREAM;
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

  const unsigned numDIO = uvan_.process(evt,  protonWindow.tstart, global, muStopUV);
  if(numDIO) {

    hwidthPCTightDIO_.fill(protonPCHits);
    hwidthDCTightDIO_.fill(protonDCHits);
    hcdio_.fill(global);

    for(int plane=32; plane <= std::min(cutMaxPlane_, gr.max); ++plane) {
      const double r = rcheckDIO_.rmax(plane, global);
      hCCRvsPlaneDIO_->Fill(plane, r);
    }
  }

  hLastPlane_->Fill(gr.max);
  if(doMCTruth_ && (evt.nmcvtx == 2)) {
    hLastPlaneVsMCPstart_->Fill(evt.mcvertex_ptot[0], gr.max);
  }

  if(gr.max > cutMaxPlane_) {
    return CUT_MAX_RANGE;
  }

  //----------------------------------------------------------------
  // CUTS_LOOSE_PROTONS are passed at this point

  hcLooseProtons_.fill(global);

  if(doMCTruth_) {
    htruthLoose_.fill(evt);
  }

  // If we use window start time here, longer proton tracks will get sharper timing
  // because of more hits.  Restrict the amount of data we use for proton time to
  // just PC7, so that all tracks get equally bad timing.
  hProtonTime_->Fill(getMinTime(clustersPC[7]));

  hpw_.fill(global, evt);

  //----------------------------------------------------------------
  if(gr.max < 30) {
    return CUT_NOPC8;
  }

  if(doMCTruth_) {
    htruthPC8_.fill(evt);
  }

  //----------------------------------------------------------------
  // the containment check requires at least 2U and 2V planes
  if(gr.max < 32) {
    return CUT_MIN_RANGE;
  }

  if(doMCTruth_) {
    htruthMinRange_.fill(evt);
  }

  const double rext = rcheckProtonCandidates_.rmax(gr.max, global);
  hCCRvsPlaneProtons_->Fill(gr.max, rext);
  if(doMCTruth_) {
    hrtruth_.fill(evt, gr.max, rext);
  }
  if(rext > cutRextMax_) {
    return CUT_REXT;
  }

  if(doMCTruth_) {
    htruthTight_.fill(evt);
  }

  if(tightProtonsOutFile_) {
    tightProtonsOutFile_<<evt.nrun<<" "<<evt.nevt<<std::endl;
  }

  hwidthPCTightProtons_.fill(protonPCHits);
  hwidthDCTightProtons_.fill(protonDCHits);
  hcTightProtons_.fill(global);

  return CUTS_TIGHT_PROTONS;
}
