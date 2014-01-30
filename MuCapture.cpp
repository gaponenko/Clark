// Andrei Gaponenko, 2013

#include "MuCapture.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <utility>
#include <stdexcept>
#include <limits>
#include <iomanip>

#include "TH1.h"
#include "TH2.h"
#include "Math/Point2D.h"

#include "TimeWindow.h"
#include "PlaneRange.h"
#include "EventList.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//================================================================
namespace {
  struct WinHitCounter {
    unsigned numPCHits;
    unsigned numDCHits;
    WinHitCounter() : numPCHits(0), numDCHits(0) {}
    WinHitCounter(unsigned p, unsigned d) : numPCHits(p), numDCHits(d) {}

    WinHitCounter operator+(const TimeWindow& w) const {
      return WinHitCounter(numPCHits + w.pcHits.size(),
                           numDCHits + w.dcHits.size());
    }
  };
}

//================================================================
bool MuCapture::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) {
  if (not Conf.read<bool>("MuCapture/Do")) {
    Log->info( "MuCapture code will not be run");
    return false ;
  }
  Log->info( "Register MuCapture module");

  //       --------- Parameters initialization ---------          //
  doDefaultTWIST_ = Conf.read<bool>("MuCapture/doDefaultTWIST");
  winPCPreTrigSeparation_ = Conf.read<double>("MuCapture/winPCPreTrigSeparation");
  maxUnassignedDCHits_ = Conf.read<int>("MuCapture/maxUnassignedDCHits");
  cutMuonFirstPlane_ = Conf.read<int>("MuCapture/cutMuonFirstPlane");
  cutMuonRangeGapsEnabled_ = Conf.read<bool>("MuCapture/cutMuonRangeGapsEnabled");

  muStopRMax_ = Conf.read<double>("MuCapture/muStopRMax");

  fillXtalkPC_ = Conf.read<bool>("MuCapture/fillXtalkPC");
  fillXtalkDC_ = Conf.read<bool>("MuCapture/fillXtalkDC");
  doMCTruth_ = Conf.read<bool>("TruthBank/Do");
  inputNumberList_ = EventList(Conf.read<std::string>("MuCapture/inputEventNumberFile"));
  gEventList = EventList(Conf.read<std::string>("MuCapture/debugEventList"));

  //----------------------------------------------------------------
  pcHitProcessor_ = makeTDCHitPreprocessor(WirePlane::PC, H, *E.geo, Conf);
  dcHitProcessor_ = makeTDCHitPreprocessor(WirePlane::DC, H, *E.geo, Conf);

  //----------------------------------------------------------------
  if(doMCTruth_) {
    anDnLateResponse_.Setup(20, 0., 200, 20, 0., 200.);
    H.Store(&anDnLateResponse_, "anDnLateResponse", "MuCapture");

    hTruthMomentum_ = H.DefineTH1D("MuCapture/LateResponse", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",20, 0., 200);
    hTruthMomentumReco_ = H.DefineTH1D("MuCapture/LateResponse", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",20, 0., 200);
    hMeasVsTruthMomentum_ = H.DefineTH2D("MuCapture/LateResponse", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",20, 0., 200,20, 0., 200);
    hTruthMomentumNotReco_ = H.DefineTH1D("MuCapture/LateResponse", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",20, 0., 200);
    hMeasuredMomentum_ = H.DefineTH1D("MuCapture/LateResponse", "MCMeasuredMomentum", "Measured momentum used in response function;Momentum [MeV/c]",20, 0., 200);
  hMuStopRadius_ = H.DefineTH1D("MuCapture", "MuStopRadius", "Muon stop R (cm)", 80, 0., 8.);

  }

  //----------------------------------------------------------------
  uvOutFileName_ = Conf.read<std::string>("MuCapture/uvOutFileName", "");
  if(!uvOutFileName_.empty()) {
    uvOutFile_.open(uvOutFileName_.c_str());
    if(!uvOutFile_) {
      throw std::runtime_error("Error opening output file "+uvOutFileName_);
    }
  }

  //       --------- Histograms initialization ---------          //
  pcWindowing_.init(H, "MuCapture/WindowingPC", *E.geo, Conf);
  dcWindowing_.init(H, "MuCapture/WindowingDC", *E.geo, Conf);

  pactCut_.init(H, Conf);
  protonWindow_.init(H, *E.geo, Conf, TimeWindow::DOWNSTREAM, 1050./*FIXME*/);

  //anUpLate_.init(H, *E.geo, Conf, TimeWindow::DOWNSTREAM, 1050./*FIXME*/);
  anDnLate_.init(H, "MuCapture/dnLate", *E.geo, Conf, &anDnLateRes_, TimeWindow::DOWNSTREAM, 1050./*FIXME*/);

  h_cuts_r = H.DefineTH1D("MuCapture", "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = H.DefineTH1D("MuCapture", "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hPCPreTrigSeparation_ = H.DefineTH1D("MuCapture", "pcPreTrigSeparation", "PC win pre-trigger separation", 600, -6000., 0.);

  hWinDCUnassignedCount_ = H.DefineTH1D("MuCapture", "winDCUnassignedCount", "Count of unassigned DC hits", 101, -0.5, 100.5);

  hMuonFirstPlane_ = H.DefineTH1D("MuCapture", "muonFirstPlane", "Muon first plane", 56, 0.5, 56.5);
  hMuonLastPlaneBeforeGaps_ = H.DefineTH1D("MuCapture", "muonLastPlaneBeforeGaps", "Muon last plane, before gaps", 56, 0.5, 56.5);
  hMuonLastPlaneAfterGaps_ = H.DefineTH1D("MuCapture", "muonLastPlaneAfterGaps", "Muon last plane", 56, 0.5, 56.5);

  hMuonRangeGaps_ = H.DefineTH2D("MuCapture", "muonRangeGaps", "Muon win gap end vs start", 57, -0.5, 56.5, 57, -0.5, 56.5);
  hMuonRangeGaps_->SetOption("colz");

  hMuonMissingPlanes_ = H.DefineTH1D("MuCapture", "muonMissingPlanes", "Muon range missing planes", 56, 0.5, 56.5);

  // Make the bin size half a cell
  hMuStopUVCell_ = H.DefineTH2D("MuCapture", "MuStopUVCell", "Muon stop V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);
  hMuStopUVPos_ = H.DefineTH2D("MuCapture", "MuStopUVPos", "Muon stop V vs U position (cm)", 201, -10.05, +10.05, 201, -10.05, +10.05 );
  hMuStopRadius_ = H.DefineTH1D("MuCapture", "MuStopRadius", "Muon stop R (cm)", 80, 0., 8.);

  //----------------------------------------------------------------
  hwidthPCall_.init("MuCapture/pcWidthAll", "pcwidth", 12, H, Conf);
  hwidthPCfiltered_.init("MuCapture/pcWidthFiltered", "pcwidth", 12, H, Conf);
  hwidthDCall_.init("MuCapture/dcWidthAll", "dcwidth", 44, H, Conf);
  hwidthDCfiltered_.init("MuCapture/dcWidthFiltered", "dcwidth", 44, H, Conf);

  if(fillXtalkPC_) {
    hAfterPulsingPCAll_.init("MuCapture/afterPulsingPCAll", 12, 160, H, Conf);
    hAfterPulsingPCFiltered_.init("MuCapture/afterPulsingPCFiltered", 12, 160, H, Conf);
    hXtalkSameWirePC_.init("MuCapture/xtalkSameWirePC", 12, 0, 0, H, Conf);
    hXtalk1PC_.init("MuCapture/xtalk1PC", 12, 1, 1, H, Conf);
    hXtalkPlanePC_.init("MuCapture/xtalkPlanePC", 12, 1, 999, H, Conf);
  }

  if(fillXtalkDC_) {
    hAfterPulsingDCAll_.init("MuCapture/afterPulsingDCAll", 44, 80, H, Conf);
    hAfterPulsingDCFiltered_.init("MuCapture/afterPulsingDCFiltered", 44, 80, H, Conf);
    hXtalkSameWireDC_.init("MuCapture/xtalkSameWireDC", 44, 0, 0, H, Conf);
    hXtalk1DC_.init("MuCapture/xtalk1DC", 44, 1, 1, H, Conf);
    hXtalkPlaneDC_.init("MuCapture/xtalkPlaneDC", 44, 1, 999, H, Conf);
  }

  //----------------------------------------------------------------
  hOccupancyPCAll_.init("MuCapture", "hitMapPCAll", 12, 160, H, Conf);
  hOccupancyDCAll_.init("MuCapture", "hitMapDCAll", 44, 80, H, Conf);
  haccidentalsTrig_.init("MuCapture/AccidentalsTrig", H, Conf);
  haccidentalsStop_.init("MuCapture/AccidentalsStop", H, Conf);

  winTimeBeforeNoTrigWin_.init("MuCapture/winTime", "beforeNoTrigWin", H, Conf);
  winTimeMuStop_.init("MuCapture/winTime", "muStop", H, Conf);

  dioUp_.init("MuCapture/DIOUp", H, Conf,
              TimeWindow::UPSTREAM,
              Conf.read<double>("MuCapture/DIOUp/cutMinTime"));

  dioDn_.init("MuCapture/DIODn", H, Conf,
              TimeWindow::DOWNSTREAM,
              Conf.read<double>("MuCapture/DIODn/cutMinTime"));

  hdriftPCAll_.init("MuCapture/driftTimePCAll", H, 12, 1000./*ns*/,
                    Conf.read<double>("MuCapture/HistDriftTime/cutEffTrackHitDtPC"),
                    Conf);
  hdriftPCFiltered_.init("MuCapture/driftTimePCFiltered", H, 12, 1000./*ns*/,
                         Conf.read<double>("MuCapture/HistDriftTime/cutEffTrackHitDtPC"),
                         Conf);

  hdriftDCAll_.init("MuCapture/driftTimeDCAll", H, 44, 5000./*ns*/,
                    Conf.read<double>("MuCapture/HistDriftTime/cutEffTrackHitDtDC"),
                    Conf);
  hdriftDCFiltered_.init("MuCapture/driftTimeDCFiltered", H, 44, 5000./*ns*/,
                         Conf.read<double>("MuCapture/HistDriftTime/cutEffTrackHitDtDC"),
                         Conf);

  if(doMCTruth_) {
    hmuStopTruthAll_.init(H, "MuCapture/MuStopTruthAll", *E.geo, Conf);
    hmuStopTruthAfterGaps_.init(H, "MuCapture/MuStopTruthAfterGaps", *E.geo, Conf);

    hTruthAll_.init(H, "MuCapture/MCTruthAll", Conf);
    hTruthAfterPreTrigHits_.init(H, "MuCapture/MCTruthAfterPreTrigHits", Conf);
    hTruthAfterUnassignedDCHits_.init(H, "MuCapture/MCTruthAfterUnassignedDCHits", Conf);
    hTruthAfterMuRangeGaps_.init(H, "MuCapture/MCTruthAfterMuRangeGapss", Conf);
    hTruthAfterMuLastPlane_.init(H, "MuCapture/MCTruthAfterMuLastPlane", Conf);
    hTruthAfterMuSingleCluster_.init(H, "MuCapture/MCTruthAfterMuSingleCluster", Conf);
    hTruthAfterMuStopUV_.init(H, "MuCapture/MCTruthAfterMuStopUV", Conf);
    hTruthMuStop_.init(H, "MuCapture/MCTruthMuStop", Conf);
  }

  //----------------------------------------------------------------

  return true;
}

//================================================================
bool MuCapture::Process(EventClass &evt, HistogramFactory &hist) {

  anDnLateRes_.accepted = false;
  EventCutNumber c = analyze(evt, hist);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  if(doMCTruth_) {
    if(evt.iCaptureMcVtxStart != -1) {  // Signal event. Fill the unfolding matrix.
      const double p_true = evt.mcvertex_ptot[evt.iCaptureMcVtxStart];
      if(anDnLateRes_.accepted) {
        anDnLateResponse_.Fill(anDnLateRes_.momentum, p_true);
        hTruthMomentum_->Fill(p_true);
        hMeasuredMomentum_->Fill(anDnLateRes_.momentum);
        hTruthMomentumReco_->Fill(p_true);
        hMeasVsTruthMomentum_->Fill(p_true,anDnLateRes_.momentum);
      }
      else {
        anDnLateResponse_.Miss(p_true);
        hTruthMomentum_->Fill(p_true);
        hTruthMomentumNotReco_->Fill(p_true);
      }
    }
  }

  if(gEventList.requested(evt)) {
    const int cbin = h_cuts_p->GetXaxis()->FindFixBin(c);
    std::cout<<__func__<<": run "<<evt.nrun<<" event "<<evt.nevt
             <<" status "<<c<<": "<<h_cuts_p->GetXaxis()->GetBinLabel(cbin)
             <<std::endl;
  }

  return doDefaultTWIST_;
}

//================================================================
MuCapture::EventCutNumber MuCapture::analyze(EventClass &evt, HistogramFactory &hist) {

  // If input event numbers are defined we want to ignore everything else.
  // Thus this cut is before filling the "all" histograms.
  if(!inputNumberList_.empty() &&  !inputNumberList_.requested(evt)) {
    return CUT_EVENT_NUMBER;
  }

  TDCHitPreprocessing::PassThrough hitpass;
  TDCHitPreprocessing::Hits allPCHitsBuf(evt.pc_hits(), hitpass);
  const TDCHitWPPtrCollection& allPCHits = allPCHitsBuf.get();

  TDCHitPreprocessing::Hits filteredPCHitsBuf(evt.pc_hits(), *pcHitProcessor_);
  const TDCHitWPPtrCollection& filteredPCHits = filteredPCHitsBuf.get();

  hwidthPCall_.fill(allPCHits);
  hwidthPCfiltered_.fill(filteredPCHits);

  if(fillXtalkPC_) {
    hAfterPulsingPCAll_.fill(allPCHits);
    hAfterPulsingPCFiltered_.fill(filteredPCHits);
    hXtalkSameWirePC_.fill(allPCHits);
    hXtalk1PC_.fill(allPCHits);
    hXtalkPlanePC_.fill(allPCHits);
  }

  TDCHitPreprocessing::Hits allDCHitsBuf(evt.dc_hits(), hitpass);
  const TDCHitWPPtrCollection& allDCHits = allDCHitsBuf.get();

  TDCHitPreprocessing::Hits filteredDCHitsBuf(evt.dc_hits(), *dcHitProcessor_);
  const TDCHitWPPtrCollection& filteredDCHits = filteredDCHitsBuf.get();

  hwidthDCall_.fill(allDCHits);
  hwidthDCfiltered_.fill(filteredDCHits);

  if(fillXtalkDC_) {
    hAfterPulsingDCAll_.fill(allDCHits);
    hAfterPulsingDCFiltered_.fill(filteredDCHits);
    hXtalkSameWireDC_.fill(allDCHits);
    hXtalk1DC_.fill(allDCHits);
    hXtalkPlaneDC_.fill(allDCHits);
  }

  hOccupancyPCAll_.fill(evt.pc_hits());
  hOccupancyDCAll_.fill(evt.dc_hits());

  if(doMCTruth_) {
    hmuStopTruthAll_.fill(evt);
    hTruthAll_.fill(evt);
  }

  //----------------------------------------------------------------
  // Sort PC hits into time windows

  if(filteredPCHits.empty()) {
    return CUT_NOPCHITS;
  }

  TimeWindowingResults wres;
  pcWindowing_.assignPCHits(filteredPCHits, &wres);
  assert(filteredPCHits.size() ==
         std::accumulate(wres.windows.begin(), wres.windows.end(), WinHitCounter()).numPCHits);

  winTimeBeforeNoTrigWin_.fill(wres);

  if(wres.iTrigWin == -1u) {
    return CUT_NOTRIGWIN;
  }

  haccidentalsTrig_.fill(wres);

  if(wres.iTrigWin > 0) {
    // Trig time is 0, dt from that rather than from less precise trigWin time
    const double dt = wres.windows[wres.iTrigWin - 1].tstart;
    hPCPreTrigSeparation_->Fill(dt);
    if(std::abs(dt) < winPCPreTrigSeparation_) {
      return CUT_PCWIN_TRIGSEPPAST;
    }
  }

  if(doMCTruth_) {
    hTruthAfterPreTrigHits_.fill(evt);
  }

  const TimeWindow& trigWin = wres.windows[wres.iTrigWin];

  //----------------
  // Process DC hits
  dcWindowing_.assignDCHits(filteredDCHits, &wres);

  hWinDCUnassignedCount_->Fill(wres.unassignedDCHits.size());
  if(wres.unassignedDCHits.size() > maxUnassignedDCHits_) {
    return CUT_UNASSIGNEDDCHITS;
  }

  if(doMCTruth_) {
    hTruthAfterUnassignedDCHits_.fill(evt);
  }

  //----------------
  const ClustersByPlane muonPCtmp = constructPlaneClusters(12, trigWin.pcHits);
  const ClustersByPlane muonDCtmp = constructPlaneClusters(44, trigWin.dcHits);
  const ClustersByPlane muonGlobalClusters = globalPlaneClusters(muonPCtmp, muonDCtmp);
  const PlaneRange muonRange = findPlaneRange(muonGlobalClusters);

  hMuonFirstPlane_->Fill(muonRange.min());
  if(cutMuonFirstPlane_ < muonRange.min()) {
    return CUT_MUON_FIRST_PLANE;
  }

  //----------------
  hMuonLastPlaneBeforeGaps_->Fill(muonRange.max());
  for(int i = 0; i + 1 < muonRange.segments().size(); ++i) {
    hMuonRangeGaps_->Fill(muonRange.segments()[i].max + 1, muonRange.segments()[i+1].min - 1);
  }
  for(int iplane = 1; iplane <= muonRange.max(); ++iplane) {
    if(muonGlobalClusters[iplane].empty()) {
      hMuonMissingPlanes_->Fill(iplane);
    }
  }
  if(cutMuonRangeGapsEnabled_ && !muonRange.noGaps()) {
    return CUT_MUON_RANGE_GAPS;
  }

  if(doMCTruth_) {
    hmuStopTruthAfterGaps_.fill(evt);
    hTruthAfterMuRangeGaps_.fill(evt);
  }

  //----------------
  hMuonLastPlaneAfterGaps_->Fill(muonRange.max());
  if(muonRange.max() != 28) {
    return CUT_MUON_LAST_PLANE;
  }

  if(doMCTruth_) {
    hTruthAfterMuLastPlane_.fill(evt);
  }

  //----------------------------------------------------------------
  if(gEventList.requested(evt)) {
    std::cout<<__func__<<": run "<<evt.nrun<<" event "<<evt.nevt
             <<": muonGlobalClusters = "<<muonGlobalClusters
             <<std::endl;
  }

  //----------------------------------------------------------------
  // CUT_MUSTOP_UV
  if((muonGlobalClusters[27].size() != 1) || (muonGlobalClusters[28].size() != 1)) {
    return CUT_MUSTOP_SINGLECLUSTER;
  }

  if(doMCTruth_) {
    hTruthAfterMuSingleCluster_.fill(evt);
  }

  const double pc5wire = muonGlobalClusters[27].front().centralCell();
  const double pc6wire = muonGlobalClusters[28].front().centralCell();

  // See dt_geo.00061 and twist-coordinate-system.uvplanes.pdf
  hMuStopUVCell_->Fill(pc6wire, pc5wire);

  ROOT::Math::XYPoint muStop =
    WirePlane::uv(evt.geo->pc(5).measurement(pc5wire),
                  evt.geo->pc(6).measurement(pc6wire));

  hMuStopUVPos_->Fill(muStop.x(), muStop.y());

  const double muStopRadius = muStop.R();
  hMuStopRadius_->Fill(muStopRadius);
  if(muStopRadius > muStopRMax_) {
    return CUT_MUSTOP_UV;
  }

  if(doMCTruth_) {
    hTruthAfterMuStopUV_.fill(evt);
  }

  //----------------------------------------------------------------
  if(1 != pactCut_.quadrant(muonGlobalClusters[27].front(), muonGlobalClusters[28].front())) {
    return CUT_MUSTOP_PACT;
  }

  //----------------------------------------------------------------
  // Here we have an accepted muon stop

  // Compute clusters for after-trigger windows
  std::vector<ClustersByPlane> afterTrigClusters;
  for(int iwin = 1 + wres.iTrigWin; iwin < wres.windows.size(); ++iwin) {
    const TimeWindow& win = wres.windows[iwin];
    const ClustersByPlane pcClusters = constructPlaneClusters(12, win.pcHits);
    const ClustersByPlane dcClusters = constructPlaneClusters(44, win.dcHits);
    afterTrigClusters.push_back(globalPlaneClusters(pcClusters, dcClusters));
  }

  // Call the subanalyses
  winTimeMuStop_.fill(wres);
  haccidentalsStop_.fill(wres);
  protonWindow_.process(muStop, wres, evt);
  anDnLate_.process(evt, wres, muStop, afterTrigClusters);

  int idio = dioUp_.process(evt, muStop);
  if(idio >= 0) {
    hdriftPCAll_.fill(evt, idio, allPCHits);
    hdriftPCFiltered_.fill(evt, idio, filteredPCHits);
    hdriftDCAll_.fill(evt, idio, allDCHits);
    hdriftDCFiltered_.fill(evt, idio, filteredDCHits);
    if(uvOutFile_) {
      uvOutFile_<<std::fixed<<std::showpos<<evt.hefit_u0[idio]<<"\t"<<evt.hefit_v0[idio]<<std::endl;
    }
  }
  idio = dioDn_.process(evt, muStop);
  if(idio >= 0) {
    hdriftPCAll_.fill(evt, idio, allPCHits);
    hdriftPCFiltered_.fill(evt, idio, filteredPCHits);
    hdriftDCAll_.fill(evt, idio, allDCHits);
    hdriftDCFiltered_.fill(evt, idio, filteredDCHits);
    if(uvOutFile_) {
      uvOutFile_<<std::fixed<<std::showpos<<evt.hefit_u0[idio]<<"\t"<<evt.hefit_v0[idio]<<std::endl;
    }
  }

  if(doMCTruth_) {
    hTruthMuStop_.fill(evt);
  }

  //----------------------------------------------------------------
  return CUTS_MUSTOP_ACCEPTED;
}

//================================================================
TDCHitPreprocessing::IProcessor *
MuCapture::makeTDCHitPreprocessor(WirePlane::DetType d,
                                  HistogramFactory& hf,
                                  const DetectorGeo& geom,
                                  ConfigFile& conf)
{
  const std::string proc = conf.read<std::string>("MuCapture/HitPreproc/"+WirePlane::detName(d)+"/processor");
  if(proc == "NarrowHitDiscarder") {
    return new TDCHitPreprocessing::NarrowHitDiscarder("MuCapture/HitPreproc", d, hf, geom, conf);
  }
  if(proc == "SameWireHitDiscarder") {
    return new TDCHitPreprocessing::SameWireHitDiscarder("MuCapture/HitPreproc", d, hf, geom, conf);
  }
  throw std::runtime_error("MuCapture::makeTDCHitPreprocessor(): unknown processor name \""+proc+"\"");
}

//================================================================
