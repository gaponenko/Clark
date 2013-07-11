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
  cutTrigPCWinGapsEnabled_ = Conf.read<bool>("MuCapture/cutTrigPCWinGapsEnabled");
  cutTrigPCWinStartPlane_ = Conf.read<int>("MuCapture/cutTrigPCWinStartPlane");
  maxUnassignedDCHits_ = Conf.read<int>("MuCapture/maxUnassignedDCHits");

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
  muStopOutFileName_ = Conf.read<std::string>("MuCapture/muStopOutFileName", "");
  if(!muStopOutFileName_.empty()) {
    muStopOutFile_.open(muStopOutFileName_.c_str());
    if(!muStopOutFile_) {
      throw std::runtime_error("Error opening output file "+muStopOutFileName_);
    }
  }

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
  anDnLate_.init(H, "MuCapture/dnLate", *E.geo, Conf, TimeWindow::DOWNSTREAM, 1050./*FIXME*/);

  h_cuts_r = H.DefineTH1D("MuCapture", "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = H.DefineTH1D("MuCapture", "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hPCPreTrigSeparation_ = H.DefineTH1D("MuCapture", "pcPreTrigSeparation", "PC win pre-trigger separation", 600, -6000., 0.);

  hTrigPCWinStartPlane_ = H.DefineTH1D("MuCapture", "trigPCWinStartPlane", "Trig win start PC plane", 12, 0.5, 12.5);

  hTrigPCWinGaps_ = H.DefineTH2D("MuCapture", "trigPCWinGaps", "Trig win PC gap end vs start", 13, -0.5, 12.5, 13, -0.5, 12.5);
  hTrigPCWinGaps_->SetOption("colz");

  hTrigPCWinMissingPlanes_ = H.DefineTH1D("MuCapture", "trigPCWinMissingPlanes", "Trig win missing PC planes", 12, 0.5, 12.5);

  hWinDCUnassignedCount_ = H.DefineTH1D("MuCapture", "winDCUnassignedCount", "Count of unassigned DC hits", 101, -0.5, 100.5);

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
  winTimeBeforeTrigPCWinType_.init("MuCapture/winTime", "beforeTrigPCWinType", H, Conf);
  winTimeBeforeTrigPCWinGaps_.init("MuCapture/winTime", "beforeTrigPCWinGaps", H, Conf);

  winTimeBeforeTrigPCWinEndPlane_.init("MuCapture/winTime", "beforeTrigPCEndPlane", H, Conf);
  winTimeBeforeTrigPCWinStartPlane_.init("MuCapture/winTime", "beforeTrigPCStartPlane", H, Conf);

  winTimeBeforeTrigDCWinType_.init("MuCapture/winTime", "beforeTrigDCWinType", H, Conf);
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
    hTruthAll_.init(H, "MuCapture/MCTruthAll", Conf);
  }

  //----------------------------------------------------------------

  return true;
}

//================================================================
bool MuCapture::Process(EventClass &evt, HistogramFactory &hist) {

  EventCutNumber c = analyze(evt, hist);
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
  const TDCHitWPPtrCollection allPCHits = allPCHitsBuf.get();

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
  const TDCHitWPPtrCollection allDCHits = allDCHitsBuf.get();

  TDCHitPreprocessing::Hits filteredDCHitsBuf(evt.dc_hits(), *dcHitProcessor_);
  const TDCHitWPPtrCollection filteredDCHits = filteredDCHitsBuf.get();

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

  winTimeBeforeTrigPCWinType_.fill(wres);

  const TimeWindow& trigWin = wres.windows[wres.iTrigWin];
  if(trigWin.stream != TimeWindow::UPSTREAM) {
    return CUT_TRIGPCWIN_TYPE;
  }

  const ClustersByPlane muonPCClusters = constructPlaneClusters(12, trigWin.pcHits);
  const PlaneRange trigPCRange = findPlaneRange(muonPCClusters);

  winTimeBeforeTrigPCWinEndPlane_.fill(wres);
  if(trigPCRange.max() != 6) {
    return CUT_TRIGPCWIN_END_PLANE;
  }

  winTimeBeforeTrigPCWinStartPlane_.fill(wres);
  hTrigPCWinStartPlane_->Fill(trigPCRange.min());
  if(cutTrigPCWinStartPlane_ < trigPCRange.min()) {
    return CUT_TRIGPCWIN_START_PLANE;
  }

  winTimeBeforeTrigPCWinGaps_.fill(wres);
  for(int i = 0; i + 1 < trigPCRange.segments().size(); ++i) {
    hTrigPCWinGaps_->Fill(trigPCRange.segments()[i].max, trigPCRange.segments()[i+1].min);
  }
  for(int iplane = 1; iplane <= 6; ++iplane) {
    if(muonPCClusters[iplane].empty()) {
      hTrigPCWinMissingPlanes_->Fill(iplane);
    }
  }
  if(cutTrigPCWinGapsEnabled_ && !trigPCRange.noGaps()) {
    return CUT_TRIGPCWIN_GAPS;
  }

  winTimeBeforeTrigDCWinType_.fill(wres);

  //----------------
  // Process DC hits
  dcWindowing_.assignDCHits(filteredDCHits, &wres);
  if(trigWin.stream != TimeWindow::UPSTREAM) {
    return CUT_TRIGDCWIN_TYPE;
  }

  hWinDCUnassignedCount_->Fill(wres.unassignedDCHits.size());
  if(wres.unassignedDCHits.size() > maxUnassignedDCHits_) {
    return CUT_UNASSIGNEDDCHITS;
  }

  //----------------------------------------------------------------
  const ClustersByPlane muonDCClusters = constructPlaneClusters(44, trigWin.dcHits);
  const ClustersByPlane muonGlobalClusters = globalPlaneClusters(muonPCClusters, muonDCClusters);

  if(gEventList.requested(evt)) {
    std::cout<<__func__<<": run "<<evt.nrun<<" event "<<evt.nevt
             <<": muonGlobalClusters = "<<muonGlobalClusters
             <<std::endl;
  }

  const PlaneRange muonRange = findPlaneRange(muonGlobalClusters);
  if(!muonRange.noGaps()) {
    return CUT_MU_RANGE_GAPS;
  }

  //----------------------------------------------------------------
  // CUT_MUSTOP_UV
  if((muonPCClusters[5].size() != 1) || (muonPCClusters[6].size() != 1)) {
    return CUT_MUSTOP_SINGLECLUSTER;
  }

  const double pc5wire = muonPCClusters[5].front().centralCell();
  const double pc6wire = muonPCClusters[6].front().centralCell();

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

  //----------------------------------------------------------------
  if(1 != pactCut_.quadrant(muonPCClusters[5].front(), muonPCClusters[6].front())) {
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

  if(muStopOutFile_) {
    muStopOutFile_<<evt.nrun<<" "<<evt.nevt<<std::endl;
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
