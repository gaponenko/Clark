// Andrei Gaponenko, 2013

#include "MuCapture.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <utility>

#include "TH1.h"
#include "TH2.h"
#include "Math/Point2D.h"

#include "TimeWindow.h"
#include "PlaneRange.h"
#include "EventList.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//================================================================
bool MuCapture::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) {
  if (not Conf.read<bool>("MuCapture/Do")) {
    Log->info( "MuCapture code will not be run");
    return false ;
  }
  Log->info( "Register MuCapture module");

  doMCTruth_ = Conf.read<bool>("TruthBank/Do");
  gEventList = EventList(Conf.read<std::string>("MuCapture/debugEventList"));

  //       --------- Parameters initialization ---------          //
  doDefaultTWIST_ = Conf.read<bool>("MuCapture/doDefaultTWIST");
  cutMinTDCWidthPC_ = Conf.read<double>("MuCapture/cutMinTDCWidthPC");
  cutMinTDCWidthDC_ = Conf.read<double>("MuCapture/cutMinTDCWidthDC");
  winPCPreTrigSeparation_ = Conf.read<double>("MuCapture/winPCPreTrigSeparation");
  maxUnassignedDCHits_ = Conf.read<int>("MuCapture/maxUnassignedDCHits");

  muStopRMax_ = Conf.read<double>("MuCapture/muStopRMax");

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

  hNumAfterTrigWindows_ = H.DefineTH1D("MuCapture", "numAfterTrigPCWindows", "numAfterTrigPCWindows", 10, -0.5, 9.5);

  hWinDCUnassignedCount_ = H.DefineTH1D("MuCapture", "winDCUnassignedCount", "Count of unassigned DC hits", 101, -0.5, 100.5);

  // Make the bin size half a cell
  hMuStopUVCell_ = H.DefineTH2D("MuCapture", "MuStopUVCell", "Muon stop V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);
  hMuStopUVPos_ = H.DefineTH2D("MuCapture", "MuStopUVPos", "Muon stop V vs U position (cm)", 201, -10.05, +10.05, 201, -10.05, +10.05 );
  hMuStopRadius_ = H.DefineTH1D("MuCapture", "MuStopRadius", "Muon stop R (cm)", 80, 0., 8.);

  //----------------------------------------------------------------
  hwidthPCall_.init("MuCapture/pcWidthAll", "pcwidth", 12, H, Conf);
  hwidthDCall_.init("MuCapture/dcWidthAll", "dcwidth", 44, H, Conf);

  //----------------------------------------------------------------
  hOccupancyPCAll_.init("MuCapture", "hitMapPCAll", 12, 160, H, Conf);
  hOccupancyDCAll_.init("MuCapture", "hitMapDCAll", 44, 80, H, Conf);

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

  // Fill the "all hits" histos
  hwidthPCall_.fill(evt.pc_hits());
  hwidthDCall_.fill(evt.dc_hits());

  hOccupancyPCAll_.fill(evt.pc_hits());
  hOccupancyDCAll_.fill(evt.dc_hits());

  if(doMCTruth_) {
    hTruthAll_.fill(evt);
  }

  //----------------------------------------------------------------
  // Sort PC hits into time windows

  const TDCHitWPPtrCollection filteredPChits = selectHits(evt.pc_hits(), cutMinTDCWidthPC_);
  if(filteredPChits.empty()) {
    return CUT_NOPCHITS;
  }

  TimeWindowingResults wres;
  pcWindowing_.assignPCHits(filteredPChits, &wres);

  if(wres.iTrigWin == -1u) {
    return CUT_NOTRIGWIN;
  }

  if(wres.iTrigWin > 0) {
    // Trig time is 0, dt from that rather than from less precise trigWin time
    const double dt = wres.windows[wres.iTrigWin - 1].tstart;
    hPCPreTrigSeparation_->Fill(dt);
    if(std::abs(dt) < winPCPreTrigSeparation_) {
      return CUT_PCWIN_TRIGSEPPAST;
    }
  }

  const TimeWindow& trigWin = wres.windows[wres.iTrigWin];
  if(trigWin.stream != TimeWindow::UPSTREAM) {
    return CUT_TRIGPCWIN_TYPE;
  }

  const ClustersByPlane muonPCClusters = constructPlaneClusters(12, trigWin.pcHits);
  const PlaneRange trigPCRange = findPlaneRange(muonPCClusters);
  if(!trigPCRange.noGaps) {
    return CUT_TRIGPCWIN_GAPS;
  }

  if((trigPCRange.min != 1) || (trigPCRange.max != 6)) {
    return CUT_TRIGPCWIN_RANGE;
  }

  //----------------
  // Process DC hits
  dcWindowing_.assignDCHits(selectHits(evt.dc_hits(), cutMinTDCWidthDC_), &wres);
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
  if(!muonRange.noGaps) {
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
  hNumAfterTrigWindows_->Fill(int(wres.windows.size() - wres.iTrigWin - 1));
  if(wres.iTrigWin + 2 != wres.windows.size()) {
    return CUT_WIN_NUMAFTERTRIG;
  }

  //----------------------------------------------------------------
  const unsigned iProtonWin = 1 + wres.iTrigWin;
  const TimeWindow& protonWin = wres.windows[iProtonWin];
  const ClustersByPlane protonPCClusters = constructPlaneClusters(12, protonWin.pcHits);
  const ClustersByPlane protonDCClusters = constructPlaneClusters(44, protonWin.dcHits);
  const ClustersByPlane protonGlobalClusters = globalPlaneClusters(protonPCClusters, protonDCClusters);

  protonWindow_.process(muStop, wres, evt);
  anDnLate_.process(evt, wres, iProtonWin, muStop, protonGlobalClusters);

  //----------------------------------------------------------------
  return CUTS_ACCEPTED;
}

//================================================================
TDCHitWPPtrCollection MuCapture::selectHits(const TDCHitWPCollection& hits, double minWidthCut) {
  TDCHitWPPtrCollection res;
  for(unsigned i=0; i<hits.size(); ++i) {
    if(hits[i].width() > minWidthCut) {
      res.push_back(TDCHitWPPtr(hits, i));
    }
  }
  return res;
}

//================================================================
