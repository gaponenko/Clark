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

  const std::string hdir="MuCapture";

  //       --------- Parameters initialization ---------          //
  doDefaultTWIST_ = Conf.read<bool>("MuCapture/doDefaultTWIST");
  winPCPreTrigSeparation_ = Conf.read<double>("MuCapture/winPCPreTrigSeparation");
  cutMuonFirstPlane_ = Conf.read<int>("MuCapture/cutMuonFirstPlane");

  muStopRMax_ = Conf.read<double>("MuCapture/muStopRMax");

  cutBeamVetoMaxPCplanes_ = Conf.read<double>("MuCapture/cutBeamVetoMaxPCplanes");

  cutWinTimeMin_ = Conf.read<double>("MuCapture/cutWinTimeMin");
  cutWinTimeMax_ = Conf.read<double>("MuCapture/cutWinTimeMax");

  cutMultiwinNextdt_ = Conf.read<double>("MuCapture/cutMultiwinNextdt");

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
    anDnLateResponse_.Setup(25, 0., 250, 25, 0., 250.);
    H.Store(&anDnLateResponse_, "anDnLateResponse", hdir);

    hTruthMomentum_ = H.DefineTH1D(hdir+"/LateResponse", "MCTruthMomentum", "True momentum used in response function;Momentum [MeV/c]",25, 0., 250);
    hTruthMomentumReco_ = H.DefineTH1D(hdir+"/LateResponse", "MCTruthMomentumReco", "True momentum of reconstructed tracks;Momentum [MeV/c]",25, 0., 250);
    hMeasVsTruthMomentum_ = H.DefineTH2D(hdir+"/LateResponse", "MCMeasVsTruthMomentum", "Measured vs. true momentum used in response function;True momentum [MeV/c];Measured momentum [MeV/c]",25, 0., 250,25, 0., 250);
    hTruthMomentumNotReco_ = H.DefineTH1D(hdir+"/LateResponse", "MCTruthMomentumNotReco", "True momentum of tracks not reconstructed;Momentum [MeV/c]",25, 0., 250);
  }

  hMeasuredMomentum_ = H.DefineTH1D(hdir+"/LateResponse", "MeasuredMomentum", "Measured momentum spectrum;Momentum [MeV/c]",25, 0., 250);

  //----------------------------------------------------------------
  uvOutFileName_ = Conf.read<std::string>(hdir+"/uvOutFileName", "");
  if(!uvOutFileName_.empty()) {
    uvOutFile_.open(uvOutFileName_.c_str());
    if(!uvOutFile_) {
      throw std::runtime_error("Error opening output file "+uvOutFileName_);
    }
  }

  //       --------- Histograms initialization ---------          //
  pcWindowing_.init(H, hdir+"/WindowingPC", *E.geo, Conf);
  dcWindowing_.init(H, hdir+"/WindowingDC", *E.geo, Conf);

  pactCut_.init(H, Conf);
  protonWindow_.init(H, *E.geo, Conf, TimeWindow::DOWNSTREAM, 1050./*FIXME*/);
  anDnLate_.init(H, hdir+"/dnLate", *E.geo, Conf, &anDnLateRes_, TimeWindow::DOWNSTREAM, 1050./*FIXME*/);

  h_cuts_r = H.DefineTH1D(hdir, "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());
  h_cuts_r->SetOption("hist text");

  h_cuts_p = H.DefineTH1D(hdir, "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);
  h_cuts_p->SetOption("hist text");

  hPCPreTrigSeparation_ = H.DefineTH1D(hdir, "pcPreTrigSeparation", "PC win pre-trigger separation", 600, -6000., 0.);

  hMuonFirstPlane_ = H.DefineTH1D(hdir, "muonFirstPlane", "Muon first plane", 56, 0.5, 56.5);
  hMuonLastPlaneAfterGaps_ = H.DefineTH1D(hdir, "muonLastPlaneAfterGaps", "Muon last plane", 56, 0.5, 56.5);

  hStoppedMuonRangeGaps_ = H.DefineTH2D(hdir+"/AcceptedMuStop", "muStopRangeGaps", "Stopped muon range gap end vs start", 57, -0.5, 56.5, 57, -0.5, 56.5);
  hStoppedMuonRangeGaps_->SetOption("colz");

  hStoppedMuonMissingPlanes_ = H.DefineTH1D(hdir+"/AcceptedMuStop", "muStopMissingPlanes", "Stopped muon missing planes", 56, 0.5, 56.5);

  // Make the bin size half a cell
  hMuStopUVCell_ = H.DefineTH2D(hdir, "MuStopUVCell", "Muon stop V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);
  hMuStopUVPos_ = H.DefineTH2D(hdir, "MuStopUVPos", "Muon stop V vs U position (cm)", 201, -10.05, +10.05, 201, -10.05, +10.05 );
  hMuStopRadius_ = H.DefineTH1D(hdir, "MuStopRadius", "Muon stop R (cm)", 80, 0., 8.);

  hBeamVetoNumHitPlanes_ = H.DefineTH1D(hdir, "beamVetoNumHitPlanes", "beamVetoNumHitPlanes", 6, -0.5, 5.5);
  hHitPCsAterBeamVeto_ = H.DefineTH1D(hdir, "hitUpsteamPCsAfterBeamVeto", "hitUpsteamPCsAfterBeamVeto", 6, -0.5, 5.5);

  hWindowTimeBefore_ = H.DefineTH1D(hdir, "windowTimeBeforeCut", "Decay window start time, before time cut", 1000, 0., 10000.);
  hWindowTimeAfter_ = H.DefineTH1D(hdir, "windowTimeAfterCut", "Decay window start time, after time cut", 1000, 0., 10000.);

  hNumAfterTrigWindows_ = H.DefineTH1D(hdir, "numAfterTrigTimeWindows", "numAfterTrigTimeWindows", 10, -0.5, 9.5);
  hWindow2Time_ = H.DefineTH1D(hdir, "window2Time", "Second window start time", 1000, 0., 10000.);
  hWindow2dt_ = H.DefineTH1D(hdir, "window2dt", "Second window time - proton win time", 1000, 0., 10000.);

  //----------------------------------------------------------------
  hwidthPCall_.init(hdir+"/pcWidthAll", "pcwidth", 12, H, Conf);
  hwidthPCfiltered_.init(hdir+"/pcWidthFiltered", "pcwidth", 12, H, Conf);
  hwidthDCall_.init(hdir+"/dcWidthAll", "dcwidth", 44, H, Conf);
  hwidthDCfiltered_.init(hdir+"/dcWidthFiltered", "dcwidth", 44, H, Conf);

  if(fillXtalkPC_) {
    hAfterPulsingPCAll_.init(hdir+"/afterPulsingPCAll", 12, 160, H, Conf);
    hAfterPulsingPCFiltered_.init(hdir+"/afterPulsingPCFiltered", 12, 160, H, Conf);
    hXtalkSameWirePC_.init(hdir+"/xtalkSameWirePC", 12, 0, 0, H, Conf);
    hXtalk1PC_.init(hdir+"/xtalk1PC", 12, 1, 1, H, Conf);
    hXtalkPlanePC_.init(hdir+"/xtalkPlanePC", 12, 1, 999, H, Conf);
  }

  if(fillXtalkDC_) {
    hAfterPulsingDCAll_.init(hdir+"/afterPulsingDCAll", 44, 80, H, Conf);
    hAfterPulsingDCFiltered_.init(hdir+"/afterPulsingDCFiltered", 44, 80, H, Conf);
    hXtalkSameWireDC_.init(hdir+"/xtalkSameWireDC", 44, 0, 0, H, Conf);
    hXtalk1DC_.init(hdir+"/xtalk1DC", 44, 1, 1, H, Conf);
    hXtalkPlaneDC_.init(hdir+"/xtalkPlaneDC", 44, 1, 999, H, Conf);
  }

  //----------------------------------------------------------------
  hOccupancyPCAll_.init(hdir, "hitMapPCAll", 12, 160, H, Conf);
  hOccupancyDCAll_.init(hdir, "hitMapDCAll", 44, 80, H, Conf);
  haccidentalsTrig_.init(hdir+"/AccidentalsTrig", H, Conf);
  haccidentalsStop_.init(hdir+"/AccidentalsStop", H, Conf);

  winDCUnassignedAfterWindowing_.init(hdir+"/dcUnassignedWindowing", H, Conf);
  winDCUnassignedMuStop_.init(hdir+"/dcUnassignedMuStop", H, Conf);
  winDCUnassignedDnDecay_.init(hdir+"/dcUnassignedDnDecay", H, Conf);

  winTimeBeforeNoTrigWin_.init(hdir+"/winTime", "beforeNoTrigWin", H, Conf);
  winTimeMuStop_.init(hdir+"/winTime", "muStop", H, Conf);

  dioUp_.init(hdir+"/DIOUp", H, Conf,
              TimeWindow::UPSTREAM,
              Conf.read<double>(hdir+"/DIOUp/cutMinTime"));

  dioDn_.init(hdir+"/DIODn", H, Conf,
              TimeWindow::DOWNSTREAM,
              Conf.read<double>(hdir+"/DIODn/cutMinTime"));

  hdriftPCAll_.init(hdir+"/driftTimePCAll", H, 12, 1000./*ns*/,
                    Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtPC"),
                    Conf);
  hdriftPCFiltered_.init(hdir+"/driftTimePCFiltered", H, 12, 1000./*ns*/,
                         Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtPC"),
                         Conf);

  hdriftDCAll_.init(hdir+"/driftTimeDCAll", H, 44, 5000./*ns*/,
                    Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtDC"),
                    Conf);
  hdriftDCFiltered_.init(hdir+"/driftTimeDCFiltered", H, 44, 5000./*ns*/,
                         Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtDC"),
                         Conf);

  if(doMCTruth_) {
    hmuStopTruthAll_.init(H, hdir+"/MuStopTruthAll", *E.geo, Conf);
    hmuStopTruthAfterGaps_.init(H, hdir+"/MuStopTruthAfterGaps", *E.geo, Conf);
    hTruthAll_.init(H, hdir+"/MCTruthAll", Conf);
    hTruthMuStop_.init(H, hdir+"/MCTruthMuStop", Conf);
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
  if(anDnLateRes_.accepted) {
    hMeasuredMomentum_->Fill(anDnLateRes_.momentum);
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

  const TimeWindow& trigWin = wres.windows[wres.iTrigWin];

  //----------------
  // Process DC hits
  dcWindowing_.assignDCHits(filteredDCHits, &wres);
  winDCUnassignedAfterWindowing_.fill(wres);

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
  if(doMCTruth_) {
    hmuStopTruthAfterGaps_.fill(evt);
  }

  //----------------
  hMuonLastPlaneAfterGaps_->Fill(muonRange.max());
  if(muonRange.max() != 28) {
    return CUT_MUON_LAST_PLANE;
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
  winDCUnassignedMuStop_.fill(wres);
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

  //  ----------------
  // More stopped muon histos
  for(int i = 0; i + 1 < muonRange.segments().size(); ++i) {
    hStoppedMuonRangeGaps_->Fill(muonRange.segments()[i].max + 1, muonRange.segments()[i+1].min - 1);
  }
  for(int iplane = 1; iplane <= muonRange.max(); ++iplane) {
    if(muonGlobalClusters[iplane].empty()) {
      hStoppedMuonMissingPlanes_->Fill(iplane);
    }
  }

  // return CUTS_MUSTOP_ACCEPTED;

  //----------------------------------------------------------------
  bool have_downstream_pchits = false;

  if(1 + wres.iTrigWin < wres.windows.size()) {
    // Look for PC hits in the first after-trig time window.
    const TimeWindow& win = wres.windows[1+wres.iTrigWin];
    for(unsigned i=0; i < win.pcHits.size(); ++i) {
      if(win.pcHits[i]->plane() > 6) {
        have_downstream_pchits = true;
        break;
      }
    }
  }
  // FIXME: histo of #dn PC, all PC, all DC hits in the window
  // For now after-trig window: num unassigned after-trig DC hits.


  if(!have_downstream_pchits) {
    return CUT_DOWNSTREAM_PCHITS;
  }

  const unsigned iDecayWin = 1 + wres.iTrigWin;
  const TimeWindow& decayWindow = wres.windows[iDecayWin];

  //----------------------------------------------------------------
  // Veto accidental beam particles (also upstream DIOs and some protons)

  const ClustersByPlane& protonGlobalClusters = afterTrigClusters[0];

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

  hWindowTimeBefore_->Fill(decayWindow.tstart);
  if( (decayWindow.tstart < cutWinTimeMin_) || (cutWinTimeMax_ < decayWindow.tstart)) {
    return CUT_WIN_TIME;
  }
  hWindowTimeAfter_->Fill(decayWindow.tstart);

  //----------------------------------------------------------------
  // Deal with multiple after-trigger time windows
  hNumAfterTrigWindows_->Fill(afterTrigClusters.size());
  if(afterTrigClusters.size() > 1) {
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
  winDCUnassignedDnDecay_.fill(wres);
  // What do events with many unassigned hits look like?
  if(false && (wres.unassignedDCHits.size()>10) ) {
    std::cout<<"Unassigned DC hits debug: run "<<evt.nrun<<" event "<<evt.nevt<<std::endl;
    std::cout<<"unassigned  ("<<wres.unassignedDCHits.size()<<"): ";
    for(unsigned i=0; i<wres.unassignedDCHits.size(); ++i) {
      std::cout<<*wres.unassignedDCHits[i]<<" ";
    }
    std::cout<<"\nwindows:\n";
    for(unsigned iwin=0; iwin<wres.windows.size(); ++iwin) {
      std::cout<<"iwin="<<iwin<<" tstart "<<wres.windows[iwin].tstart<<": PC = "<<wres.windows[iwin].pcHits<<std::endl;
      std::cout<<"iwin="<<iwin<<" tstartDC "<<wres.windows[iwin].tstartDC<<": DC = "<<wres.windows[iwin].dcHits<<std::endl;
    }
  }

  // Other analyses....


  return CUTS_DOWNSTREAM_ACCEPTED;
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
