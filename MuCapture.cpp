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
#include "MuCapUtilities.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//  Backdoor for the hit-based hotspot debugging
TimeWindowingResults *gwres;

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
    std::cout<<"  STARTING !!!!!!!!!!!!!!!!!"<<std::endl;

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

  //       --------- Histograms initialization ---------          //
  const std::string chanstr = Conf.read<std::string>("MuCapture/channels/versions");
  const std::vector<std::string> channelNames = StrToStrVect(chanstr);
  for(unsigned i=0; i<channelNames.size(); ++i) {
    std::cout<<"MuCapture: adding channel "<<channelNames[i]<<std::endl;
    channels_.push_back(HistMuCapAnalysisChannels());
    channels_.back().init(H, "MuCapture/channels", channelNames[i], *E.geo, Conf);
  }

//minhist:  rooUnfHits_.init(H, "MuCapture", Conf);

  //----------------------------------------------------------------
  dnPosTracks_.init(hdir+"/dnPosTracks", H, *E.geo, Conf, "pos", TimeWindow::DOWNSTREAM, &anDnLateRes_);
  //dnDIOVetoTracks_.init(hdir+"/dnDIOVetoTracks", H, *E.geo, Conf, "dioVeto", TimeWindow::DOWNSTREAM);
  dnDIONormTracks_.init(hdir+"/dnDIONormTracks", H, *E.geo, Conf, "dioNorm", TimeWindow::DOWNSTREAM);

//minhist:  dnPosTrkContainment_.init(hdir+"/dnPosTrkContainment", H, *E.geo, Conf);
//minhist:  dnPosTrkRangeStudies_.init(H, hdir+"/dnPosTrkRangeStudies", *E.geo, Conf);
//minhist:  dnPosTrkClassification_.init(H, hdir+"/dnPosTrkClassification", *E.geo, Conf);
//minhist:
//minhist:  hProtonPID_.init(hdir+"/posTrackPID", H, Conf);
//minhist:  hContainedProtonPID_.init(hdir+"/posContainedTrackPID", H, Conf);

  //----------------------------------------------------------------
  pcXTalkProcessor_ = makeXTalkPreprocessor(hdir+"/hitLevel/XTalk", WirePlane::PC, H, *E.geo, Conf);
  dcXTalkProcessor_ = makeXTalkPreprocessor(hdir+"/hitLevel/XTalk", WirePlane::DC, H, *E.geo, Conf);

  pcHitProcessor_ = makeTDCHitPreprocessor(hdir+"/hitLevel/HitPreProc", WirePlane::PC, H, *E.geo, Conf);
  dcHitProcessor_ = makeTDCHitPreprocessor(hdir+"/hitLevel/HitPreProc", WirePlane::DC, H, *E.geo, Conf);

//minhist:  hwidthPCall_.init(hdir+"/hitLevel/pcWidthAll", "pcwidth", 12, H, Conf);
//minhist:  hwidthPCfiltered_.init(hdir+"/hitLevel/pcWidthFiltered", "pcwidth", 12, H, Conf);
//minhist:  hwidthDCall_.init(hdir+"/hitLevel/dcWidthAll", "dcwidth", 44, H, Conf);
//minhist:  hwidthDCfiltered_.init(hdir+"/hitLevel/dcWidthFiltered", "dcwidth", 44, H, Conf);
//minhist:
//minhist:  hwidthMuHits_.init(H, hdir+"/hitLevel/muonHitWidth", *E.geo, Conf);
//minhist:  hwidthDIOHits_.init(H, hdir+"/hitLevel/dioHitWidth", *E.geo, Conf);

  if(fillXtalkPC_) {
    hAfterPulsingPCAll_.init(hdir+"/hitLevel/afterPulsingPCAll", 12, 160, H, Conf);
    hAfterPulsingPCFiltered_.init(hdir+"/hitLevel/afterPulsingPCFiltered", 12, 160, H, Conf);
    hXtalkSameWirePC_.init(hdir+"/hitLevel/xtalkSameWirePC", 12, 0, 0, H, Conf);
    hXtalk1PC_.init(hdir+"/hitLevel/xtalk1PC", 12, 1, 1, H, Conf);
    hXtalkPlanePC_.init(hdir+"/hitLevel/xtalkPlanePC", 12, 1, 999, H, Conf);
    hXT2PlanePC_.init(hdir+"/hitLevel/xt2PlanePC", 12, 1, 999, 40., H, Conf);
  }

  if(fillXtalkDC_) {
    hAfterPulsingDCAll_.init(hdir+"/hitLevel/afterPulsingDCAll", 44, 80, H, Conf);
    hAfterPulsingDCFiltered_.init(hdir+"/hitLevel/afterPulsingDCFiltered", 44, 80, H, Conf);
    hXtalkSameWireDC_.init(hdir+"/hitLevel/xtalkSameWireDC", 44, 0, 0, H, Conf);
    hXtalk1DC_.init(hdir+"/hitLevel/xtalk1DC", 44, 1, 1, H, Conf);
    hXtalkPlaneDC_.init(hdir+"/hitLevel/xtalkPlaneDC", 44, 1, 999, H, Conf);
    hXT2PlaneDC_.init(hdir+"/hitLevel/xt2PlaneDC", 44, 1, 999, 50., H, Conf);
  }

  hXT3_.init(hdir+"/hitLevel/xt3", H, Conf);
  hXT4_.init(hdir+"/hitLevel/xt4", H, Conf);
  hXT5_.init(hdir+"/hitLevel/xt5", H, *E.geo, Conf);
//minhist:  h200ns_.init(hdir+"/hitLevel/peak200ns", H, *E.geo, Conf);
//minhist:  h200nsDIO_.init(hdir+"/hitLevel/peak200nsDIO", H, *E.geo, Conf);

  hOccupancyPCAll_.init(hdir+"/hitLevel", "hitMapPCAll", 12, 160, H, Conf);
  hOccupancyDCAll_.init(hdir+"/hitLevel", "hitMapDCAll", 44, 80, H, Conf);

  //----------------------------------------------------------------
  uvOutFileName_ = Conf.read<std::string>(hdir+"/uvOutFileName", "");
  if(!uvOutFileName_.empty()) {
    uvOutFile_.open(uvOutFileName_.c_str());
    if(!uvOutFile_) {
      throw std::runtime_error("Error opening output file "+uvOutFileName_);
    }
  }

  //----------------------------------------------------------------
  commonSkimOutFileName_ = Conf.read<std::string>("MuCapture/commonSkimOutFileName", "");
  if(!commonSkimOutFileName_.empty()) {
    commonSkimOutFile_.open(commonSkimOutFileName_.c_str());
    if(!commonSkimOutFile_) {
      throw std::runtime_error("Error opening output file "+commonSkimOutFileName_);
    }
  }

  //----------------------------------------------------------------
  pcWindowing_.init(H, hdir+"/WindowingPC", *E.geo, Conf);
  dcWindowing_.init(H, hdir+"/WindowingDC", *E.geo, Conf);

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

  hMuonRMaxDCvsPC_ = H.DefineTH2D(hdir, "muonRMaxDCvsPC", "Muon Rmax DC vs PC", 45, 0., 18., 45, 0., 18.);
  hMuonRMaxDCvsPC_->SetOption("colz");

  hMuonRMaxVsRStop_ = H.DefineTH2D(hdir, "muonRMaxVsRStop", "Muon Rmax vs Rstop", 45, 0., 18., 45, 0., 18.);
  hMuonRMaxVsRStop_->SetOption("colz");

  hMuonRangeUnrestricted_.init(H, hdir+"/muonRangeUR", *E.geo, Conf);
  hMuonRangeRestricted_.init(H, hdir+"/muonRangeRestricted", *E.geo, Conf);

  // Make the bin size half a cell
  hMuStopUVCell_ = H.DefineTH2D(hdir, "MuStopUVCell", "Muon stop V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);
  hMuStopUVCell_->SetOption("colz");
  hMuStopUVPos_ = H.DefineTH2D(hdir, "MuStopUVPos", "Muon stop V vs U position (cm)", 201, -10.05, +10.05, 201, -10.05, +10.05 );
  hMuStopUVPos_->SetOption("colz");
  hMuStopRadius_ = H.DefineTH1D(hdir, "MuStopRadius", "Muon stop R (cm)", 80, 0., 8.);

  pactCut_.init(H, *E.geo, Conf);

  hStoppedMuonRangeGaps_ = H.DefineTH2D(hdir+"/AcceptedMuStop", "muStopRangeGaps", "Stopped muon range gap end vs start", 57, -0.5, 56.5, 57, -0.5, 56.5);
  hStoppedMuonRangeGaps_->SetOption("colz");

  hStoppedMuonMissingPlanes_ = H.DefineTH1D(hdir+"/AcceptedMuStop", "muStopMissingPlanes", "Stopped muon missing planes", 56, 0.5, 56.5);

//mincuts:  hBeamVetoNumHitPlanes_ = H.DefineTH1D(hdir, "beamVetoNumHitPlanes", "beamVetoNumHitPlanes", 6, -0.5, 5.5);
//mincuts:  hHitPCsAterBeamVeto_ = H.DefineTH1D(hdir, "hitUpsteamPCsAfterBeamVeto", "hitUpsteamPCsAfterBeamVeto", 6, -0.5, 5.5);

  hWindowTimeFirst_ = H.DefineTH1D(hdir, "windowTimeFirst", "Window start time, first after-trig window", 1000, 0., 10000.);
  hWindowTimeAll_ = H.DefineTH1D(hdir, "windowTimeAll", "Window start time, all non-trig windows", 16000, -6000, 10000.);

  hBestTrackWindDt_ = H.DefineTH1D(hdir, "bestTrackWinDt", "Best track-win time", 501, -250.5, 250.5);

  hNumAfterTrigWindows_ = H.DefineTH1D(hdir, "numAfterTrigTimeWindows", "numAfterTrigTimeWindows", 10, -0.5, 9.5);
//mincuts:  hWindow2Time_ = H.DefineTH1D(hdir, "window2Time", "Second window start time", 1000, 0., 10000.);
//mincuts:  hWindow2dt_ = H.DefineTH1D(hdir, "window2dt", "Second window time - proton win time", 1000, 0., 10000.);


//mincuts:  hPosNegMom_ = H.DefineTH2D(hdir, "posnegmom", "p(-) vs p(+)", 300, 0., 300., 300, 0., 300.);
//mincuts:  hPosNegMom_->SetOption("colz");
//mincuts:
//mincuts:  hPosNegCosth_ = H.DefineTH2D(hdir, "posnegcosth", "cos(-) vs cos(+)", 100, -1., 1., 100, -1., 1.);
//mincuts:  hPosNegCosth_->SetOption("colz");
//mincuts:
//mincuts:  hVetoedPCosth_ = H.DefineTH2D(hdir, "vetoedpcosth", "cos(theta) vs p of vetoed tracks", 300, 0., 300., 100, -1., 1.);
//mincuts:  hVetoedPCosth_->SetOption("colz");
//mincuts:
//mincuts:  hVetoingPCosth_ = H.DefineTH2D(hdir, "vetoingpcosth", "cos(theta) vs p of veto tracks", 300, 0., 300., 100, -1., 1.);
//mincuts:  hVetoingPCosth_->SetOption("colz");

  //----------------------------------------------------------------
//minhist:  haccidentalsTrig_.init(hdir+"/accidentals/trig", H, Conf);
//minhist:  haccidentalsStop_.init(hdir+"/accidentals/stop", H, Conf);
//minhist:
//minhist:  winDCUnassignedAfterWindowing_.init(hdir+"/unassignedDC/windowing", H, Conf);
//minhist:  winDCUnassignedMuStop_.init(hdir+"/unassignedDC/mustop", H, Conf);
//minhist:  winDCUnassignedDnDecay_.init(hdir+"/unassignedDC/dndecay", H, Conf);
//minhist:
//minhist:  dioUp_.init(hdir+"/chamberEff/DIOUp", H, Conf,
//minhist:              TimeWindow::UPSTREAM,
//minhist:              Conf.read<double>(hdir+"/DIOUp/cutMinTime"));
//minhist:
//minhist:  dioDn_.init(hdir+"/chamberEff/DIODn", H, Conf,
//minhist:              TimeWindow::DOWNSTREAM,
//minhist:              Conf.read<double>(hdir+"/DIODn/cutMinTime"));

  hdriftPCAll_.init(hdir+"/chamberEff/driftTimePCAll", H, 12, 200./*ns*/,
                    Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtPC"),
                    Conf);
  hdriftPCFiltered_.init(hdir+"/chamberEff/driftTimePCFiltered", H, 12, 200./*ns*/,
                         Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtPC"),
                         Conf);

  hdriftDCAll_.init(hdir+"/chamberEff/driftTimeDCAll", H, 44, 4000./*ns*/,
                    Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtDC"),
                    Conf);
  hdriftDCFiltered_.init(hdir+"/chamberEff/driftTimeDCFiltered", H, 44, 4000./*ns*/,
                         Conf.read<double>(hdir+"/HistDriftTime/cutEffTrackHitDtDC"),
                         Conf);

  if(doMCTruth_) {
    hTruthAll_.init(H, hdir+"/MCTruth/all", *E.geo, Conf);
    hTruthLastPlane_.init(H, hdir+"/MCTruth/lastPlane", *E.geo, Conf);
    hTruthPACT_.init(H, hdir+"/MCTruth/pact", *E.geo, Conf);
//mincuts:    hTruthDnPC_.init(H, hdir+"/MCTruth/dnPC", *E.geo, Conf);
//minhist:    hTruthDnCandidate_.init(H, hdir+"/MCTruth/dnCandidate", *E.geo, Conf);
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
//minhist:  if(doMCTruth_ && c > CUT_EVENT_NUMBER) {
//minhist:    rooUnfHits_.FillAndMiss();
//minhist:  }

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

  for(unsigned i=0; i<channels_.size(); i++) {
    channels_[i].fillReferenceSample(evt);
  }

//minhist:  rooUnfHits_.SaveEventVariables(evt);

  TDCHitPreprocessing::Hits allPCHitsBuf(evt.pc_hits());
  const TDCHitWPPtrCollection& allPCHits = allPCHitsBuf.get();

//debug:  if(true) {
//debug:    TDCHitWPPtrCollection localhits = allPCHits;
//debug:    std::sort(localhits.begin(), localhits.end(), TDCHitWPCmpTime());
//debug:    std::stable_sort(localhits.begin(), localhits.end(), TDCHitWPCmpGeom());
//debug:    // here we have localhits sorted by cell then time
//debug:      std::cout<<"#================================================================\n";
//debug:      std::cout<<"# event "<<evt.nrun<<" "<<evt.nevt<<std::endl;
//debug:      std::cout<<"all PC hits = "<<localhits<<std::endl;
//debug:  }

  TDCHitWPPtrCollection nxtPCHits;
  pcXTalkProcessor_->process(&nxtPCHits, allPCHits);

  TDCHitWPPtrCollection filteredPCHits;
  pcHitProcessor_->process(&filteredPCHits, nxtPCHits);

//minhist:  hwidthPCall_.fill(allPCHits);
//minhist:  hwidthPCfiltered_.fill(filteredPCHits);

  if(fillXtalkPC_) {
    hAfterPulsingPCAll_.fill(allPCHits);
    hAfterPulsingPCFiltered_.fill(filteredPCHits);
    hXtalkSameWirePC_.fill(allPCHits);
    hXtalk1PC_.fill(allPCHits);
    hXtalkPlanePC_.fill(allPCHits);
    hXT2PlanePC_.fill(allPCHits);
  }

  TDCHitPreprocessing::Hits allDCHitsBuf(evt.dc_hits());
  const TDCHitWPPtrCollection& allDCHits = allDCHitsBuf.get();

//debug:  if(true) {
//debug:    TDCHitWPPtrCollection localhits = allDCHits;
//debug:    std::sort(localhits.begin(), localhits.end(), TDCHitWPCmpTime());
//debug:    std::stable_sort(localhits.begin(), localhits.end(), TDCHitWPCmpGeom());
//debug:    // here we have localhits sorted by cell then time
//debug:      std::cout<<"#================================================================\n";
//debug:      std::cout<<"# event "<<evt.nrun<<" "<<evt.nevt<<std::endl;
//debug:      std::cout<<"all DC hits = "<<localhits<<std::endl;
//debug:  }

  TDCHitWPPtrCollection nxtDCHits;
  dcXTalkProcessor_->process(&nxtDCHits, allDCHits);

  TDCHitWPPtrCollection filteredDCHits;
  dcHitProcessor_->process(&filteredDCHits, nxtDCHits);

//minhist:  hwidthDCall_.fill(allDCHits);
//minhist:  hwidthDCfiltered_.fill(filteredDCHits);

  if(fillXtalkDC_) {
    hAfterPulsingDCAll_.fill(allDCHits);
    hAfterPulsingDCFiltered_.fill(filteredDCHits);
    hXtalkSameWireDC_.fill(allDCHits);
    hXtalk1DC_.fill(allDCHits);
    hXtalkPlaneDC_.fill(allDCHits);
    hXT2PlaneDC_.fill(allDCHits);
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
  gwres = &wres;
  pcWindowing_.assignPCHits(filteredPCHits, &wres);
  assert(filteredPCHits.size() ==
         std::accumulate(wres.windows.begin(), wres.windows.end(), WinHitCounter()).numPCHits);

  if(wres.iTrigWin == -1u) {
    return CUT_NOTRIGWIN;
  }

//minhist:  haccidentalsTrig_.fill(wres);

  if(wres.iTrigWin > 0) {
    // Trig time is 0, dt from that rather than from less precise trigWin time
    const double dt = wres.windows[wres.iTrigWin - 1].tstart;
    hPCPreTrigSeparation_->Fill(dt);
    if(std::abs(dt) < winPCPreTrigSeparation_) {
      // Without this cut a pre-trigger particle can create a time
      // window that would "eat" hits belonging to a capture proton in
      // our windowing procedure.  Then the proton containment cut
      // would produce a wrong result.
      return CUT_PCWIN_TRIGSEPPAST;
    }
  }

  const TimeWindow& trigWin = wres.windows[wres.iTrigWin];

  //----------------
  // Process DC hits
  dcWindowing_.assignDCHits(filteredDCHits, &wres);
//minhist:  winDCUnassignedAfterWindowing_.fill(wres);

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
  const double muonRMaxPC = MuCapUtilities::computeHitRMax(trigWin.pcHits, WirePlane::PC, *evt.geo);
  const double muonRMaxDC = MuCapUtilities::computeHitRMax(trigWin.dcHits, WirePlane::DC, *evt.geo);
  const double muonRMax = std::max(muonRMaxPC, muonRMaxDC);
  hMuonRMaxDCvsPC_->Fill(muonRMaxPC, muonRMaxDC);
  hMuonRangeUnrestricted_.fill(muonRange, muonGlobalClusters, evt);
  if(muonRMax < 4.0 /*cm, hardcoded same value as for MOFIA's H34159 in dplot.f90  */) {
    hMuonRangeRestricted_.fill(muonRange, muonGlobalClusters, evt);
  }

  if(muonRange.max() != 28) {
    return CUT_MUON_LAST_PLANE;
  }

  if(doMCTruth_) {
    hTruthLastPlane_.fill(evt);
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
  hMuonRMaxVsRStop_->Fill(muStopRadius, muonRMax);
  if(muStopRadius > muStopRMax_) {
    return CUT_MUSTOP_UV;
  }

  //----------------------------------------------------------------
  if(1 != pactCut_.quadrant(muonGlobalClusters[27].front(), muonGlobalClusters[28].front(), evt)) {
    return CUT_MUSTOP_PACT;
  }

  if(doMCTruth_) {
    hTruthPACT_.fill(evt);
  }

  //----------------------------------------------------------------
  // Here we have an accepted muon stop

  // Call the subanalyses

//minhist:  winDCUnassignedMuStop_.fill(wres);
//minhist:  haccidentalsStop_.fill(wres);

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

//mincuts:  //----------------------------------------------------------------
//mincuts:  bool have_downstream_pchits = false;
//mincuts:
//mincuts:  if(1 + wres.iTrigWin < wres.windows.size()) {
//mincuts:    // Look for PC hits in the first after-trig time window.
//mincuts:    const TimeWindow& win = wres.windows[1+wres.iTrigWin];
//mincuts:    for(unsigned i=0; i < win.pcHits.size(); ++i) {
//mincuts:      if(win.pcHits[i]->plane() > 6) {
//mincuts:        have_downstream_pchits = true;
//mincuts:        break;
//mincuts:      }
//mincuts:    }
//mincuts:  }
//mincuts:  // FIXME: histo of #dn PC, all PC, all DC hits in the window
//mincuts:  // For now after-trig window: num unassigned after-trig DC hits.
//mincuts:
//mincuts:
//mincuts:  if(!have_downstream_pchits) {
//mincuts:    return CUT_DOWNSTREAM_PCHITS;
//mincuts:  }
//mincuts:  if(doMCTruth_) {
//mincuts:    hTruthDnPC_.fill(evt);
//mincuts:  }
//mincuts:
//mincuts:  const unsigned iDecayWin = 1 + wres.iTrigWin;
//mincuts:  const TimeWindow& decayWindow = wres.windows[iDecayWin];
//mincuts:
//mincuts:  //----------------------------------------------------------------
//mincuts:  // Veto accidental beam particles (also upstream DIOs and some protons)
//mincuts:
//mincuts:  const ClustersByPlane& afterTrig0GlobalCl = afterTrigClusters[0];
//mincuts:
//mincuts:  int numBeamVetoHitPlanes(0);
//mincuts:  for(int i=1; i<=4; ++i) {
//mincuts:    if(!afterTrig0GlobalCl[i].empty()) {
//mincuts:      ++numBeamVetoHitPlanes;
//mincuts:    }
//mincuts:  }
//mincuts:
//mincuts:  hBeamVetoNumHitPlanes_->Fill(numBeamVetoHitPlanes);
//mincuts:  if(numBeamVetoHitPlanes > cutBeamVetoMaxPCplanes_) {
//mincuts:    return CUT_BEAM_VETO;
//mincuts:  }
//mincuts:
//mincuts:  for(int i=1; i<=4; ++i) {
//mincuts:    if(!afterTrig0GlobalCl[i].empty()) {
//mincuts:      hHitPCsAterBeamVeto_->Fill(i);
//mincuts:    }
//mincuts:  }
//mincuts:

  //----------------------------------------------------------------
  // Trig time is 0, dt from that rather than from less precise trigWin time

  for(int i=0; i < wres.windows.size(); ++i) {
    if(i != wres.iTrigWin) {
      hWindowTimeAll_->Fill(wres.windows[i].tstart);
    }
  }
  if(wres.iTrigWin + 1 < wres.windows.size()) {
    hWindowTimeFirst_->Fill(wres.windows[1 + wres.iTrigWin].tstart);
  }

  //----------------------------------------------------------------
  // Deal with multiple after-trigger time windows
  hNumAfterTrigWindows_->Fill(wres.windows.size() - wres.iTrigWin - 1);
//mincuts:  if(afterTrigClusters.size() > 1) {
//mincuts:    const TimeWindow& win1 = wres.windows[wres.iTrigWin + 1];
//mincuts:    const TimeWindow& win2 = wres.windows[wres.iTrigWin + 2];
//mincuts:    const double dt2 = win2.tstart - win1.tstart;
//mincuts:    hWindow2Time_->Fill(win2.tstart);
//mincuts:    hWindow2dt_->Fill(dt2);
//mincuts:    if(dt2 < cutMultiwinNextdt_) {
//mincuts:      return CUT_MULTIWIN_NEXTDT;
//mincuts:    }
//mincuts:  }

  //----------------------------------------------------------------
  //  The pre-selection for downstream decays/captures is passed here
  // select the normalization sample
  const int iDIONorm = dnDIONormTracks_.process(evt, muStop);

  //const int iNegTrack = dnDIOVetoTracks_.process(evt, muStop);
  const int iNegTrack = -1; // disable the veto

  const int iPosTrack = dnPosTracks_.process(evt, muStop);

  const int iWinPosTrack = findTrackWindow(evt, iPosTrack, wres.windows);
  if((iPosTrack > -1 ) && (iWinPosTrack - wres.iTrigWin < 1)) {
    //throw std::runtime_error("MuCapture: iWinPosTrack - wres.iTrigWin < 1");
    std::cerr<<"MuCapture: iWinPosTrack - wres.iTrigWin < 1: iWinPosTrack="<<iWinPosTrack
             <<", iTrigWin = "<<wres.iTrigWin
             <<", numWindows = "<<wres.windows.size()
             <<std::endl;

    std::cerr<<"track time = "<<evt.hefit_time[iPosTrack] <<std::endl;
    std::cerr<<"wres = "<<wres<<std::endl;

  }

  ClustersByPlane protonGlobalClusters;
  if(iWinPosTrack > -1) {
    const TimeWindow& win = wres.windows[iWinPosTrack];
    const ClustersByPlane pcClusters = constructPlaneClusters(12, win.pcHits);
    const ClustersByPlane dcClusters = constructPlaneClusters(44, win.dcHits);
    protonGlobalClusters = globalPlaneClusters(pcClusters, dcClusters);
  }

//minhist:  const bool isPosTrackContained = dnPosTrkContainment_.contained(evt, iPosTrack, protonGlobalClusters);
//minhist:  const double rangePIDVar = ((iPosTrack != -1)&& isPosTrackContained) ? hContainedProtonPID_.fill(evt, iPosTrack, protonGlobalClusters) : 0.;

  for(unsigned i=0; i<channels_.size(); i++) {
    channels_[i].fill(evt, iPosTrack, iNegTrack, protonGlobalClusters);
  }

//minhist:  rooUnfHits_.FillMeasured(evt, iPosTrack, iNegTrack, isPosTrackContained, rangePIDVar);
//minhist:  if(channels_[0].referenceSampleAccepted()) {
//minhist:    rooUnfHits_.Fill(evt, iPosTrack, iNegTrack, isPosTrackContained, rangePIDVar);
//minhist:  }

  //----------------------------------------------------------------
  // Fill extra distributions

//minhist:  if(iNegTrack == -1) { // Veto DIO events
//minhist:    if(iPosTrack != -1) { // Got a reconstructed capture track
//minhist:      hProtonPID_.fill(evt, iPosTrack, protonGlobalClusters);
//minhist:    }
//minhist:  }
//minhist:
//minhist:  if(doMCTruth_) {
//minhist:    hTruthDnCandidate_.fill(evt);
//minhist:  }
//minhist:
//minhist:  hwidthMuHits_.fill(evt, muonGlobalClusters);
//minhist:  if(iDIONorm != -1) {
//minhist:    const int iWinDIONorm = findTrackWindow(evt, iDIONorm, wres.windows);
//minhist:    if(iWinDIONorm - wres.iTrigWin < 1) {
//minhist:      throw std::runtime_error("MuCapture: iWinDIONorm - wres.iTrigWin < 1");
//minhist:    }
//minhist:
//minhist:    ClustersByPlane dioGlobalClusters;
//minhist:    const TimeWindow& win = wres.windows[iWinDIONorm];
//minhist:    const ClustersByPlane pcClusters = constructPlaneClusters(12, win.pcHits);
//minhist:    const ClustersByPlane dcClusters = constructPlaneClusters(44, win.dcHits);
//minhist:    dioGlobalClusters = globalPlaneClusters(pcClusters, dcClusters);
//minhist:
//minhist:    hwidthDIOHits_.fill(evt, dioGlobalClusters);
//minhist:    h200nsDIO_.fill(evt, iDIONorm, dioGlobalClusters);
//minhist:  }
//minhist:
//minhist:  // Do we have any ambiguous events (both "DIO" and "proton")?
//minhist:  if((iNegTrack != -1)&&(iPosTrack != -1)) {
//minhist:    hPosNegMom_->Fill(evt.ptot[iPosTrack], evt.ptot[iNegTrack]);
//minhist:    hPosNegCosth_->Fill(evt.costh[iPosTrack], evt.costh[iNegTrack]);
//minhist:    hVetoedPCosth_->Fill(evt.ptot[iPosTrack], evt.costh[iPosTrack]);
//minhist:    hVetoingPCosth_->Fill(evt.ptot[iNegTrack], evt.costh[iNegTrack]);
//minhist:  }
//minhist:
//minhist:  winDCUnassignedDnDecay_.fill(wres);
//minhist:  // What do events with many unassigned hits look like?
//minhist:  if(false && (wres.unassignedDCHits.size()>10) ) {
//minhist:    std::cout<<"Unassigned DC hits debug: run "<<evt.nrun<<" event "<<evt.nevt<<std::endl;
//minhist:    std::cout<<"unassigned  ("<<wres.unassignedDCHits.size()<<"): ";
//minhist:    for(unsigned i=0; i<wres.unassignedDCHits.size(); ++i) {
//minhist:      std::cout<<*wres.unassignedDCHits[i]<<" ";
//minhist:    }
//minhist:    std::cout<<"\nwindows:\n";
//minhist:    for(unsigned iwin=0; iwin<wres.windows.size(); ++iwin) {
//minhist:      std::cout<<"iwin="<<iwin<<" tstart "<<wres.windows[iwin].tstart<<": PC = "<<wres.windows[iwin].pcHits<<std::endl;
//minhist:      std::cout<<"iwin="<<iwin<<" tstartDC "<<wres.windows[iwin].tstartDC<<": DC = "<<wres.windows[iwin].dcHits<<std::endl;
//minhist:    }
//minhist:  }
//minhist:
//minhist:  if(isPosTrackContained) {
//minhist:    dnPosTrkRangeStudies_.fill(evt, iPosTrack, protonGlobalClusters);
//minhist:  }
//minhist:
//minhist:  dnPosTrkClassification_.fill(evt, iPosTrack, protonGlobalClusters);

//minhist:
//minhist:  int idio = dioUp_.process(evt, muStop); // NOTE: upstream fit disabled in  MOFIA
//minhist:  if(idio >= 0) {
//minhist:    hdriftPCAll_.fill(evt, idio, allPCHits);
//minhist:    hdriftPCFiltered_.fill(evt, idio, filteredPCHits);
//minhist:    hdriftDCAll_.fill(evt, idio, allDCHits);
//minhist:    hdriftDCFiltered_.fill(evt, idio, filteredDCHits);
//minhist:    if(uvOutFile_) {
//minhist:      uvOutFile_<<std::fixed<<std::showpos<<evt.hefit_u0[idio]<<"\t"<<evt.hefit_v0[idio]<<std::endl;
//minhist:    }
//minhist:  }
//minhist: idio = dioDn_.process(evt, muStop);

  int idio = iDIONorm;
  if(idio >= 0) {
    hdriftPCAll_.fill(evt, idio, allPCHits);
    hdriftPCFiltered_.fill(evt, idio, filteredPCHits);
    hdriftDCAll_.fill(evt, idio, allDCHits);
    hdriftDCFiltered_.fill(evt, idio, filteredDCHits);
    if(uvOutFile_) {
      uvOutFile_<<std::fixed<<std::showpos<<evt.hefit_u0[idio]<<"\t"<<evt.hefit_v0[idio]<<std::endl;
    }
  }

  if(iPosTrack != -1) {
    hXT3_.fill(evt, iPosTrack, protonGlobalClusters);
    hXT4_.fill(evt, iPosTrack);
    hXT5_.fill(evt, iPosTrack);
//minhist:    h200ns_.fill(evt, iPosTrack, protonGlobalClusters);
  }

  //----------------
  if(commonSkimOutFile_) {
    commonSkimOutFile_<<evt.nrun<<" "<<evt.nevt<<std::endl;
  }

  return CUTS_MUSTOP_ACCEPTED; // CUTS_DOWNSTREAM_ACCEPTED
}

//================================================================
int MuCapture::findTrackWindow(const EventClass& evt,
                               int itrack,
                               const TimeWindowCollection& windows)
{
  int iwin = -1;

  if(itrack > -1) {
    double best_dt = std::numeric_limits<double>::max();
    for(int i=0; i<windows.size(); ++i) {
      double dt = evt.hefit_time[itrack] -  windows[i].tstart;
      if(std::abs(dt) < best_dt) {
        best_dt = dt;
        iwin = i;
      }
    }
    hBestTrackWindDt_->Fill(best_dt);
  }

  return iwin;
}

//================================================================
TDCHitPreprocessing::IProcessor *
MuCapture::makeXTalkPreprocessor(const std::string& topdir,
                                 WirePlane::DetType d,
                                 HistogramFactory& hf,
                                 const DetectorGeo& geom,
                                 ConfigFile& conf)
{
  const bool doXTalk = conf.read<bool>("MuCapture/HitPreproc/"+WirePlane::detName(d)+"/applyXTalk");
  if(doXTalk) {
    return new TDCHitPreprocessing::
      MOFIA_XTalkDiscarder(topdir, d, hf, geom, conf);
  }
  else {
    return new TDCHitPreprocessing::PassThrough();
  }
}

//================================================================
TDCHitPreprocessing::IProcessor *
MuCapture::makeTDCHitPreprocessor(const std::string& topdir,
                                  WirePlane::DetType d,
                                  HistogramFactory& hf,
                                  const DetectorGeo& geom,
                                  ConfigFile& conf)
{
  const std::string proc = conf.read<std::string>("MuCapture/HitPreproc/"+WirePlane::detName(d)+"/processor");
  if(proc == "NarrowHitDiscarder") {
    return new TDCHitPreprocessing::NarrowHitDiscarder(topdir, d, hf, geom, conf);
  }
  if(proc == "MOFIA_XTalkDiscarder") {
    return new TDCHitPreprocessing::MOFIA_XTalkDiscarder(topdir, d, hf, geom, conf);
  }
  if(proc == "SameWireHitDiscarder") {
    return new TDCHitPreprocessing::SameWireHitDiscarder(topdir, d, hf, geom, conf);
  }
  if(proc == "PassThrough") {
    return new TDCHitPreprocessing::PassThrough();
  }
  throw std::runtime_error("MuCapture::makeTDCHitPreprocessor(): unknown processor name \""+proc+"\"");
}

//================================================================
