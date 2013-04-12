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

#include "PlaneRange.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//================================================================
namespace { // local helpers
  //----------------------------------------------------------------
  std::pair<int,int> uvminmax(const TDCHitWPPtrCollection& hits, const std::set<int>& onlyPlanes = std::set<int>()) {
    int cmin=999, cmax=-1;
    for(unsigned i=0; i<hits.size(); ++i) {
      if(onlyPlanes.empty() || (onlyPlanes.find(hits[i]->plane) != onlyPlanes.end())) {
        if(hits[i]->cell < cmin) {
          cmin = hits[i]->cell;
        }
        if(cmax < hits[i]->cell) {
          cmax = hits[i]->cell;
        }
      }
    }
    return std::make_pair(cmin,cmax);
  }

  //----------------------------------------------------------------
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
  cutMinTDCWidthPC_ = Conf.read<double>("MuCapture/cutMinTDCWidthPC");
  cutMinTDCWidthDC_ = Conf.read<double>("MuCapture/cutMinTDCWidthDC");
  winPCLength_ = Conf.read<double>("MuCapture/winPCLength");
  winPCSeparation_ = Conf.read<double>("MuCapture/winPCSeparation");
  winDCLength_ = Conf.read<double>("MuCapture/winDCLength");
  winDCEarlyMargin_ = Conf.read<double>("MuCapture/winDCEarlyMargin");
  maxUnassignedDCHits_ = Conf.read<int>("MuCapture/maxUnassignedDCHits");

  muUVPCCellMin_ = Conf.read<int>("MuCapture/muUVPCCellMin");
  muUVPCCellMax_ = Conf.read<int>("MuCapture/muUVPCCellMax");
  muUVDCCellMin_ = Conf.read<int>("MuCapture/muUVDCCellMin");
  muUVDCCellMax_ = Conf.read<int>("MuCapture/muUVDCCellMax");

  muStopRMax_ = Conf.read<double>("MuCapture/muStopRMax");

  //       --------- Histograms initialization ---------          //
  pactCut_.init(H, Conf);
  protonWindow_.init(H, Conf);

  h_cuts_r = H.DefineTH1D("MuCapture", "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());

  h_cuts_p = H.DefineTH1D("MuCapture", "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);

  hNumPCWin_ = H.DefineTH1D( "MuCapture", "numPCWin",   "Number of PC time windows", 10, -0.5, 9.5);
  hWinPCTimeAll_ = H.DefineTH1D( "MuCapture", "winPCTimeAll",   "PC time window start, all", 1600, -6000., 10000.);
  hWinPCTimeTrig_ = H.DefineTH1D( "MuCapture", "winPCTimeTrig",   "PC time window start, trig", 200, -50., 50.);

  hWinPCTStartBeforeTrig_ = H.DefineTH1D( "MuCapture", "winPCTStartBeforeTrig",   "PC 1 win start before trigger", 1200, -6000., 0.);
  hWinPCTStartAfterTrig_ = H.DefineTH1D( "MuCapture", "winPCTStartAfter",   "PC 1 win start after trigger", 2000, 0., 10000.);

  hWinDCUnassignedEarly_ = H.DefineTH1D( "MuCapture", "winDCUnassignedEarly",   "Unassigned DC hits, early", 1000, -1000., 0.);
  hWinDCUnassignedLate_ = H.DefineTH1D( "MuCapture", "winDCUnassignedLate",   "Unassigned DC hits, late", 1000, 0., 1000.);
  hWinDCUnassignedCount_ = H.DefineTH1D( "MuCapture", "winDCUnassignedCount",   "Count of unassigned DC hits", 101, -0.5, 100.5);

  hMuRangePCFirst_ = H.DefineTH1D("MuCapture", "MuRangePCFirst",   "Muon range PC first", 13, -0.5, 12.5);
  hMuRangePCLast_ = H.DefineTH1D("MuCapture", "MuRangePCLast",   "Muon range PC last", 13, -0.5, 12.5);
  hMuRangeDCFirst_ = H.DefineTH1D("MuCapture", "MuRangeDCFirst",   "Muon range DC first", 45, -0.5, 44.5);
  hMuRangeDCLast_ = H.DefineTH1D("MuCapture", "MuRangeDCLast",   "Muon range DC last", 45, -0.5, 44.5);

  hMuUVLimitsPCUp_ = H.DefineTH2D("MuCapture", "MuUVLimitsPCUp", "Muon PC1234 max vs min coordinate", 160, 0.5, 160.5, 160, 0.5, 160.5);
  hMuUVLimitsDC_ = H.DefineTH2D("MuCapture", "MuUVLimitsDC", "Muon DC up max vs min coordinate", 80, 0.5, 80.5, 80, 0.5, 80.5);

  // Make the bin size half a cell
  hMuStopUVCell_ = H.DefineTH2D("MuCapture", "MuStopUVCell", "Muon stop V vs U position (cell units)", 107, 53.75, 107.25,  107, 53.75, 107.25);
  hMuStopUVPos_ = H.DefineTH2D("MuCapture", "MuStopUVPos", "Muon stop V vs U position (cm)", 201, -10.05, +10.05, 201, -10.05, +10.05 );
  hMuStopRadius_ = H.DefineTH1D("MuCapture", "MuStopRadius", "Muon stop R (cm)", 80, 0., 8.);

  //----------------------------------------------------------------
  hwidthPCall_.init("MuCapture/pcWidthAll", "pcwidth", 12, H, Conf);
  hwidthDCall_.init("MuCapture/dcWidthAll", "dcwidth", 44, H, Conf);

  hwidthPCMuWin_.init("MuCapture/pcWidthMuWin", "pcmuwidth", 12, H, Conf);
  hwidthDCMuWin_.init("MuCapture/dcWidthMuWin", "dcmuwidth", 44, H, Conf);

  //----------------------------------------------------------------
  hOccupancyPCAll_.init("MuCapture", "hitMapPCAll", 12, 160, H, Conf);
  hOccupancyDCAll_.init("MuCapture", "hitMapDCAll", 44, 80, H, Conf);

  hOccupancyDCUnassigned_.init("MuCapture", "hitMapDCUnassigned", 44, 80, H, Conf);

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

  return doDefaultTWIST_;
}

//================================================================
MuCapture::EventCutNumber MuCapture::analyze(EventClass &evt, HistogramFactory &hist) {

  // Fill the "all hits" histos
  hwidthPCall_.fill(evt.pc_hits_by_time());
  hwidthDCall_.fill(evt.dc_hits_by_time());

  hOccupancyPCAll_.fill(evt.pc_hits_by_time());
  hOccupancyDCAll_.fill(evt.dc_hits_by_time());

  //----------------------------------------------------------------
  // Sort PC hits into time windows

  const TDCHitWPPtrCollection filtered_pc_hits_by_time = selectHits(evt.pc_hits_by_time(), cutMinTDCWidthPC_);
  const TimeWindowCollection winpcs = constructTimeWindows(filtered_pc_hits_by_time, winPCLength_);

  const int iPCTrigWin = findTriggerWindow(winpcs);

  if(iPCTrigWin < 0) {
    return CUT_NOPCWIN;
  }

  //----------------
  // PC window histograms
  hNumPCWin_->Fill(winpcs.size());
  hWinPCTimeTrig_->Fill(winpcs[iPCTrigWin].tstart);
  for(unsigned i=0; i<winpcs.size(); ++i) {
    hWinPCTimeAll_->Fill(winpcs[i].tstart);
  }

  //----------------
  if(iPCTrigWin > 0) {
    hWinPCTStartBeforeTrig_->Fill(winpcs[iPCTrigWin - 1].tstart);
    if( std::abs(winpcs[iPCTrigWin - 1].tstart) < winPCSeparation_ ) {
      return CUT_PCWIN_TRIGSEPPAST;
    }
  }

  //----------------
  if(iPCTrigWin + 2 != winpcs.size()) {
    return CUT_PCWIN_NUMAFTERTRIG;
  }

  //----------------
  hWinPCTStartAfterTrig_->Fill(winpcs[iPCTrigWin + 1].tstart);
  if( std::abs(winpcs[iPCTrigWin + 1].tstart) < winPCSeparation_ ) {
    return CUT_PCWIN_TRIGSEPFUTURE;
  }

  //----------------
  // Process DC hits
  TDCHitWPPtrCollection unassignedDCHits;
  const TDCHitWPPtrCollection filtered_dc_hits_by_time = selectHits(evt.dc_hits_by_time(), cutMinTDCWidthDC_);
  const TimeWindowCollection windcs = assignDCHits(&unassignedDCHits, filtered_dc_hits_by_time, winpcs);
  hWinDCUnassignedCount_->Fill(unassignedDCHits.size());
  hOccupancyDCUnassigned_.fill(unassignedDCHits);

  if(unassignedDCHits.size() > maxUnassignedDCHits_) {
    return CUT_UNASSIGNEDDCHITS;
  }

  //----------------------------------------------------------------
  const ClustersByPlane muonPCClusters = constructPlaneClusters(12, winpcs[iPCTrigWin].hits);
  //std::cout<<"iPCTrigWin hits = \n"<<winpcs[iPCTrigWin].hits<<std::endl;
  //std::cout<<"muonPCClusters = \n"<<muonPCClusters<<std::endl;

  const ClustersByPlane muonDCClusters = constructPlaneClusters(44, windcs[iPCTrigWin].hits);

  const PlaneRange pcRange = findPlaneRange(muonPCClusters);
  const PlaneRange dcRange = findPlaneRange(muonDCClusters);
  if(!(pcRange.noGaps  && dcRange.noGaps)) {
    return CUT_MU_RANGE_GAPS;
  }

  hMuRangePCFirst_->Fill(pcRange.min);
  hMuRangePCLast_->Fill(pcRange.max);
  hMuRangeDCFirst_->Fill(dcRange.min);
  hMuRangeDCLast_->Fill(dcRange.max);

  if((pcRange.min != 1)||(pcRange.max != 6)) {
    return CUT_MU_PC_RANGE;
  }

  if((dcRange.min != 1)||(dcRange.max != 22)) {
    return CUT_MU_DC_RANGE;
  }

  //----------------------------------------------------------------
  // Muon UV cuts

  std::set<int> pc1234;
  pc1234.insert(1); pc1234.insert(2); pc1234.insert(3); pc1234.insert(4);
  std::pair<int,int> uvpcuplimits = uvminmax(winpcs[iPCTrigWin].hits, pc1234);
  hMuUVLimitsPCUp_->Fill(uvpcuplimits.first, uvpcuplimits.second);

  if( (uvpcuplimits.first < muUVPCCellMin_) || (muUVPCCellMax_ < uvpcuplimits.second)) {
    return CUT_MU_UV_PC;
  }

  std::pair<int,int> uvdclimits = uvminmax(windcs[iPCTrigWin].hits);
  hMuUVLimitsDC_->Fill(uvdclimits.first, uvdclimits.second);
  if((uvdclimits.first < muUVDCCellMin_) || (muUVDCCellMax_ < uvdclimits.second)) {
    return CUT_MU_UV_DC;
  }

  //----------------------------------------------------------------
  // CUT_MUSTOP_UV

  if((muonPCClusters[5].size() != 1) || (muonPCClusters[6].size() != 1)) {
    return CUT_MUSTOP_SINGLECLUSTER;
  }

  // See dt_geo.00061 and twist-coordinate-system.uvplanes.pdf
  const double muStopVCell = muonPCClusters[5].front().centralCell();
  const double muStopUCell = muonPCClusters[6].front().centralCell();
  hMuStopUVCell_->Fill(muStopUCell, muStopVCell);

  // convert to cm.
  // PC5 is V-plane rot=+135, V increases with cell number
  const double muStopV = +0.2*(muStopVCell - 80.5);
  // PC6 is U-plane rot=-135, V decreases with cell number
  const double muStopU = -0.2*(muStopUCell - 80.5);
  hMuStopUVPos_->Fill(muStopU, muStopV);

  const double muStopRadius = sqrt(std::pow(muStopV, 2) + std::pow(muStopU, 2));
  hMuStopRadius_->Fill(muStopRadius);
  if(muStopRadius > muStopRMax_) {
    return CUT_MUSTOP_UV;
  }

  //----------------------------------------------------------------
  if(1 != pactCut_.quadrant(muonPCClusters[5].front(), muonPCClusters[6].front())) {
    return CUT_MUSTOP_PACT;
  }

  //----------------------------------------------------------------
  hwidthPCMuWin_.fill(winpcs[iPCTrigWin].hits);
  hwidthDCMuWin_.fill(windcs[iPCTrigWin].hits);

  const int iProtonWin = iPCTrigWin + 1;
  protonWindow_.process(muStopU, muStopV, winpcs[iProtonWin], windcs[iProtonWin], unassignedDCHits, evt);

  //----------------------------------------------------------------
  return CUTS_ACCEPTED;
}

//================================================================
TDCHitWPPtrCollection MuCapture::selectHits(const TDCHitWPCollection& hits, double minWidthCut) {
  TDCHitWPPtrCollection res;
  for(unsigned i=0; i<hits.size(); ++i) {
    if(hits[i].width > minWidthCut) {
      res.push_back(TDCHitWPPtr(hits, i));
    }
  }
  return res;
}

//================================================================
TimeWindowCollection MuCapture::constructTimeWindows(const TDCHitWPPtrCollection& pchits, double winLength) {

  TimeWindowCollection winpcs;

  for(unsigned i=0; i< pchits.size(); ++i) {

    // This hit starts a new window
    TimeWindow win;
    win.tstart = pchits[i]->time;
    win.tend = win.tstart + winLength;

    // Put all hits falling in the  given time interval into the same window
    while((i < pchits.size()) && (pchits[i]->time < win.tend)) {
      win.hits.push_back(pchits[i]);
      ++i;
    }

    // Order hits in each window by plane/cell
    std::sort(win.hits.begin(), win.hits.end(), TDCHitWPCmpGeom());

    winpcs.push_back(win);
  }

  return winpcs;
}

//================================================================
int MuCapture::findTriggerWindow(const TimeWindowCollection& windows) {
  int itrig = -1;

  if(!windows.empty()) {
    itrig = 0;
    for(int i=1; i<windows.size(); ++i) {
      if(std::abs(windows[i].tstart) < std::abs(windows[itrig].tstart)) {
        itrig = i;
      }
    }
  }

  return itrig;
}

//================================================================
TimeWindowCollection MuCapture::assignDCHits(TDCHitWPPtrCollection *unassignedDCHits,
                                             const TDCHitWPPtrCollection& timeSortedDCHits,
                                             const TimeWindowCollection& winpcs)
{
  TimeWindowCollection windcs(winpcs.size());

  if(!timeSortedDCHits.empty()) {

    unsigned idchit=0;

    for(unsigned ipcwin=0; ipcwin<winpcs.size(); ++ipcwin) {

      TimeWindow dcwin;
      dcwin.tstart = winpcs[ipcwin].tstart - winDCEarlyMargin_;
      dcwin.tend = dcwin.tstart + winDCLength_;

      for(; idchit < timeSortedDCHits.size() && (timeSortedDCHits[idchit]->time < dcwin.tend); ++idchit) {
        TDCHitWPPtr phit = timeSortedDCHits[idchit];
        if(timeSortedDCHits[idchit]->time < dcwin.tstart) {

          unassignedDCHits->push_back(phit);

          hWinDCUnassignedEarly_->Fill(phit->time - dcwin.tstart);
          if(ipcwin > 0) {
            hWinDCUnassignedLate_->Fill(phit->time - windcs[ipcwin-1].tend);
          }

        }
        else {
          dcwin.hits.push_back(phit);
        }
      }

      // Filled DC window corresponding to the current PC window with hits
      // Order then and move to the next pc window
      std::sort(dcwin.hits.begin(), dcwin.hits.end(), TDCHitWPCmpGeom());

      windcs[ipcwin] = dcwin;

    } // for(pcwin)

    // Record any leftover hits
    for(; idchit < timeSortedDCHits.size(); ++idchit) {
      TDCHitWPPtr phit = timeSortedDCHits[idchit];
      unassignedDCHits->push_back(phit);
      hWinDCUnassignedLate_->Fill(phit->time - windcs.back().tend);
    }

  } // !timeSortedDCHits.empty()

  return windcs;
}

//================================================================
