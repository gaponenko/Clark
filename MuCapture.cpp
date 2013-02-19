// Andrei Gaponenko, 2013

#include "MuCapture.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "TH1.h"
#include "TH2.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//================================================================
namespace { // local helpers
  struct PlaneRange {
    int min;
    int max;
    bool noGaps;
    PlaneRange() : min(-1), max(-1), noGaps(false) {}
  };

  PlaneRange findRange(const ClustersByPlane& cp) {
    PlaneRange res;
    typedef ClustersByPlane::const_iterator Iter;

    Iter i = cp.begin();
    if(i != cp.end()) {
      assert(!i->second.empty());
      res.min = i->first;
      res.max = i->first;

      for(++i; i!=cp.end(); ++i) {
        if(res.min > i->first) {
          res.min = i->first;
        }
        if(i->first > res.max) {
          res.max= i->first;
        }
      }
    }

    res.noGaps = ((res.max - res.min + 1) == cp.size());

    return res;
  }
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
  winPCLength_ = Conf.read<double>("MuCapture/winPCLength");
  winPCSeparation_ = Conf.read<double>("MuCapture/winPCSeparation");
  winDCLength_ = Conf.read<double>("MuCapture/winDCLength");
  winDCEarlyMargin_ = Conf.read<double>("MuCapture/winDCEarlyMargin");
  maxUnassignedDCHits_ = Conf.read<double>("MuCapture/maxUnassignedDCHits");

  //       --------- Histograms initialization ---------          //

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

  //----------------------------------------------------------------
  // Sort PC hits into time windows

  const TimeWindowCollection winpcs = constructTimeWindows(evt.pc_hits_by_time(), winPCLength_);

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
  const TimeWindowCollection windcs = assignDCHits(&unassignedDCHits, evt.dc_hits_by_time(), winpcs);
  hWinDCUnassignedCount_->Fill(unassignedDCHits.size());

  if(unassignedDCHits.size() > maxUnassignedDCHits_) {
    return CUT_UNASSIGNEDDCHITS;
  }

  //----------------
  const ClustersByPlane muonPCClusters = constructPlaneClusters(winpcs[iPCTrigWin].hits);
  const ClustersByPlane muonDCClusters = constructPlaneClusters(windcs[iPCTrigWin].hits);

  const PlaneRange pcRange = findRange(muonPCClusters);
  const PlaneRange dcRange = findRange(muonDCClusters);
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

  //----------------

  return CUTS_ACCEPTED;
}

//================================================================
TimeWindowCollection MuCapture::constructTimeWindows(const TDCHitWPCollection& pchits, double winLength) {

  TimeWindowCollection winpcs;

  for(unsigned i=0; i< pchits.size(); ++i) {

    // This hit starts a new window
    TimeWindow win;
    win.tstart = pchits[i].time;
    win.tend = win.tstart + winLength;

    // Put all hits falling in the  given time interval into the same window
    while((i < pchits.size()) && (pchits[i].time < win.tend)) {
      win.hits.push_back(TDCHitWPPtr(pchits, i));
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
                                             const TDCHitWPCollection& timeSortedDCHits,
                                             const TimeWindowCollection& winpcs)
{
  TimeWindowCollection windcs(winpcs.size());

  if(!timeSortedDCHits.empty()) {

    unsigned idchit=0;

    for(unsigned ipcwin=0; ipcwin<winpcs.size(); ++ipcwin) {

      TimeWindow dcwin;
      dcwin.tstart = winpcs[ipcwin].tstart - winDCEarlyMargin_;
      dcwin.tend = dcwin.tstart + winDCLength_;

      for(; idchit < timeSortedDCHits.size() && (timeSortedDCHits[idchit].time < dcwin.tend); ++idchit) {
        TDCHitWPPtr phit(timeSortedDCHits, idchit);
        if(timeSortedDCHits[idchit].time < dcwin.tstart) {

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
      TDCHitWPPtr phit(timeSortedDCHits, idchit);
      unassignedDCHits->push_back(phit);
      hWinDCUnassignedLate_->Fill(phit->time - windcs.back().tend);
    }

  } // !timeSortedDCHits.empty()

  return windcs;
}

//================================================================
ClustersByPlane MuCapture::constructPlaneClusters(const TDCHitWPPtrCollection& hits) {
  ClustersByPlane res;
  for(unsigned i=0; i<hits.size(); ++i) {

    WireCluster cl;
    cl.plane = hits[i]->plane;
    cl.hits.push_back(hits[i]);

    for(; i<hits.size();++i) {
      if(cl.plane != hits[i]->plane) break;
      if(cl.hits.back()->cell + 1 != hits[i]->cell) break;
      cl.hits.push_back(hits[i]);
    }

    res[cl.plane].push_back(cl);
  }

  return res;
}

//================================================================
