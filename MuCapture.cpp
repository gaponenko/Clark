// Andrei Gaponenko, 2013

#include "MuCapture.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TH1.h"
#include "TH2.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

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
  TimeWindowCollection windcs;
  windcs.reserve(winpcs.size());

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

      windcs.push_back(dcwin);

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
