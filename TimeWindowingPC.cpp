#include "TimeWindowingPC.h"

#include <cmath>
#include <cstdlib>

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"

//================================================================
void TimeWindowingPC::init(HistogramFactory &hf,
                           const std::string& hdir,
                           const DetectorGeo& geom,
                           const ConfigFile &conf)
{
  winPCLength_ = conf.read<double>("MuCapture/winPCLength");
  winTrigMaxdt_ = conf.read<double>("MuCapture/winTrigMaxdt");

  hNumPCWin_ = hf.DefineTH1D(hdir, "numPCWin",   "Number of PC time windows", 10, -0.5, 9.5);
  hWinPCTimeAll_ = hf.DefineTH1D(hdir, "winPCTimeAll",   "PC time window start, all", 1600, -6000., 10000.);

  hWinPCTimeTrig_ = hf.DefineTH1D(hdir, "winPCTimeTrig",   "PC time window start, trig", 200, -50., 50.);

  hWinPCTStartBeforeTrig_ = hf.DefineTH1D(hdir, "winPCTStartBeforeTrig",   "PC 1 win start before trigger", 1200, -6000., 0.);
  hWinPCTStartAfterTrig_ = hf.DefineTH1D(hdir, "winPCTStartAfter",   "PC 1 win start after trigger", 2000, 0., 10000.);

  hWinTimeTypes_ = hf.DefineTH2D(hdir, "winTimeTypes", "iwin-itrig vs PC win type", 3, -1.5, 1.5, 20, -9.5, 10.5);
  hWinTimeTypes_->SetOption("colz");

  hWinMultTypes_ = hf.DefineTH2D(hdir, "winMultTypes", "Num PC win vs win type", 3, -1.5, 1.5, 20, -0.5, 19.5);
  hWinMultTypes_->SetOption("colz");
}

//================================================================
void TimeWindowingPC::assignPCHits(const TDCHitWPPtrCollection& pchits, TimeWindowingResults *out) {
  TimeWindowCollection winpcs;

  for(unsigned i=0; i < pchits.size(); ++i) {

    // This hit starts a new window
    TimeWindow win;
    win.tstart = pchits[i]->time();
    win.tendPC = win.tstart + winPCLength_;

    // Put all hits falling in the  given time interval into the same window
    while((i < pchits.size()) && (pchits[i]->time() < win.tendPC)) {
      win.pcHits.push_back(pchits[i]);
      ++i;
    }

    // Order hits in each window by plane/cell
    std::sort(win.pcHits.begin(), win.pcHits.end(), TDCHitWPCmpGeom());

    // figure out window stream type
    const int pmin = win.pcHits.front()->plane();
    const int pmax = win.pcHits.back()->plane();
    win.stream = (pmax <= 6) ? TimeWindow::UPSTREAM :
      ((7 <= pmin) ? TimeWindow::DOWNSTREAM : TimeWindow::MIXED);

    winpcs.push_back(win);
  }

  const unsigned iPCTrigWin = findTriggerWindow(winpcs);

  //----------------------------------------------------------------
  // PC window histograms
  hNumPCWin_->Fill(winpcs.size());
  for(unsigned i=0; i<winpcs.size(); ++i) {
    hWinPCTimeAll_->Fill(winpcs[i].tstart);
  }

  // Histos that require trigger and more to be defined
  if(iPCTrigWin != -1u) {
    hWinPCTimeTrig_->Fill(winpcs[iPCTrigWin].tstart);
    if(0 < iPCTrigWin) {
      hWinPCTStartBeforeTrig_->Fill(winpcs[iPCTrigWin - 1].tstart);
    }
    if(iPCTrigWin + 1 < winpcs.size()) {
      hWinPCTStartAfterTrig_->Fill(winpcs[iPCTrigWin + 1].tstart);
    }

    for(unsigned i=0; i<winpcs.size(); ++i) {
      hWinTimeTypes_->Fill(winpcs[i].stream, int(i - iPCTrigWin));
      hWinMultTypes_->Fill(winpcs[i].stream, winpcs.size());
    }
  }

  //----------------------------------------------------------------
  out->windows.swap(winpcs);
  out->iTrigWin = iPCTrigWin;
}

//================================================================
unsigned TimeWindowingPC::findTriggerWindow(const TimeWindowCollection& windows) {
  unsigned itrig = -1u;

  if(!windows.empty()) {
    itrig = 0;
    for(int i=1; i<windows.size(); ++i) {
      if(std::abs(windows[i].tstart) < std::abs(windows[itrig].tstart)) {
        itrig = i;
      }
    }

    // If the closes window is too far from t=0 don't call it the trigger
    if(std::abs(windows[itrig].tstart) > winTrigMaxdt_) {
      itrig = -1u;
    }
  }

  return itrig;
}

//================================================================
