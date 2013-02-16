// Andrei Gaponenko, 2013

#include "MuCapture.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TH1.h"
#include "TH2.h"

#include "TDCHitWP.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//================================================================
struct TimeWindow {
  double tstart;
  double tend;
  std::vector<TDCHitWPPtr> hits;
  TimeWindow() : tstart(), tend() {}
};

typedef std::vector<TimeWindow> TimeWindowCollection;

//================================================================
bool MuCapture::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *Log) {
  if (not Conf.read<bool>("MuCapture/Do")) {
    Log->info( "MuCapture code will not be run");
    return false ;
  }
  Log->info( "Register MuCapture module");

  //       --------- Histograms initialization ---------          //
  hNumPCWin_ = H.DefineTH1D( "MuCapture", "numPCWin",   "Number of PC time windows", 10, -0.5, 9.5);
  hWinPCTimeAll_ = H.DefineTH1D( "MuCapture", "winPCTimeAll",   "PC time window start, all", 1600, -6000., 10000.);
  hWinPCTimeTrig_ = H.DefineTH1D( "MuCapture", "winPCTimeTrig",   "PC time window start, trig", 200, -50., 50.);

  //       --------- Parameters initialization ---------          //
  doDefaultTWIST_ = Conf.read<bool>("MuCapture/doDefaultTWIST");
  windtpc_ = Conf.read<double>("MuCapture/win_dt_pc");

  return true;
}

bool MuCapture::Process(EventClass &evt, HistogramFactory &hist) {

  TimeWindowCollection winpcs;
  const TDCHitWPCollection& pchits = evt.pc_hits_by_time();
  for(unsigned i=0; i< pchits.size(); ++i) {

    // This hit starts a new window
    TimeWindow win;
    win.tstart = pchits[i].time;
    win.tend = win.tstart + windtpc_;

    // Put all hits falling in the windtpc_ time interval into the same window
    while((i < pchits.size()) && (pchits[i].time < win.tend)) {
      win.hits.push_back(TDCHitWPPtr(pchits, i));
      ++i;
    }

    winpcs.push_back(win);
  }

  hNumPCWin_->Fill(winpcs.size());

  int iPCTrigWin = -1;
  if(!winpcs.empty()) {
    double tpctrig = winpcs[0].tstart;
    iPCTrigWin = 0;
    for(unsigned i=0; i<winpcs.size(); ++i) {
      hWinPCTimeAll_->Fill(winpcs[i].tstart);
      if(std::abs(winpcs[i].tstart) < std::abs(tpctrig)) {
        tpctrig = winpcs[i].tstart;
        iPCTrigWin = i;
      }
    }

    hWinPCTimeTrig_->Fill(tpctrig);
  }

  return doDefaultTWIST_;
}
