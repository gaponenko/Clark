#include "TimeWindowingDC.h"

#include <cmath>
#include <cstdlib>
#include <limits>

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"

#include <iostream>

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

//================================================================
void TimeWindowingDC::init(HistogramFactory &hf,
                           const std::string& hdir,
                           const DetectorGeo& geom,
                           const ConfigFile &conf)
{
  winDCEnd_ = conf.read<double>("MuCapture/winDCEnd");
  winDCStart_ = conf.read<double>("MuCapture/winDCStart");
  winDCDoHistos_ = conf.read<bool>("MuCapture/winDCDoHistos");

  hWinDCMuMixStreamMap_.init(hdir, "winDCMuMixStreamMap", 44, 80, hf, conf);
  hWinDCProtonMixStreamMap_.init(hdir, "winDCProtonMixStreamMap", 44, 80, hf, conf);

  if(winDCDoHistos_) {
    hOccupancyDCUnassigned_.init(hdir, "hitMapDCUnassigned", 44, 80, hf, conf);
    hWinDCUnassignedEarly_ = hf.DefineTH1D(hdir, "winDCUnassignedEarly", "Unassigned DC hits, early", 1000, -1000., 0.);
    hWinDCUnassignedLate_ = hf.DefineTH1D(hdir, "winDCUnassignedLate", "Unassigned DC hits, late", 1000, 0., 1000.);
  }
}

//================================================================
void TimeWindowingDC::assignDCHits(const TDCHitWPPtrCollection& dcHits,
                                   TimeWindowingResults *inout)
{
  TimeWindowCollection& windows = inout->windows;
  for(unsigned i=0; i < windows.size(); ++i) {
    windows[i].tstartDC = windows[i].tstart + winDCStart_;
    windows[i].tendDC = windows[i].tstart + winDCEnd_;
  }

  //----------------
  for(unsigned ihit=0; ihit < dcHits.size(); ++ihit) {
    TDCHitWPPtr phit = dcHits[ihit];

    const TimeWindow::StreamType dcStream = (phit->plane() <= 22) ?
      TimeWindow::UPSTREAM : TimeWindow::DOWNSTREAM;

    // Look for a time-compatible PC "analysis" (not before trigger)
    // window of the same stream type first
    bool assigned = false;
    for(unsigned iwin = inout->iTrigWin; !assigned && (iwin < windows.size()); ++iwin) {
      TimeWindow& win = windows[iwin];
      if(win.stream == dcStream) {
        if((win.tstartDC < phit->time()) && (phit->time() <= win.tendDC)) {
          windows[iwin].dcHits.push_back(phit);
          assigned = true;
        }
      }
    }

    // Look for any time-compatible PC window here
    for(unsigned iwin = 0; !assigned && (iwin < windows.size()); ++iwin) {
      TimeWindow& win = windows[iwin];
      if((win.tstartDC < phit->time()) && (phit->time() <= win.tendDC)) {

        if(iwin == inout->iTrigWin) {
          if(win.stream != TimeWindow::MIXED) {
            hWinDCMuMixStreamMap_.fill(*phit);
          }
        }
        else if(iwin == 1 + inout->iTrigWin) {
          if(win.stream != TimeWindow::MIXED) {
            hWinDCProtonMixStreamMap_.fill(*phit);
          }
        }

        windows[iwin].dcHits.push_back(phit);
        win.stream = TimeWindow::MIXED;
        assigned = true;
      }
    }

    if(!assigned) {
      inout->unassignedDCHits.push_back(phit);
    }
  }

  if(winDCDoHistos_) {
    fillDiagnostics(*inout);
  }
}

//================================================================
void TimeWindowingDC::fillDiagnostics(const TimeWindowingResults& wres) {

  hOccupancyDCUnassigned_.fill(wres.unassignedDCHits);

  for(unsigned ihit = 0; ihit < wres.unassignedDCHits.size(); ++ihit) {

    TDCHitWPPtr phit = wres.unassignedDCHits[ihit];

    double dtMinEarly = std::numeric_limits<double>::max();
    double dtMinLate = std::numeric_limits<double>::max();

    for(unsigned iwin = 0; iwin < wres.windows.size(); ++iwin) {
      if(phit->time() < wres.windows[iwin].tstart) {
        dtMinEarly = std::min(dtMinEarly, wres.windows[iwin].tstartDC - phit->time());
      }
      else {
        dtMinLate = std::min(dtMinLate, phit->time() -  wres.windows[iwin].tendDC);
      }
    }

    if(dtMinEarly < dtMinLate) {
      hWinDCUnassignedEarly_->Fill(-dtMinEarly);
    }
    else {
      hWinDCUnassignedLate_->Fill(dtMinLate);
    }

  } // for (unassigned hits)
}

//================================================================
