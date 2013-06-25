// Andrei Gaponenko, 2013

#include "HistWinTime.h"

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistWinTime::init(const std::string& hdir,
                       const std::string& hname,
                       HistogramFactory &hf,
                       const ConfigFile& conf)
{
  hWinTime_ = hf.DefineTH1D(hdir,
                            (hname.empty()? "winTime" : hname),
                            "window time" + (hname.empty() ? "" : ", "+hname),
                            1600, -6000., 10000.);
}

//================================================================
void HistWinTime::fill(const TimeWindowingResults& wres) {
  for(unsigned iwin = 0; iwin < wres.windows.size(); ++iwin) {
    hWinTime_->Fill(wres.windows[iwin].tstart);
  }
}

//================================================================
