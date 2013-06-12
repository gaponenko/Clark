// Andrei Gaponenko, 2013

#include "HistAccidentals.h"

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistAccidentals::init(const std::string& hdir,
                           HistogramFactory &hf,
                           const ConfigFile& conf)
{
  cutPreTrigTimeMax_ = conf.read<double>("MuCapture/Accidentals/tmax");
  const double cycleLength = conf.read<double>("MuCapture/Accidentals/cycleLength");
  const double numCycles = conf.read<double>("MuCapture/Accidentals/numCycles");
  cutPreTrigTimeMin_ = cutPreTrigTimeMax_ - numCycles * cycleLength;

  htstartAll_ = hf.DefineTH1D(hdir, "tstartAll", "t, all windows",
                              8*numCycles, cutPreTrigTimeMin_, cutPreTrigTimeMax_);

  htstartDn_ = hf.DefineTH1D(hdir, "tstartDn", "t, Downstream windows",
                             8*numCycles, cutPreTrigTimeMin_, cutPreTrigTimeMax_);

  hnumwinAll_ = hf.DefineTH1D(hdir, "numWinAll", "number of all windows in range",11, -0.5, 10.5);
  hnumwinDn_ = hf.DefineTH1D(hdir, "numWinDn", "number of downstream windows in range",11, -0.5, 10.5);
}

//================================================================
void HistAccidentals::fill(const TimeWindowingResults& wres) {
  int numAll(0), numDn(0);
  for(unsigned iwin = 0; iwin < wres.iTrigWin; ++iwin) {
    htstartAll_->Fill(wres.windows[iwin].tstart);
    if(wres.windows[iwin].stream == TimeWindow::DOWNSTREAM) {
      htstartDn_->Fill(wres.windows[iwin].tstart);
    }

    if((cutPreTrigTimeMin_ < wres.windows[iwin].tstart ) && (wres.windows[iwin].tstart <= cutPreTrigTimeMax_)) {

      ++numAll;

      if(wres.windows[iwin].stream == TimeWindow::DOWNSTREAM) {
        ++numDn;
      }
    }
  }

  hnumwinAll_->Fill(numAll);
  hnumwinDn_->Fill(numDn);
}

//================================================================
