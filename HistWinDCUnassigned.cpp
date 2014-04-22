// Andrei Gaponenko, 2014

#include "HistWinDCUnassigned.h"

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistWinDCUnassigned::init(const std::string& hdir,
                               HistogramFactory &hf,
                               const ConfigFile& conf)
{
  hWinDCUnassignedAll_ = hf.DefineTH2D(hdir,
                                       "unassignedAll",
                                       "Num time windows vs per-stream unassigned DC hits",
                                       201, -100.5, +100.5,
                                       6, -0.5, 5.5
                                       );

  hWinDCUnassignedAll_->SetOption("colz");

  hWinDCUnassignedAfterTrig_ = hf.DefineTH2D(hdir,
                                             "unassignedAfterTrig",
                                             "Num after-trig time windows vs per-stream unassigned DC after-trig hits",
                                             201, -100.5, +100.5,
                                             6, -0.5, 5.5
                                             );

  hWinDCUnassignedAfterTrig_->SetOption("colz");

  hWinDCNumUnassignedAfterTrig_ = hf.DefineTH1D(hdir,
                                                "numUnassignedAfterTrig",
                                                "Num unassigned DC after-trig hits",
                                                51, -0.5, +50.5);
}

//================================================================
void HistWinDCUnassigned::fill(const TimeWindowingResults& wres) {

  int unassignedAllUp=0, unassignedAllDn=0;
  for(unsigned i = 0; i < wres.unassignedDCHits.size(); ++i) {
    ++( (wres.unassignedDCHits[i]->plane() < 23 ) ? unassignedAllUp : unassignedAllDn);
  }

  if(unassignedAllUp) hWinDCUnassignedAll_->Fill(-unassignedAllUp, wres.windows.size());
  if(unassignedAllDn) hWinDCUnassignedAll_->Fill(+unassignedAllDn, wres.windows.size());

  if(wres.iTrigWin != -1) {
    const double trigTime = wres.windows[wres.iTrigWin].tstart;
    int unassignedAfterTrigUp=0, unassignedAfterTrigDn=0;
    for(unsigned i = 0; i < wres.unassignedDCHits.size(); ++i) {
      if(wres.unassignedDCHits[i]->time() > trigTime) {
        ++( (wres.unassignedDCHits[i]->plane() < 23 ) ? unassignedAfterTrigUp : unassignedAfterTrigDn);
      }
    }
    const int numAfterTrigWindows =  wres.windows.size() - wres.iTrigWin - 1;
    if(unassignedAfterTrigUp) hWinDCUnassignedAfterTrig_->Fill(-unassignedAfterTrigUp, numAfterTrigWindows);
    if(unassignedAfterTrigDn) hWinDCUnassignedAfterTrig_->Fill(+unassignedAfterTrigDn, numAfterTrigWindows);
    hWinDCNumUnassignedAfterTrig_->Fill(unassignedAfterTrigDn + unassignedAfterTrigUp);
  }
}

//================================================================
