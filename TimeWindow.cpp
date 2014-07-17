#include "TimeWindow.h"

#include <algorithm>

std::ostream& operator<<(std::ostream& os, const TimeWindow& win) {
  os<<"TimeWindow: "<<"tstart = "<<win.tstart<<", stream="<<win.stream
    <<"\n      pcHits = "<<win.pcHits
    <<"\n      dcHits = "<<win.dcHits
    <<std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const TimeWindowingResults& wres) {

  os<<"TimeWindowingResults: numWin="<<wres.windows.size()<<", iTrigWin="<<wres.iTrigWin<<std::endl;

  for(unsigned i=0; i<wres.windows.size(); ++i) {
    os<<"  "<<i<<" "<<wres.windows[i]<<std::endl;
  }

  os<<"   unassignedDCHits = "<<wres.unassignedDCHits<<std::endl;

  return os;
}
