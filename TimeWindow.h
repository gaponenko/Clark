// Andrei Gaponenko, 2013

#ifndef TimeWindow_h
#define TimeWindow_h

#include <vector>
#include "TDCHitWP.h"

struct TimeWindow {
  enum StreamType { UPSTREAM=-1, MIXED=0, DOWNSTREAM=+1 };

  StreamType stream;
  double tstart; // PC
  double tendPC;
  double tstartDC;
  double tendDC;

  TDCHitWPPtrCollection pcHits;
  TDCHitWPPtrCollection dcHits;

  TimeWindow() : stream(MIXED), tstart(), tendPC(), tstartDC(), tendDC() {}
};


typedef std::vector<TimeWindow> TimeWindowCollection;

struct TimeWindowingResults {
  TimeWindowCollection windows;
  unsigned iTrigWin;
  TDCHitWPPtrCollection unassignedDCHits;
  TimeWindowingResults(): iTrigWin(-1u) {}
};



#endif/*TimeWindow_h*/
