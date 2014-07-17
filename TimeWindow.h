// Andrei Gaponenko, 2013

#ifndef TimeWindow_h
#define TimeWindow_h

#include <vector>
#include <ostream>

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

std::ostream& operator<<(std::ostream& os, const TimeWindow& win);
std::ostream& operator<<(std::ostream& os, const TimeWindowingResults& wres);

#endif/*TimeWindow_h*/
