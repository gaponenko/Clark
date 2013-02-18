// Andrei Gaponenko, 2013

#ifndef TimeWindow_h
#define TimeWindow_h

#include <vector>
#include "TDCHitWP.h"

struct TimeWindow {
  double tstart;
  double tend;
  std::vector<TDCHitWPPtr> hits;
  TimeWindow() : tstart(), tend() {}
};

typedef std::vector<TimeWindow> TimeWindowCollection;

#endif/*TimeWindow_h*/
