// Andrei Gaponenko, 2013

#include "HistDriftTime.h"

#include <sstream>
#include <vector>
#include <utility>

#include "TH1.h"
#include "TProfile2D.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistDriftTime::init(const std::string& hdir,
                         HistogramFactory &hf,
                         unsigned maxPlaneNumber,
                         double driftTimeHistLimit,
                         const ConfigFile& conf)
{
  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"driftTime"<<std::setw(2)<<std::setfill('0')<<plane;
    hdriftTime_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 200, 0., driftTimeHistLimit));
  }
}

//================================================================
void HistDriftTime::fill(const EventClass& evt, int idio, const TDCHitWPPtrCollection& hits) {
  if(idio >= 0) { // Only look at events with a good DIO track

    // The selected DIO tracks are supposed to start at the stopping target.
    // Look at the planes that should be crossed by the found track.
    std::pair<int,int> planeRange = (evt.costh[idio] < 0.) ?
      std::make_pair(1UL, 1 + hdriftTime_.size()/2) :
      std::make_pair(1 + hdriftTime_.size()/2, 1 + hdriftTime_.size());

    const double t0 = evt.hefit_time[idio];
    for(unsigned i=0; i<hits.size(); ++i) {
      if((planeRange.first <= hits[i]->plane()) && (hits[i]->plane() < planeRange.second)) {
        hdriftTime_.at(hits[i]->plane() - 1)->Fill(hits[i]->time() - t0);
      }
    }
  }
}

//================================================================
