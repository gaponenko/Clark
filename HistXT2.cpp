// Andrei Gaponenko, 2013

#include "HistXT2.h"

#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "TH1.h"
#include "TProfile2D.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistXT2::init(const std::string& hdir,
                   unsigned maxPlaneNumber,
                   unsigned cutNeighborDistanceMin,
                   unsigned cutNeighborDistanceMax,
                   double cutHitWidth,
                   HistogramFactory &hf,
                   const ConfigFile& conf)
{
  cutNeighborDistanceMin_ = cutNeighborDistanceMin;
  cutNeighborDistanceMax_ = cutNeighborDistanceMax;
  cutHitWidth_ = cutHitWidth;

  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"wireDist_"<<std::setw(2)<<std::setfill('0')<<plane;
    hWireDistance_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 81, -0.5, 80.5));
  }

  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"dt_"<<std::setw(2)<<std::setfill('0')<<plane;
    hdt_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 300, -150., 150.));
  }

  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"dwidth_"<<std::setw(2)<<std::setfill('0')<<plane;
    hdw_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 200, 0., 1000.));
  }

  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"dtvsdw_"<<std::setw(2)<<std::setfill('0')<<plane;
    hdtvsdw_.push_back(hf.DefineTH2D(hdir, os.str(), "dt vs dwidth", 200, 0., 1000., 400, -1000., 1000.));
    hdtvsdw_.back()->SetOption("colz");
  }
}

//================================================================
namespace {
  struct PlaneCmp {
    bool operator()(const TDCHitWPPtr& h1, const TDCHitWPPtr& h2) const {
      return h1->plane() < h2->plane();
    }
  };
}

void HistXT2::fill(TDCHitWPPtrCollection hits) {
  std::sort(hits.begin(), hits.end(), TDCHitWPCmpGeom());

  typedef TDCHitWPPtrCollection::iterator Iter;
  typedef std::pair<Iter,Iter> IterRange;

  Iter current = hits.begin();
  while(current != hits.end()) {
    const int currentPlane = (*current)->plane();
    IterRange rr = std::equal_range(current, hits.end(), *current, PlaneCmp());
    for(Iter ih1 = rr.first; ih1 != rr.second; ++ih1) {
      for(Iter ih2 = ih1 + 1; ih2 != rr.second; ++ih2) {
        const int dist = std::abs((*ih1)->cell() - (*ih2)->cell());
        if((cutNeighborDistanceMin_ <= dist) && (dist <= cutNeighborDistanceMax_)) {

          // the first hit is the wider one:
          const bool ordered = ((*ih1)->width() > (*ih2)->width());
          const TDCHitWP& hit1 = ordered ? **ih1 : **ih2;
          const TDCHitWP& hit2 = ordered ? **ih2 : **ih1;

          if(hit1.width() > cutHitWidth_) {
            hWireDistance_[currentPlane-1]->Fill(dist);
            hdt_[currentPlane-1]->Fill(hit2.time() - hit1.time());
            hdw_[currentPlane-1]->Fill(hit1.width() - hit2.width());
            hdtvsdw_[currentPlane-1]->Fill(hit1.width() - hit2.width(), hit2.time() - hit1.time());
          }
        }
      }
    }
    current = rr.second;
  }
}

//================================================================
