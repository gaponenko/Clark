// Andrei Gaponenko, 2013

#include "HistXtalk.h"

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
void HistXtalk::init(const std::string& hdir,
                     unsigned maxPlaneNumber,
                     unsigned cutNeighborDistanceMin,
                     unsigned cutNeighborDistanceMax,
                     HistogramFactory &hf,
                     const ConfigFile& conf)
{
  cutNeighborDistanceMin_ = cutNeighborDistanceMin;
  cutNeighborDistanceMax_ = cutNeighborDistanceMax;

  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"TDCWidth2VsWidth1_"<<std::setw(2)<<std::setfill('0')<<plane;
    hTDCWidthVsWidth_.push_back(hf.DefineTH2D(hdir, os.str(), os.str(), 200, 0., 1000., 200, 0., 1000.));
    hTDCWidthVsWidth_.back()->SetOption("colz");
  }

  for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"TDCWidth2VsDt"<<std::setw(2)<<std::setfill('0')<<plane;
    hTDCWidthVsDt_.push_back(hf.DefineTH2D(hdir, os.str(), os.str(), 200, 0., 1000., 200, 0., 1000.));
    hTDCWidthVsDt_.back()->SetOption("colz");
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

void HistXtalk::fill(TDCHitWPPtrCollection hits) {
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
          const bool ordered = ((*ih1)->time() < (*ih2)->time());
          const TDCHitWP& hit1 = ordered ? **ih1 : **ih2;
          const TDCHitWP& hit2 = ordered ? **ih2 : **ih1;
          hTDCWidthVsWidth_[currentPlane-1]->Fill(hit1.width(), hit2.width());
          hTDCWidthVsDt_[currentPlane-1]->Fill(hit2.time() - hit1.time(), hit2.width());
        }
      }
    }
    current = rr.second;
  }
}

//================================================================
