// Andrei Gaponenko, 2013

#include "HistTDCWidth.h"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
namespace {
  class MinWidthFinder {
  public:
    typedef std::map<int,double> MinWidthMap;
    const MinWidthMap& getMinWidths() const { return mw_; }
    void addHit(const TDCHitWP& hit);
  private:
    MinWidthMap mw_;
  };

  void MinWidthFinder::addHit(const TDCHitWP& hit) {
    MinWidthMap::iterator p = mw_.find(hit.plane);
    if(p == mw_.end()) {
      mw_.insert(std::make_pair(hit.plane, double(hit.width)));
    }
    else {
      p->second = std::min(double(hit.width), p->second);
    }
  }
}

//================================================================
void HistTDCWidth::init(const std::string& hdir,
                        const std::string& namePrefix,
                        unsigned maxPlaneNumber,
                        HistogramFactory &hf,
                        const ConfigFile& conf)
{
  hWidth_ = hf.DefineTH1D(hdir, namePrefix, namePrefix + " TDC width", 1000, 0., 1000.);

  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC Width, plane "<<std::setw(2)<<std::setfill('0')<<i;
    byPlane_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 1000, 0., 1000.));
  }

  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_minWidth_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC MIN Width in plane "<<std::setw(2)<<std::setfill('0')<<i;
    minWidth_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 1000, 0., 1000.));
  }
}

//================================================================
void HistTDCWidth::fill(const TDCHitWP& hit) {
  hWidth_->Fill(hit.width);
  byPlane_[hit.plane]->Fill(hit.width);
}

//================================================================
void HistTDCWidth::fill(const TDCHitWPCollection& hits) {
  MinWidthFinder mw;
  for(TDCHitWPCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(*i);
    mw.addHit(*i);
  }

  typedef MinWidthFinder::MinWidthMap PM;
  const PM& mm = mw.getMinWidths();
  for(PM::const_iterator i=mm.begin(); i!=mm.end(); ++i) {
    minWidth_[i->first]->Fill(i->second);
  }
}

//================================================================
void HistTDCWidth::fill(const TDCHitWPPtrCollection& hits) {
  MinWidthFinder mw;
  for(TDCHitWPPtrCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(**i);
    mw.addHit(**i);
  }

  typedef MinWidthFinder::MinWidthMap PM;
  const PM& mm = mw.getMinWidths();
  for(PM::const_iterator i=mm.begin(); i!=mm.end(); ++i) {
    minWidth_[i->first]->Fill(i->second);
  }
}

//================================================================
