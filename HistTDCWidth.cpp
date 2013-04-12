// Andrei Gaponenko, 2013

#include "HistTDCWidth.h"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
namespace {
  struct Extrema {
    double min;
    double max;

    Extrema()
      : min(std::numeric_limits<double>::max())
      , max(std::numeric_limits<double>::min())
    {}

    void fill(double width) {
      min = std::min(min, width);
      max = std::max(max, width);
    }
  };

  class ExtremeWidthFinder {
  public:
    typedef std::map<int,Extrema> WidthMap;
    const WidthMap& getWidths() const { return mw_; }
    void addHit(const TDCHitWP& hit);
  private:
    WidthMap mw_;
  };

  void ExtremeWidthFinder::addHit(const TDCHitWP& hit) {
    mw_[hit.plane].fill(hit.width);
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

  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_maxWidth_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC MAX Width in plane "<<std::setw(2)<<std::setfill('0')<<i;
    maxWidth_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 1000, 0., 1000.));
  }
}

//================================================================
void HistTDCWidth::fill(const TDCHitWP& hit) {
  hWidth_->Fill(hit.width);
  byPlane_[hit.plane]->Fill(hit.width);
}

//================================================================
void HistTDCWidth::fill(const TDCHitWPCollection& hits) {
  ExtremeWidthFinder mw;
  for(TDCHitWPCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(*i);
    mw.addHit(*i);
  }

  typedef ExtremeWidthFinder::WidthMap PM;
  const PM& mm = mw.getWidths();
  for(PM::const_iterator i=mm.begin(); i!=mm.end(); ++i) {
    minWidth_[i->first]->Fill(i->second.min);
    maxWidth_[i->first]->Fill(i->second.max);
  }
}

//================================================================
void HistTDCWidth::fill(const TDCHitWPPtrCollection& hits) {
  ExtremeWidthFinder mw;
  for(TDCHitWPPtrCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(**i);
    mw.addHit(**i);
  }

  typedef ExtremeWidthFinder::WidthMap PM;
  const PM& mm = mw.getWidths();
  for(PM::const_iterator i=mm.begin(); i!=mm.end(); ++i) {
    minWidth_[i->first]->Fill(i->second.min);
    maxWidth_[i->first]->Fill(i->second.max);
  }
}

//================================================================