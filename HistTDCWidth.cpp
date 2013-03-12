// Andrei Gaponenko, 2013

#include "HistTDCWidth.h"

#include <algorithm>
#include <sstream>
#include <iomanip>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

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
}

//================================================================
void HistTDCWidth::fill(const TDCHitWP& hit) {
  hWidth_->Fill(hit.width);
  byPlane_[hit.plane]->Fill(hit.width);
}

//================================================================
void HistTDCWidth::fill(const TDCHitWPCollection& hits) {
  for(TDCHitWPCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(*i);
  }
}

//================================================================
void HistTDCWidth::fill(const TDCHitWPPtrCollection& hits) {
  for(TDCHitWPPtrCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(**i);
  }
}

//================================================================
