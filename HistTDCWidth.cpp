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
#include "MuCapUtilities.h"

//================================================================
namespace {
  class ExtremeWidthFinder {
  public:
    typedef std::map<int,MuCapUtilities::Stats> WidthMap;
    const WidthMap& getWidths() const { return mw_; }
    void addHit(const TDCHitWP& hit);
  private:
    WidthMap mw_;
  };

  void ExtremeWidthFinder::addHit(const TDCHitWP& hit) {
    mw_[hit.plane()].fill(hit.width());
  }
}

//================================================================
void HistTDCWidth::init(const std::string& hdir,
                        const std::string& namePrefix,
                        unsigned maxPlaneNumber,
                        HistogramFactory &hf,
                        const ConfigFile& conf)
{
  hWidth_ = hf.DefineTH1D(hdir, namePrefix, namePrefix + " TDC width", 2000, 0., 2000.);

  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC Width, plane "<<std::setw(2)<<std::setfill('0')<<i;
    byPlane_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 2000, 0., 2000.));
    byPlane_.back()->GetXaxis()->SetTitle("TDC width [ns]");
  }

  {
    std::string hn(namePrefix+"_minWidth");
    std::string ht(namePrefix+"TDC MIN Width");
    hminWidth_ = hf.DefineTH1D(hdir, hn.c_str(), ht.c_str(), 2000, 0., 2000.);
  }
  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_minWidth_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC MIN Width in plane "<<std::setw(2)<<std::setfill('0')<<i;
    minWidth_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 2000, 0., 2000.));
    minWidth_.back()->GetXaxis()->SetTitle("TDC width [ns]");
  }

  {
    std::string hn(namePrefix+"_maxWidth");
    std::string ht(namePrefix+"TDC MAX Width");
    hmaxWidth_ = hf.DefineTH1D(hdir, hn.c_str(), ht.c_str(), 2000, 0., 2000.);
  }
  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_maxWidth_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC MAX Width in plane "<<std::setw(2)<<std::setfill('0')<<i;
    maxWidth_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 2000, 0., 2000.));
    maxWidth_.back()->GetXaxis()->SetTitle("TDC width [ns]");
  }

  {
    std::string hn(namePrefix+"_meanWidth");
    std::string ht(namePrefix+"TDC mean width");
    hmeanWidth_ = hf.DefineTH1D(hdir, hn.c_str(), ht.c_str(), 2000, 0., 2000.);
  }
  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_meanWidth_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC mean width in plane "<<std::setw(2)<<std::setfill('0')<<i;
    meanWidth_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 2000, 0., 2000.));
    meanWidth_.back()->GetXaxis()->SetTitle("TDC width [ns]");
  }

  {
    std::string hn(namePrefix+"_medianWidth");
    std::string ht(namePrefix+"TDC median width");
    hmedianWidth_ = hf.DefineTH1D(hdir, hn.c_str(), ht.c_str(), 2000, 0., 2000.);
  }
  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_medianWidth_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" TDC median width in plane "<<std::setw(2)<<std::setfill('0')<<i;
    medianWidth_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 2000, 0., 2000.));
    medianWidth_.back()->GetXaxis()->SetTitle("TDC width [ns]");
  }

  {
    std::string hn(namePrefix+"_hitsPerSet");
    std::string ht(namePrefix+" Num hites per set");
    hhitsPerSet_ = hf.DefineTH1D(hdir, hn.c_str(), ht.c_str(), 10, -0.5, 9.5);
  }
  for(unsigned i=0; i<=maxPlaneNumber; ++i) {
    std::ostringstream osname;
    osname<<namePrefix<<"_hitsPerSet_"<<std::setw(2)<<std::setfill('0')<<i;
    std::ostringstream ostitle;
    ostitle<<namePrefix<<" Num hits per set in plane "<<std::setw(2)<<std::setfill('0')<<i;
    hitsPerSet_.push_back(hf.DefineTH1D(hdir, osname.str(), ostitle.str(), 10, -0.5, 9.5));
    hitsPerSet_.back()->GetXaxis()->SetTitle("Hit multiplicity");
  }
}

//================================================================
void HistTDCWidth::fill(const TDCHitWP& hit) {
  hWidth_->Fill(hit.width());
  byPlane_[hit.plane()]->Fill(hit.width());
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
    minWidth_[i->first]->Fill(i->second.min());
    maxWidth_[i->first]->Fill(i->second.max());
    meanWidth_[i->first]->Fill(i->second.mean());
    medianWidth_[i->first]->Fill(i->second.median());
    hitsPerSet_[i->first]->Fill(i->second.numEntries());

    hminWidth_->Fill(i->second.min());
    hmaxWidth_->Fill(i->second.max());
    hmeanWidth_->Fill(i->second.mean());
    hmedianWidth_->Fill(i->second.median());
    hhitsPerSet_->Fill(i->second.numEntries());
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
    minWidth_[i->first]->Fill(i->second.min());
    maxWidth_[i->first]->Fill(i->second.max());
    meanWidth_[i->first]->Fill(i->second.mean());
    medianWidth_[i->first]->Fill(i->second.median());
    hitsPerSet_[i->first]->Fill(i->second.numEntries());

    hminWidth_->Fill(i->second.min());
    hmaxWidth_->Fill(i->second.max());
    hmeanWidth_->Fill(i->second.mean());
    hmedianWidth_->Fill(i->second.median());
    hhitsPerSet_->Fill(i->second.numEntries());
  }
}

//================================================================
