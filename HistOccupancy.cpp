// Andrei Gaponenko, 2013

#include "HistOccupancy.h"

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
void HistOccupancy::init(const std::string& hdir,
                         const std::string& namePrefix,
                         unsigned maxPlaneNumber,
                         unsigned maxCellNumber,
                         HistogramFactory &hf,
                         const ConfigFile& conf)
{
  hitMap_ = hf.DefineTH2D(hdir, namePrefix, namePrefix,
                          1+maxPlaneNumber, -0.5, maxPlaneNumber+0.5,
                          1+maxCellNumber, -0.5, maxCellNumber+0.5);

  hitMap_->SetOption("colz");
}

//================================================================
void HistOccupancy::fill(const TDCHitWP& hit) {
  hitMap_->Fill(hit.plane(), hit.cell());
}

//================================================================
void HistOccupancy::fill(const TDCHitWPCollection& hits) {
  for(TDCHitWPCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(*i);
  }
}

//================================================================
void HistOccupancy::fill(const TDCHitWPPtrCollection& hits) {
  for(TDCHitWPPtrCollection::const_iterator i = hits.begin(); i!=hits.end(); ++i) {
    fill(**i);
  }
}

//================================================================
