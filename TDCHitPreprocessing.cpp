// Andrei Gaponenko, 2013

#include "TDCHitPreprocessing.h"

#include <algorithm>

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

namespace TDCHitPreprocessing {

  //================================================================
  void PassThrough::process(TDCHitWPPtrCollection *res,
                            TDCHitWPCollection *buf,
                            const TDCHitWPCollection& hits)
  {
    res->clear();
    res->reserve(hits.size());
    for(unsigned i=0; i<hits.size(); ++i) {
      res->push_back(TDCHitWPPtr(hits, i));
    }
  }

  //================================================================
  NarrowHitDiscarder::NarrowHitDiscarder(const std::string& topdir,
                                         WirePlane::DetType det,
                                         HistogramFactory& hf,
                                         const DetectorGeo& geom,
                                         const ConfigFile& conf)
    : cutMinTDCWidth_(conf.read<float>("MuCapture/TDCHitPreproc/NarrowHitDiscarder/"+WirePlane::detName(det)+"/cutMinTDCWidth"))
  {}

  //================================================================
  void NarrowHitDiscarder::process(TDCHitWPPtrCollection *res,
                                   TDCHitWPCollection * /*buf*/,
                                   const TDCHitWPCollection& hits)
  {
    res->clear();
    for(unsigned i=0; i<hits.size(); ++i) {
      if(hits[i].width() > cutMinTDCWidth_) {
        res->push_back(TDCHitWPPtr(hits, i));
      }
    }
  }

  //================================================================

} // namespace TDCHitPreprocessing
