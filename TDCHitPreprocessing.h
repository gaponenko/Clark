// Andrei Gaponenko, 2013

#ifndef TDCHitPreprocessing_h
#define TDCHitPreprocessing_h

#include <string>
#include <vector>

#include "TDCHitWP.h"
#include "DetectorGeo.h"

#include "HistOccupancy.h"

class HistogramFactory;
class ConfigFile;

class TH1;
class TH2;

namespace  TDCHitPreprocessing {

  //================================================================
  // Interface for a "processor" class

  class IProcessor {
  public:
    // res and buf arguments should point to valid empty collections.
    // After the call *res will contain the result.  It may or may not
    // use *buf for storage.
    virtual void process(TDCHitWPPtrCollection *res,
                         TDCHitWPCollection *buf,
                         const TDCHitWPCollection& inputs) = 0;

    virtual ~IProcessor() {}
  };

  //----------------------------------------------------------------
  class PassThrough : public IProcessor {
  public:
    virtual void process(TDCHitWPPtrCollection *res,
                         TDCHitWPCollection *buf,
                         const TDCHitWPCollection& inputs);
  };

  //----------------------------------------------------------------
  class NarrowHitDiscarder : public IProcessor {
  public:
    NarrowHitDiscarder(const std::string& topdir,
                       WirePlane::DetType det,
                       HistogramFactory& hf,
                       const DetectorGeo& geom,
                       const ConfigFile& conf);

    virtual void process(TDCHitWPPtrCollection *res,
                         TDCHitWPCollection *buf,
                         const TDCHitWPCollection& inputs);

  private:
    float cutMinTDCWidth_;
  };

  //================================================================
  // The "data" class

  class Hits {
  public:
    const TDCHitWPPtrCollection& get() const { return phits_; }

    Hits(const TDCHitWPCollection& in, IProcessor& proc) {
      proc.process(&phits_, &buf_, in);
    }

  private:
    TDCHitWPCollection buf_;
    TDCHitWPPtrCollection phits_;
  };

  //================================================================
}

#endif/*TDCHitPreprocessing_h*/
