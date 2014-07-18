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
    // res should point to valid empty collections at the call time.
    // After the call *res will contain the result.
    virtual void process(TDCHitWPPtrCollection *res,
                         const TDCHitWPCollection& inputs) = 0;

    virtual ~IProcessor() {}
  };

  //----------------------------------------------------------------
  class PassThrough : public IProcessor {
  public:
    virtual void process(TDCHitWPPtrCollection *res,
                         const TDCHitWPCollection& inputs);
  };

  //----------------------------------------------------------------
  // drop hits that MOFIA flagged as cross talk
  class MOFIA_XTalkDiscarder : public IProcessor {
  public:
    MOFIA_XTalkDiscarder(const std::string& topdir,
                         WirePlane::DetType det,
                         HistogramFactory& hf,
                         const DetectorGeo& geom,
                         const ConfigFile& conf);

    virtual void process(TDCHitWPPtrCollection *res,
                         const TDCHitWPCollection& inputs);
  private:
    TH1 *hxt_;
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
                         const TDCHitWPCollection& inputs);

  private:
    float cutMinTDCWidth_;
  };

  //----------------------------------------------------------------
  class SameWireHitDiscarder : public IProcessor {
  public:
    SameWireHitDiscarder(const std::string& topdir,
                         WirePlane::DetType det,
                         HistogramFactory& hf,
                         const DetectorGeo& geom,
                         const ConfigFile& conf);

    virtual void process(TDCHitWPPtrCollection *res,
                         const TDCHitWPCollection& inputs);

  private:
    float cutSameWireDt_;
    std::vector<TH1*> hSameCellDtAll_;
    std::vector<TH1*> hSameCellDtDropped_;
    std::vector<TH1*> hSameCellWidthDropped_;
    std::vector<TH1*> hSameCellWidthKept_;
    HistOccupancy hSameCellOccupancyDropped_;
    HistOccupancy hSameCellOccupancyKept_;
  };

  //================================================================
  // The "data" class

  class Hits {
  public:
    const TDCHitWPPtrCollection& get() const { return phits_; }

    Hits(const TDCHitWPCollection& in, IProcessor& proc) {
      proc.process(&phits_, in);
    }

  private:
    TDCHitWPPtrCollection phits_;
  };

  //================================================================
}

#endif/*TDCHitPreprocessing_h*/
