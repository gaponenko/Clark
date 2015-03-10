// Combinations of variables to try for the contained channel.
//
// Andrei Gaponenko, 2015

#ifndef MuCapContainedVars_h
#define MuCapContainedVars_h

#include <string>

#include "WireCluster.h"

#include "MuCapTrkContainmentCut.h"

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class EventClass;


namespace MuCapContainedVars {

  struct Result {
    double xvar;
    double yvar;
    bool contained;
    Result(double x, double y, bool c) : xvar(x), yvar(y), contained(c) {}
  };

  class IVarProcessor {
  public:

    // return ref to make for easier use of c_str() when calling ROOT API.
    virtual const std::string& xtitle() const = 0;
    virtual const std::string& ytitle() const = 0;

    virtual void init(const std::string& hdir,
                      HistogramFactory &hf,
                      const DetectorGeo&,
                      const ConfigFile& conf) = 0;


    virtual Result compute(const EventClass& evt, int iPosTrack, const ClustersByPlane& globalPlaneClusters) = 0;

    virtual ~IVarProcessor() {}
  };

  //================================================================
  class RangeCosVsP: virtual public IVarProcessor {
    static std::string xtitle_;
    static std::string ytitle_;
    MuCapTrkContainmentCut ccut_;
  public:
    virtual const std::string& xtitle() const { return xtitle_; }
    virtual const std::string& ytitle() const { return ytitle_; }

    virtual void init(const std::string& hdir,
                      HistogramFactory &hf,
                      const DetectorGeo&,
                      const ConfigFile& conf);

    virtual Result compute(const EventClass& evt, int iPosTrack, const ClustersByPlane& globalPlaneClusters);
  };

  //================================================================
}


#endif/*MuCapContainedVars_h*/
