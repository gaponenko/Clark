// Andrei Gaponenko, 2013

#ifndef HistXtalk_h
#define HistXtalk_h

#include <string>
#include <vector>

#include "TDCHitWP.h"
#include "HistOccupancy.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class ConfigFile;


//================================================================
class HistXtalk {
public:
  // Note that setting cutNeighborDistanceMin=cutNeighborDistanceMax=0
  // will look at hits on the same wire, but in a different way than
  // HistAfterPulsing: in Xtalk we look at all combinations of hits
  // not just at the pairs of hits closest in time to each other.
  void init(const std::string& hdir,
            unsigned maxPlaneNumber,
            unsigned cutNeighborDistanceMin,
            unsigned cutNeighborDistanceMax,
            HistogramFactory &hf,
            const ConfigFile &conf);

  // pass by value because the method needs a copy anyway
  void fill(TDCHitWPPtrCollection hits);

  HistXtalk() {}

private :
  unsigned cutNeighborDistanceMin_;
  unsigned cutNeighborDistanceMax_;
  std::vector<TH2*> hTDCWidthVsWidth_;
  std::vector<TH2*> hTDCWidthVsDt_;
};

#endif/*HistXtalk_h*/
