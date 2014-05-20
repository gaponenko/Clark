// Unlike HistXtalk we use a pre-defined cut to identify "wide" hits.
//
// Andrei Gaponenko, 2014

#ifndef HistXT2_h
#define HistXT2_h

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
class HistXT2 {
public:
  // Note that setting cutNeighborDistanceMin=cutNeighborDistanceMax=0
  // will look at hits on the same wire, but in a different way than
  // HistAfterPulsing: in Xtalk we look at all combinations of hits
  // not just at the pairs of hits closest in time to each other.
  void init(const std::string& hdir,
            unsigned maxPlaneNumber,
            unsigned cutNeighborDistanceMin,
            unsigned cutNeighborDistanceMax,
            double cutHitWidth,
            HistogramFactory &hf,
            const ConfigFile &conf);

  // pass by value because the method needs a copy anyway
  void fill(TDCHitWPPtrCollection hits);

  HistXT2() {}

private :
  unsigned cutNeighborDistanceMin_;
  unsigned cutNeighborDistanceMax_;
  double cutHitWidth_;

  std::vector<TH1*> hWireDistance_;
  std::vector<TH1*> hdt_;
  std::vector<TH1*> hdw_;
  std::vector<TH2*> hdtvsdw_;
};

#endif/*HistXT2_h*/
