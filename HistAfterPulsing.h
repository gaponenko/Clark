// Andrei Gaponenko, 2013

#ifndef HistAfterPulsing_h
#define HistAfterPulsing_h

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
class HistAfterPulsing {
public:
  void init(const std::string& hdir,
            unsigned maxPlaneNumber,
            unsigned maxCellNumber,
            HistogramFactory &hf,
            const ConfigFile &conf);

  // pass by value because the method needs a copy anyway
  void fill(TDCHitWPPtrCollection hits);
  HistAfterPulsing() : hcellHitMultiplicity_() {}

private :
  std::vector<TH1*> hSameCellDt_;
  std::vector<TH2*> hSameCellTDCWidthVsWidth_;
  std::vector<TH2*> hSameCellTDCWidthVsDt_;
  HistOccupancy hOccupancyMultiHit_;
  TProfile2D *hcellHitMultiplicity_;
};

#endif/*HistAfterPulsing_h*/
