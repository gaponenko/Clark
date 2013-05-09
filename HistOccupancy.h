// Andrei Gaponenko, 2013

#ifndef HistOccupancy_h
#define HistOccupancy_h

#include <vector>

#include "TDCHitWP.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;

//================================================================
class HistOccupancy {
public:
  void init(const std::string& hdir,
            const std::string& namePrefix,
            unsigned maxPlaneNumber,
            unsigned maxCellNumber,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const TDCHitWP& hit);
  void fill(const TDCHitWPCollection& hits);
  void fill(const TDCHitWPPtrCollection& hits);

  HistOccupancy() : hitMap_() {}

private :
  TH2 *hitMap_;
};

#endif/*HistOccupancy_h*/
