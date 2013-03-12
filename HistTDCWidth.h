// Andrei Gaponenko, 2013

#ifndef HistTDCWidth_h
#define HistTDCWidth_h

#include <vector>

#include "TDCHitWP.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;

//================================================================
class HistTDCWidth {
public:
  void init(const std::string& hdir,
            const std::string& namePrefix,
            unsigned maxPlaneNumber,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const TDCHitWP& hit);
  void fill(const TDCHitWPCollection& hits);
  void fill(const TDCHitWPPtrCollection& hits);

  HistTDCWidth() : hWidth_() {}

private :
  TH1 *hWidth_;
  std::vector<TH1*> byPlane_;
};

#endif/*HistTDCWidth_h*/
