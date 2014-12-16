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

  void fill(const TDCHitWPCollection& hits);
  void fill(const TDCHitWPPtrCollection& hits);

  HistTDCWidth() : hWidth_() {}

private :
  TH1 *hWidth_;
  std::vector<TH1*> byPlane_;

  TH1 *hminWidth_;
  TH1 *hmaxWidth_;
  TH1 *hmeanWidth_;
  TH1 *hmedianWidth_;
  TH1 *hhitsPerSet_;

  std::vector<TH1*> minWidth_;
  std::vector<TH1*> maxWidth_;
  std::vector<TH1*> meanWidth_;
  std::vector<TH1*> medianWidth_;
  std::vector<TH1*> hitsPerSet_;

  void fill(const TDCHitWP& hit);
};

#endif/*HistTDCWidth_h*/
