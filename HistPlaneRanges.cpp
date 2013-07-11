// Andrei Gaponenko, 2013

#include "HistPlaneRanges.h"
#include <cassert>

#include "PlaneRange.h"
#include "TimeWindow.h"

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

//================================================================
void HistPlaneRanges::init(const std::string& hdir,
                           HistogramFactory &hf,
                           const DetectorGeo& geom,
                           const ConfigFile& conf)
{
  geom_ = &geom;

  //----------------------------------------------------------------
  hNumRanges_ = hf.DefineTH1D(hdir, "numRanges", "num ranges", 10, -0.5, 9.5);

  hSingleRange_ = hf.DefineTH2D(hdir, "singleRange", "End vs begin for single range",
                                57, -0.5, 56.5, 57, -0.5, 56.5);
  hSingleRange_->SetOption("colz");

  hDoubleRangeGap_ = hf.DefineTH2D(hdir, "doubleRangeGap", "Double range gap",
                                   57, -0.5, 56.5, 57, -0.5, 56.5);
  hDoubleRangeGap_->SetOption("colz");

  hDoubleRangeMissingPlanes_ = hf.DefineTH1D(hdir, "doubleRangeMissingPlanes", "doubleRangeMissingPlanes", 57, -0.5, 56.5);

  hDoubleRangeSizes_ = hf.DefineTH2D(hdir, "doubleRangeSizes", "doubleRangeSizes",
                                     57, -0.5, 56.5, 57, -0.5, 56.5);

  hDoubleRangeSizes_->SetOption("colz");
}

//================================================================
void HistPlaneRanges::fill(const PlaneRange& gr) {

  hNumRanges_->Fill(gr.segments().size());
  switch(gr.segments().size()) {
  case 1:
    hSingleRange_->Fill(gr.min(), gr.max());
    break;
  case 2:
    hDoubleRangeGap_->Fill(gr.segments()[0].max, gr.segments()[1].min);
    assert(gr.segments()[0].max + 1 < gr.segments()[1].min);
    for(int i = gr.segments()[0].max + 1; i < gr.segments()[1].min; ++i) {
      hDoubleRangeMissingPlanes_->Fill(i);
    }

    hDoubleRangeSizes_->Fill(1 + gr.segments()[0].max - gr.segments()[0].min,
                             1 + gr.segments()[1].max - gr.segments()[1].min);

    break;
  default:
    break;
  }
}

//================================================================
