// Andrei Gaponenko, 2013

#include "HistHitStream.h"
#include <cassert>

#include "PlaneRange.h"
#include "TimeWindow.h"

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

//================================================================
void HistHitStream::init(const std::string& hdir,
                         HistogramFactory &hf,
                         const DetectorGeo& geom,
                         const ConfigFile& conf)
{
  geom_ = &geom;

  //----------------
  hNumPlanesUpVsDn_ = hf.DefineTH2D(hdir, "planesUpVsDn", "Number of hit planes upstream vs downstream",
                                    29, -0.5, 28.5, 29, -0.5, 28.5);
  hNumPlanesUpVsDn_->SetOption("colz");

  hNumPlanesDnMinusUp_ = hf.DefineTH1D(hdir, "planesDnMinusUp", "Dn - Up hit planes", 57, -28.5, 28.5);

  //----------------
  hNumPCsUpVsDn_ = hf.DefineTH2D(hdir, "pcsUpVsDn", "Number of PC planes upstream vs downstream",
                                 7, -0.5, 6.5, 7, -0.5, 6.5);
  hNumPCsUpVsDn_->SetOption("colz");

  hNumPCsDnMinusUp_ = hf.DefineTH1D(hdir, "pcsDnMinusUp", "Dn - Up PC planes", 13, -6.5, 6.5);

  //----------------
  hNumDCsUpVsDn_ = hf.DefineTH2D(hdir, "dcsUpVsDn", "Number of DC planes upstream vs downstream",
                                 23, -0.5, 22.5, 23, -0.5, 22.5);
  hNumDCsUpVsDn_->SetOption("colz");

  hNumDCsDnMinusUp_ = hf.DefineTH1D(hdir, "dcsDnMinusUp", "Dn - Up DC planes", 45, -22.5, 22.5);

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

  //----------------------------------------------------------------
  hWinStream_ = hf.DefineTH1D(hdir, "winStream", "Time window stream type", 3, -1.5, +1.5);
}
//================================================================
void HistHitStream::fill(const ClustersByPlane& gc) {

  // Count the number of planes hit up and dn of the target
  int numPlanesUp(0), numPlanesDn(0);
  int numPCsUp(0), numPCsDn(0);
  int numDCsUp(0), numDCsDn(0);

  // Using the 1-based TWIST convention for plane numbers
  for(unsigned i=1; i<=geom_->numGlobal(); ++i) {
    if(!gc[i].empty()) {
      ++(i <= geom_->numGlobal()/2 ? numPlanesUp : numPlanesDn);
      if(geom_->global(i).planeType() == WirePlane::PC) {
        ++(i <= geom_->numGlobal()/2 ? numPCsUp : numPCsDn);
      }
      if(geom_->global(i).planeType() == WirePlane::DC) {
        ++(i <= geom_->numGlobal()/2 ? numDCsUp : numDCsDn);
      }
    }
  }

  hNumPlanesUpVsDn_->Fill(numPlanesDn, numPlanesUp);
  if(numPlanesUp + numPlanesDn) {
    hNumPlanesDnMinusUp_->Fill(numPlanesDn - numPlanesUp);
  }

  hNumPCsUpVsDn_->Fill(numPCsDn, numPCsUp);
  if(numPCsUp + numPCsDn) {
    hNumPCsDnMinusUp_->Fill(numPCsDn - numPCsUp);
  }

  hNumDCsUpVsDn_->Fill(numDCsDn, numDCsUp);
  if(numDCsUp + numDCsDn) {
    hNumDCsDnMinusUp_->Fill(numDCsDn - numDCsUp);
  }

  //----------------------------------------------------------------
  PlaneRange gr = findPlaneRange(gc);
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

  //----------------------------------------------------------------
  const TimeWindow::StreamType stream =
    (gr.max() <= 28) ? TimeWindow::UPSTREAM :
    (29 <= gr.min()) ? TimeWindow::DOWNSTREAM :
    TimeWindow::MIXED;

  hWinStream_->Fill(stream);
}

//================================================================
