// Andrei Gaponenko, 2013

#include "HistTDCParticleClassifier.h"
#include "TDCHitStats.h"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "PlaneRange.h"

//================================================================
void HistTDCParticleClassifier::init(HistogramFactory &hf,
                                     const std::string& hdir,
                                     const DetectorGeo& geom,
                                     const ConfigFile &conf)
{
  hpc8vs7maxWidth_ = hf.DefineTH2D(hdir, "pc8vs7maxWidth","pc8vs7maxWidth", 200, 0., 2000., 200, 0., 2000.);
  hpc8vs7maxWidth_->SetOption("colz");

  hpc8vs7meanWidth_ = hf.DefineTH2D(hdir, "pc8vs7meanWidth","pc8vs7meanWidth", 200, 0., 2000., 200, 0., 2000.);
  hpc8vs7meanWidth_->SetOption("colz");

  hpc8vs7medianWidth_ = hf.DefineTH2D(hdir, "pc8vs7medianWidth","pc8vs7medianWidth", 200, 0., 2000., 200, 0., 2000.);
  hpc8vs7medianWidth_->SetOption("colz");

  hpc8vs7maxHits_ = hf.DefineTH2D(hdir, "pc8vs7maxHits","pc8vs7maxHits", 10, -0.5, 9.5, 10, -0.5, 9.5);
  hpc8vs7maxHits_->SetOption("colz");

  hLastVsRestMaxWidthDC_ = hf.DefineTH2D(hdir, "LastVsRestMaxWidthDC","LastVsRestMaxWidthDC",
                                       200, 0., 2000., 200, 0., 2000.);
  hLastVsRestMaxWidthDC_->SetOption("colz");

  hLastVsRestMeanWidthDC_ = hf.DefineTH2D(hdir, "LastVsRestMeanWidthDC","LastVsRestMeanWidthDC",
                                        200, 0., 2000., 200, 0., 2000.);
  hLastVsRestMeanWidthDC_->SetOption("colz");

  hLastVsRestMedianWidthDC_ = hf.DefineTH2D(hdir, "LastVsRestMedianWidthDC","LastVsRestMedianWidthDC",
                                          200, 0., 2000., 200, 0., 2000.);
  hLastVsRestMedianWidthDC_->SetOption("colz");

  hMaxWidthDC_ = hf.DefineTH1D(hdir, "maxWidthDC", "maxWidthDC", 200, 0., 2000);
  hMeanWidthDC_ = hf.DefineTH1D(hdir, "meanWidthDC", "meanWidthDC", 200, 0., 2000);
  hMedianWidthDC_ = hf.DefineTH1D(hdir, "medianWidthDC", "medianWidthDC", 200, 0., 2000);

  hMaxWidthPC_ = hf.DefineTH1D(hdir, "maxWidthPC", "maxWidthPC", 200, 0., 2000);
  hMeanWidthPC_ = hf.DefineTH1D(hdir, "meanWidthPC", "meanWidthPC", 200, 0., 2000);
  hMedianWidthPC_ = hf.DefineTH1D(hdir, "medianWidthPC", "medianWidthPC", 200, 0., 2000);
}

//================================================================
void HistTDCParticleClassifier::fill(const ClustersByPlane& gc) {
  const PlaneRange gr = findPlaneRange(gc);

  if(gr.max >= 30) { // PC7 and 8
    TDCHitStats stat7, stat8;
    stat7.fill(gc[29]);
    stat8.fill(gc[30]);

    hpc8vs7maxWidth_->Fill(stat7.widthStats().max(), stat8.widthStats().max());
    hpc8vs7meanWidth_->Fill(stat7.widthStats().mean(), stat8.widthStats().mean());
    hpc8vs7medianWidth_->Fill(stat7.widthStats().median(), stat8.widthStats().median());
    hpc8vs7maxHits_->Fill(stat7.maxHitsPerWire(), stat8.maxHitsPerWire());
  }

  if(gr.max >= 32) { // at least 2 DC planes hit
    TDCHitStats last, rest;
    last.fill(gc[gr.max]);
    for(int i=31; i<gr.max; ++i) {
      rest.fill(gc[i]);
    }

    hLastVsRestMaxWidthDC_->Fill(rest.widthStats().max(), last.widthStats().max());
    hLastVsRestMeanWidthDC_->Fill(rest.widthStats().mean(), last.widthStats().mean());
    hLastVsRestMedianWidthDC_->Fill(rest.widthStats().median(), last.widthStats().median());
  }

  if(gr.max >= 31) { // have DC hits
    TDCHitStats dcstats;
    for(int i=31; i<=gr.max; ++i) {
      dcstats.fill(gc[i]);
    }
    hMaxWidthDC_->Fill(dcstats.widthStats().max());
    hMeanWidthDC_->Fill(dcstats.widthStats().mean());
    hMedianWidthDC_->Fill(dcstats.widthStats().median());
  }

  if(gr.max >= 29) { // have PC hits
    TDCHitStats pcstats;
    for(int i=29; i<=30; ++i) {
      pcstats.fill(gc[i]);
    }

    hMaxWidthPC_->Fill(pcstats.widthStats().max());
    hMeanWidthPC_->Fill(pcstats.widthStats().mean());
    hMedianWidthPC_->Fill(pcstats.widthStats().median());
  }
}

//================================================================
