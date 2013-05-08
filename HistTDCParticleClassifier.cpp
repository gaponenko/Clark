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

  hLastVsRestMaxWidth_ = hf.DefineTH2D(hdir, "LastVsRestMaxWidthDC","LastVsRestMaxWidthDC",
                                       200, 0., 2000., 200, 0., 2000.);
  hLastVsRestMaxWidth_->SetOption("colz");

  hLastVsRestMeanWidth_ = hf.DefineTH2D(hdir, "LastVsRestMeanWidthDC","LastVsRestMeanWidthDC",
                                        200, 0., 2000., 200, 0., 2000.);
  hLastVsRestMeanWidth_->SetOption("colz");

  hLastVsRestMedianWidth_ = hf.DefineTH2D(hdir, "LastVsRestMedianWidthDC","LastVsRestMedianWidthDC",
                                          200, 0., 2000., 200, 0., 2000.);
  hLastVsRestMedianWidth_->SetOption("colz");
}

//================================================================
void HistTDCParticleClassifier::fill(const ClustersByPlane& gc) {
  TDCHitStats stat7, stat8;
  stat7.fill(gc[29]);
  stat8.fill(gc[30]);

  hpc8vs7maxWidth_->Fill(stat7.widthStats().max(), stat8.widthStats().max());
  hpc8vs7meanWidth_->Fill(stat7.widthStats().mean(), stat8.widthStats().mean());
  hpc8vs7medianWidth_->Fill(stat7.widthStats().median(), stat8.widthStats().median());
  hpc8vs7maxHits_->Fill(stat7.maxHitsPerWire(), stat8.maxHitsPerWire());

  const PlaneRange gr = findPlaneRange(gc);
  if(gr.max >= 32) { // at least 2 DC planes hit
    TDCHitStats last, rest;
    last.fill(gc[gr.max]);
    for(int i=31; i<gr.max; ++i) {
      rest.fill(gc[i]);
    }

    hLastVsRestMaxWidth_->Fill(rest.widthStats().max(), last.widthStats().max());
    hLastVsRestMeanWidth_->Fill(rest.widthStats().mean(), last.widthStats().mean());
    hLastVsRestMedianWidth_->Fill(rest.widthStats().median(), last.widthStats().median());
  }
}

//================================================================
