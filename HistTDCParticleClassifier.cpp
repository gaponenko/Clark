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
                                     const ConfigFile &conf,
                                     TimeWindow::StreamType stream)
{
  geom_ = &geom;
  stream_ = stream;

  hpc8vs7maxWidth_ = hf.DefineTH2D(hdir, "tgtpc2vs1maxWidth","Target PC 2nd vs 1st max width", 200, 0., 2000., 200, 0., 2000.);
  hpc8vs7maxWidth_->SetOption("colz");

  hpc8vs7meanWidth_ = hf.DefineTH2D(hdir, "tgtpc2vs1meanWidth","Target PC 2nd vs 1st mean width", 200, 0., 2000., 200, 0., 2000.);
  hpc8vs7meanWidth_->SetOption("colz");

  hpc8vs7medianWidth_ = hf.DefineTH2D(hdir, "tgtpc2vs1medianWidth","Target PC 2nd vs 1st median width", 200, 0., 2000., 200, 0., 2000.);
  hpc8vs7medianWidth_->SetOption("colz");

  hpc8vs7maxHits_ = hf.DefineTH2D(hdir, "tgtpc2vs1maxHits","tgtpc2vs1maxHits per wire", 10, -0.5, 9.5, 10, -0.5, 9.5);
  hpc8vs7maxHits_->SetOption("colz");

  htgtpcmaxWidth_ = hf.DefineTH1D(hdir, "tgtpcmaxWidth","Target PC max width", 200, 0., 2000.);
  htgtpcmeanWidth_ = hf.DefineTH1D(hdir, "tgtpcmeanWidth","Target PC mean width", 200, 0., 2000.);
  htgtpcmedianWidth_ = hf.DefineTH1D(hdir, "tgtpcmedianWidth","Target PC median width", 200, 0., 2000.);

  hMaxWidthDC_ = hf.DefineTH1D(hdir, "maxWidthDC", "maxWidthDC", 200, 0., 2000);
  hMeanWidthDC_ = hf.DefineTH1D(hdir, "meanWidthDC", "meanWidthDC", 200, 0., 2000);
  hMedianWidthDC_ = hf.DefineTH1D(hdir, "medianWidthDC", "medianWidthDC", 200, 0., 2000);

  hMaxWidthPC_ = hf.DefineTH1D(hdir, "maxWidthPC", "maxWidthPC", 200, 0., 2000);
  hMeanWidthPC_ = hf.DefineTH1D(hdir, "meanWidthPC", "meanWidthPC", 200, 0., 2000);
  hMedianWidthPC_ = hf.DefineTH1D(hdir, "medianWidthPC", "medianWidthPC", 200, 0., 2000);
}

//================================================================
void HistTDCParticleClassifier::fill(const ClustersByPlane& gc) {

  TDCHitStats tgt1, tgt2, tgt12; // PCs closest to the target and the next one
  switch(stream_) {
  case TimeWindow::DOWNSTREAM:
    tgt1.fill(gc[29]);
    tgt2.fill(gc[30]);
    tgt12.fill(gc[29]); tgt12.fill(gc[30]);
    break;
  case TimeWindow::UPSTREAM:
    tgt1.fill(gc[28]);
    tgt2.fill(gc[27]);
    tgt12.fill(gc[28]); tgt12.fill(gc[27]);
    break;
  default:
    throw std::runtime_error("HistTDCParticleClassifier does not support MIXED windows");
  }

  if(tgt1.widthStats().numEntries() && tgt2.widthStats().numEntries()) {
    hpc8vs7maxWidth_->Fill(tgt1.widthStats().max(), tgt2.widthStats().max());
    hpc8vs7meanWidth_->Fill(tgt1.widthStats().mean(), tgt2.widthStats().mean());
    hpc8vs7medianWidth_->Fill(tgt1.widthStats().median(), tgt2.widthStats().median());
    hpc8vs7maxHits_->Fill(tgt1.maxHitsPerWire(), tgt2.maxHitsPerWire());
  }
  if(tgt12.widthStats().numEntries()) {
    htgtpcmaxWidth_->Fill(tgt12.widthStats().max());
    htgtpcmeanWidth_->Fill(tgt12.widthStats().mean());
    htgtpcmedianWidth_->Fill(tgt12.widthStats().median());
  }

  TDCHitStats dcstats;
  TDCHitStats pcstats;
  const unsigned first = (stream_ == TimeWindow::DOWNSTREAM) ? 1+geom_->numGlobal()/2 : 1;
  const unsigned last  = (stream_ == TimeWindow::DOWNSTREAM) ? geom_->numGlobal() : geom_->numGlobal()/2;

  for(unsigned i=first; i<=last; ++i) {
    ((geom_->global(i).planeType() == WirePlane::PC) ? pcstats : dcstats).fill(gc[i]);
  }

  if(dcstats.widthStats().numEntries()) { // have DC hits
    hMaxWidthDC_->Fill(dcstats.widthStats().max());
    hMeanWidthDC_->Fill(dcstats.widthStats().mean());
    hMedianWidthDC_->Fill(dcstats.widthStats().median());
  }

  if(pcstats.widthStats().numEntries()) { // have PC hits
    hMaxWidthPC_->Fill(pcstats.widthStats().max());
    hMeanWidthPC_->Fill(pcstats.widthStats().mean());
    hMedianWidthPC_->Fill(pcstats.widthStats().median());
  }
}

//================================================================
