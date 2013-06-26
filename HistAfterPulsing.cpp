// Andrei Gaponenko, 2013

#include "HistAfterPulsing.h"

#include <sstream>
#include <vector>
#include <algorithm>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistAfterPulsing::init(const std::string& hdir,
                            unsigned maxPlaneNumber,
                            unsigned maxCellNumber,
                            HistogramFactory &hf,
                            const ConfigFile& conf)
{
  for(int plane = 0; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"numHitsPerCell"<<std::setw(2)<<std::setfill('0')<<plane;
    hNumHitsPerCell_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 10, -0.5, 9.5));
  }


  for(int plane = 0; plane <= maxPlaneNumber; ++plane) {
    std::ostringstream os;
    os<<"sameCellDt"<<std::setw(2)<<std::setfill('0')<<plane;
    hSameCellDt_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 1600, 0., 16000.));
  }

  hOccupancyMultiHit_.init(hdir, "mutiHitCells", maxPlaneNumber, maxCellNumber, hf, conf);
}

//================================================================
void HistAfterPulsing::fill(TDCHitWPPtrCollection hits) {
  std::sort(hits.begin(), hits.end(), TDCHitWPCmpTime());
  std::stable_sort(hits.begin(), hits.end(), TDCHitWPCmpGeom());
  // here we have hits sorted by cell then time
  //std::cout<<"hits = "<<hits<<std::endl;

  int extraHitCount = 0;
  for(int ihit = 0; ihit + 1 < hits.size(); ++ihit) {
    if(hits[ihit]->cid() == hits[ihit+1]->cid()) {
      ++extraHitCount;
      //std::cout<<"hits = "<<hits<<std::endl;
      hOccupancyMultiHit_.fill(hits[ihit]->cid());
      hSameCellDt_[hits[ihit]->plane()]->Fill(hits[ihit+1]->time() - hits[ihit]->time());
    }
    else {
      hNumHitsPerCell_[hits[ihit]->plane()]->Fill(extraHitCount);
      extraHitCount = 0;
    }
  }

  if(!hits.empty()) {
    hNumHitsPerCell_[hits.back()->plane()]->Fill(extraHitCount);
  }
}

//================================================================
