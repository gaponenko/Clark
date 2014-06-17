// Andrei Gaponenko, 2013

#include "TDCHitPreprocessing.h"

#include <algorithm>

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

namespace TDCHitPreprocessing {

  //================================================================
  void PassThrough::process(TDCHitWPPtrCollection *res,
                            TDCHitWPCollection *buf,
                            const TDCHitWPCollection& hits)
  {
    res->clear();
    res->reserve(hits.size());
    for(unsigned i=0; i<hits.size(); ++i) {
      res->push_back(TDCHitWPPtr(hits, i));
    }
  }

  //================================================================
  MOFIA_XTalkDiscarder::MOFIA_XTalkDiscarder(const std::string& topdir,
                                             WirePlane::DetType det,
                                             HistogramFactory& hf,
                                             const DetectorGeo& geom,
                                             const ConfigFile& conf)
  {
    const std::string hdir = topdir + "/" + WirePlane::detName(det);
    hxt_ = hf.DefineTH1D(hdir, "XTFlag", "XTFlag", 2, -0.5, 1.5);
  }

  //----------------------------------------------------------------
  void MOFIA_XTalkDiscarder::process(TDCHitWPPtrCollection *res,
                                     TDCHitWPCollection *buf,
                                     const TDCHitWPCollection& hits)
  {
    res->clear();
    res->reserve(hits.size());
    for(unsigned i=0; i<hits.size(); ++i) {
      hxt_->Fill(hits[i].xtalk());
      if(!hits[i].xtalk()) {
        res->push_back(TDCHitWPPtr(hits, i));
      }
    }
  }

  //================================================================
  NarrowHitDiscarder::NarrowHitDiscarder(const std::string& topdir,
                                         WirePlane::DetType det,
                                         HistogramFactory& hf,
                                         const DetectorGeo& geom,
                                         const ConfigFile& conf)
    : cutMinTDCWidth_(conf.read<float>("MuCapture/HitPreproc/"+WirePlane::detName(det)+"/NarrowHitDiscarder/cutMinTDCWidth"))
  {}

  //================================================================
  void NarrowHitDiscarder::process(TDCHitWPPtrCollection *res,
                                   TDCHitWPCollection * /*buf*/,
                                   const TDCHitWPCollection& hits)
  {
    res->clear();
    for(unsigned i=0; i<hits.size(); ++i) {
      if(hits[i].width() > cutMinTDCWidth_) {
        res->push_back(TDCHitWPPtr(hits, i));
      }
    }
  }

  //================================================================
  SameWireHitDiscarder::SameWireHitDiscarder(const std::string& topdir,
                                             WirePlane::DetType det,
                                             HistogramFactory& hf,
                                             const DetectorGeo& geom,
                                             const ConfigFile& conf)
    : cutSameWireDt_(conf.read<float>("MuCapture/HitPreproc/"+WirePlane::detName(det)+"/SameWireHitDiscarder/cutSameWireDt"))
  {
    const std::string hdir = topdir + "/" + WirePlane::detName(det);
    const unsigned maxPlaneNumber = geom.planes(det).size();
    const unsigned maxCellNumber = geom.maxCellNumber(det);

    const double dtHistoUpperLimit = 2000;
    const int dtHistoNbins = 200;

    for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
      std::ostringstream os;
      os<<"sameCellDtAll"<<std::setw(2)<<std::setfill('0')<<plane;
      hSameCellDtAll_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), dtHistoNbins, 0., dtHistoUpperLimit));
    }
    for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
      std::ostringstream os;
      os<<"sameCellDtDropped"<<std::setw(2)<<std::setfill('0')<<plane;
      hSameCellDtDropped_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), dtHistoNbins, 0., dtHistoUpperLimit));
    }

    for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
      std::ostringstream os;
      os<<"sameCellWidthDropped"<<std::setw(2)<<std::setfill('0')<<plane;
      hSameCellWidthDropped_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 200, 0., 1000.));
    }

    for(int plane = 1; plane <= maxPlaneNumber; ++plane) {
      std::ostringstream os;
      os<<"sameCellWidthKept"<<std::setw(2)<<std::setfill('0')<<plane;
      hSameCellWidthKept_.push_back(hf.DefineTH1D(hdir, os.str(), os.str(), 200, 0., 1000.));
    }

    hSameCellOccupancyDropped_.init(hdir, "occDropped", maxPlaneNumber, maxCellNumber, hf, conf);
    hSameCellOccupancyKept_.init(hdir, "occKept", maxPlaneNumber, maxCellNumber, hf, conf);
  }

  //================================================================
  void SameWireHitDiscarder::process(TDCHitWPPtrCollection *res,
                                     TDCHitWPCollection *buf,
                                     const TDCHitWPCollection& inputs)
  {
    PassThrough pt;
    TDCHitWPPtrCollection hits;
    pt.process(&hits, buf, inputs);

    std::sort(hits.begin(), hits.end(), TDCHitWPCmpTime());
    std::stable_sort(hits.begin(), hits.end(), TDCHitWPCmpGeom());
    // here we have hits sorted by cell then time

    TDCHitWPPtrCollection out;
    out.reserve(inputs.size());

    if(!hits.empty()) {
      out.push_back(hits[0]);
      for(int ihit = 1; ihit < hits.size(); ++ihit) {
        if(hits[ihit]->cid() == hits[ihit-1]->cid()) {
          const int plane = hits[ihit]->plane();
          const float dt = hits[ihit]->time() - hits[ihit-1]->time();
          hSameCellDtAll_[plane-1]->Fill(dt);

          if(dt < cutSameWireDt_) {
            hSameCellDtDropped_[plane-1]->Fill(dt);
            hSameCellWidthDropped_[plane-1]->Fill(hits[ihit]->width());
            hSameCellOccupancyDropped_.fill(*hits[ihit]);
          }
          else {
            out.push_back(hits[ihit]);
            hSameCellWidthKept_[plane-1]->Fill(hits[ihit]->width());
            hSameCellOccupancyKept_.fill(*hits[ihit]);
          }
        }
        else {
          out.push_back(hits[ihit]);
        }
      }
    }

    res->swap(out);
  }

  //================================================================

} // namespace TDCHitPreprocessing
