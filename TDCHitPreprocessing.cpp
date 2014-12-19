// Andrei Gaponenko, 2013

#include "TDCHitPreprocessing.h"

#include <algorithm>
#include <iterator>
#include <sstream>

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"

namespace TDCHitPreprocessing {

  //================================================================
  void PassThrough::process(TDCHitWPPtrCollection *res,
                            const TDCHitWPPtrCollection& hits)
  {
    res->clear();
    res->reserve(hits.size());
    std::copy(hits.begin(), hits.end(), std::back_inserter(*res));
  }

  //================================================================
  MOFIA_XTalkDiscarder::MOFIA_XTalkDiscarder(const std::string& topdir,
                                             WirePlane::DetType det,
                                             HistogramFactory& hf,
                                             const DetectorGeo& geom,
                                             const ConfigFile& conf)
  {
    const std::string hdir = topdir + "/" + WirePlane::detName(det) + "_XTalkDiscarder";
    hxt_ = hf.DefineTH1D(hdir, "XTFlag", "XTFlag", 2, -0.5, 1.5);
  }

  //----------------------------------------------------------------
  void MOFIA_XTalkDiscarder::process(TDCHitWPPtrCollection *res,
                                     const TDCHitWPPtrCollection& hits)
  {
    res->clear();
    res->reserve(hits.size());
    for(unsigned i=0; i<hits.size(); ++i) {
      hxt_->Fill(hits[i]->xtalk());
      if(!hits[i]->xtalk()) {
        res->push_back(hits[i]);
      }
    }
  }

  //================================================================
  NarrowHitDiscarder::NarrowHitDiscarder(const std::string& topdir,
                                         WirePlane::DetType det,
                                         HistogramFactory& hf,
                                         const DetectorGeo& geom,
                                         const ConfigFile& conf)
  {
    const std::string hdir = topdir + "/" + WirePlane::detName(det) + "_NarrowHitDiscarder";

    cutMinTDCWidth_.resize(1 + ((det == WirePlane::PC) ? geom.numPCs() : geom.numDCs()));
    hwidth_.resize(cutMinTDCWidth_.size());

    const std::string globalCutKey = "MuCapture/HitPreproc/"+WirePlane::detName(det)+"/NarrowHitDiscarder/cutMinTDCWidth";
    const std::string arrayCutKey = globalCutKey + "Array";

    if(conf.keyExists(arrayCutKey) && conf.keyExists(globalCutKey)) {
      throw std::runtime_error("Config error: both "+globalCutKey+" and "
                               +arrayCutKey+" are set.  Only one of them must be specified.");
    }

    if(conf.keyExists(globalCutKey)) {
      const double val = conf.read<float>(globalCutKey);
      for(int i=0; i<cutMinTDCWidth_.size(); ++i) {
        cutMinTDCWidth_[i] = val;
      }
      std::cout<<"NarrowHitDiscarder: using global "
               <<(det==WirePlane::PC ? "PC":"DC")
               <<" TDC width cut = "<<val
               <<std::endl;
    }
    else {
      const std::string inputs = conf.read<std::string>(arrayCutKey);
      std::istringstream iss(inputs);
      for(int i=1; i<cutMinTDCWidth_.size(); ++i) {
        if(!(iss>>cutMinTDCWidth_[i])) {
          std::ostringstream os;
          os<<"Error reading cut for plane "<<i<<" from the input string: \""<<inputs<<"\" (config key "<<arrayCutKey<<")";
          throw std::runtime_error(os.str());
        }
      }

      double dummy=1234.;
      if(iss>>dummy) {
        std::cerr<<"dummy = "<<dummy<<" cutMinTDCWidth_.size() = "<<cutMinTDCWidth_.size()<<std::endl;
        std::copy(cutMinTDCWidth_.begin(), cutMinTDCWidth_.end(), std::ostream_iterator<double>(std::cerr, " "));
        throw std::runtime_error("Error: too many cut values in the input string: \""
                                 +inputs+"\" (config key "+arrayCutKey+")");
      }

      std::cout<<"NarrowHitDiscarder: using per-plane "
               <<(det==WirePlane::PC ? "PC":"DC")
               <<" TDC width cuts = ";
      std::copy(/*skip the 0-th entry*/ ++cutMinTDCWidth_.begin(),
                cutMinTDCWidth_.end(), std::ostream_iterator<double>(std::cout, " "));
      std::cout<<std::endl;
    }

    for(int i=1; i<hwidth_.size(); ++i) {
      std::ostringstream osn;
      osn<<"width"<<std::setw(2)<<std::setfill('0')<<i;

      std::ostringstream ost;
      ost<<"TDC width for "<<(det==WirePlane::PC ? "PC":"DC")<<" plane "<<i;
      hwidth_[i] = hf.DefineTH1D(hdir, osn.str(), "TDC width", 800, -0.5, 399.5);
    }
  }

  //================================================================
  void NarrowHitDiscarder::process(TDCHitWPPtrCollection *res,
                                   const TDCHitWPPtrCollection& hits)
  {
    res->clear();
    for(unsigned i=0; i<hits.size(); ++i) {
      hwidth_[hits[i]->plane()]->Fill(hits[i]->width());
      if(hits[i]->width() > cutMinTDCWidth_[hits[i]->plane()]) {
        res->push_back(hits[i]);
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
    const std::string hdir = topdir + "/" + WirePlane::detName(det)+"_SameWireHitDiscarder";
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
                                     const TDCHitWPPtrCollection& inputs)
  {
    PassThrough pt;
    TDCHitWPPtrCollection hits;
    pt.process(&hits, inputs);

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
  Hits::Hits(const TDCHitWPCollection& in) {
    phits_.reserve(in.size());
    for(unsigned i=0; i<in.size(); ++i) {
      phits_.push_back(TDCHitWPPtr(in, i));
    }
  }

  //================================================================

} // namespace TDCHitPreprocessing
