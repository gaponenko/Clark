// Andrei Gaponenko, 2014

#include "HistTDCBCSWidth.h"
#include <limits>
#include <cmath>
#include <algorithm>

#include "HitBasedObservables.h"
#include "EventClass.h"

#include "TH1.h"
#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

#include "TimeWindow.h"

extern TimeWindowingResults *gwres;

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

//================================================================
void HistTDCBCSWidth::init(HistogramFactory& hf,
                           const std::string& hdir,
                           const DetectorGeo& geom,
                           const ConfigFile& conf)
{
  geom_ = &geom;

  clusterSizePC_ = hf.DefineTH2D(hdir, "clusterSizePC", "PC plane vs cluster size", 10, -0.5, 9.5, 12, 0.5, 12.5);
  clusterSizePC_->SetOption("colz");
  clusterSizePC_->GetXaxis()->SetTitle("PC cluster size");
  clusterSizePC_->GetYaxis()->SetTitle("PC plane number");

  clusterSizeDC_ = hf.DefineTH2D(hdir, "clusterSizeDC", "DC plane vs cluster size", 10, -0.5, 9.5, 44, 0.5, 44.5);
  clusterSizeDC_->SetOption("colz");
  clusterSizeDC_->GetXaxis()->SetTitle("DC cluster size");
  clusterSizeDC_->GetYaxis()->SetTitle("DC plane number");

  clusterSizeBinsPC_.resize(maxClusterSize);
  for(int i=1; i<=maxClusterSize; ++i) {
    std::ostringstream os;
    os<<hdir<<"/pcClusterSize"<<i;
    clusterSizeBinsPC_[i-1].init(os.str(), "pcwidth", 12, hf, conf);
  }

  clusterSizeBinsDC_.resize(maxClusterSize);
  for(int i=1; i<=maxClusterSize; ++i) {
    std::ostringstream os;
    os<<hdir<<"/dcClusterSize"<<i;
    clusterSizeBinsDC_[i-1].init(os.str(), "dcwidth", 44, hf, conf);
  }

  perWireHitDistsPC_.init(hdir+"/pcPerWireDists", "pcwidth", 12, hf, conf);
  perWireHitDistsDC_.init(hdir+"/dcPerWireDists", "dcwidth", 44, hf, conf);
}

//================================================================
void HistTDCBCSWidth::fill(const EventClass& evt, const ClustersByPlane& protonGlobalClusters) {
  for(int iplane=1; iplane < protonGlobalClusters.size(); ++iplane) {
    if(protonGlobalClusters[iplane].size() == 1) {
      fill(evt, geom_->global(iplane).planeType(), protonGlobalClusters[iplane][0]);
    }
  }
}

//================================================================
void HistTDCBCSWidth::fill(const EventClass& evt, WirePlane::DetType pt, const WireCluster& cluster) {
  int sz = std::min(cluster.numCells(), maxClusterSize);
  switch(pt) {
  case WirePlane::PC:
    {
      clusterSizePC_->Fill(cluster.numCells(), cluster.plane());
      HistTDCWidth& hbin = clusterSizeBinsPC_[sz - 1];
      hbin.fill(cluster.hits());

      // Sort cluster hits into per-wire collections and fill the per-wire distributions
      TDCHitWPPtrCollection clusterhits(cluster.hits());
      std::sort(clusterhits.begin(), clusterhits.end(), TDCHitWPCmpGeom());
      for(int wirestart = 0; wirestart < clusterhits.size(); ) {
        TDCHitWPPtrCollection wirehits;
        wirehits.push_back(clusterhits[wirestart]);
        for(int i=wirestart+1; (i < clusterhits.size())&&(clusterhits[i]->cell() == wirehits.front()->cell()) ; ++i) {
          wirehits.push_back(clusterhits[i]);
        }

        perWireHitDistsPC_.fill(wirehits);

        wirestart += wirehits.size();
      }

      break;
    }
  case WirePlane::DC:
    {
      clusterSizeDC_->Fill(cluster.numCells(), cluster.plane());
      HistTDCWidth& hbin = clusterSizeBinsDC_[sz - 1];
      hbin.fill(cluster.hits());

      // Sort cluster hits into per-wire collections and fill the per-wire distributions
      TDCHitWPPtrCollection clusterhits(cluster.hits());
      std::sort(clusterhits.begin(), clusterhits.end(), TDCHitWPCmpGeom());
      for(int wirestart = 0; wirestart < clusterhits.size(); ) {
        TDCHitWPPtrCollection wirehits;
        wirehits.push_back(clusterhits[wirestart]);
        for(int i=wirestart+1; (i < clusterhits.size())&&(clusterhits[i]->cell() == wirehits.front()->cell()) ; ++i) {
          wirehits.push_back(clusterhits[i]);
        }

        perWireHitDistsDC_.fill(wirehits);

        wirestart += wirehits.size();
      }


      break;
    }
  }
}

//================================================================
