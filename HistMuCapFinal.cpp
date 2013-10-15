// Andrei Gaponenko, 2013

#include "HistMuCapFinal.h"

#include "TH1.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapFinal::init(HistogramFactory& hf,
                          const std::string& hdir,
                          const DetectorGeo& geom,
                          const ConfigFile& conf)
{
  hLastPlane_ = hf.DefineTH1D(hdir, "lastPlane", "Last plane final", 56, 0.5, 56.5);

  hNumPlanesVsWires_ = hf.DefineTH2D(hdir, "numPlanesVsWires", "num planes vs wires",
                                     200, -0.5, 199.5,
                                     1+geom.numGlobal()/2, -0.5, geom.numGlobal()/2+0.5);

  hNumPlanesVsWires_->SetOption("colz");
}

//================================================================
void HistMuCapFinal::fill(const ClustersByPlane& globalPlaneClusters,
                          const EventClass& evt)
{
  const PlaneRange dn = findDownstreamPlaneRange(globalPlaneClusters);

  hLastPlane_->Fill(dn.max());

  int numWires(0);
  int numPlanes(0);
  for(int plane=29; plane<globalPlaneClusters.size(); ++plane) {
    if(!globalPlaneClusters[plane].empty()) {
      ++numPlanes;
    }
    numWires += MuCapUtilities::numWires(globalPlaneClusters[plane]);
  }

  hNumPlanesVsWires_->Fill(numWires, numPlanes);

  // Investigate the peak that is not in MC
  if(false/*this block is normally disabled*/
     && (numPlanes == 2)&&(numWires > 90)) {
    ///static std::ofstream dataPeak("dataPeak.txt");
    //dataPeak<<evt.nrun<<" "<<evt.nevt<<std::endl;
    std::cout<<"Data peak: run = "<<evt.nrun<<", event = "<<evt.nevt<<std::endl;

    for(int plane=29; plane<globalPlaneClusters.size(); ++plane) {
      if(!globalPlaneClusters[plane].empty()) {
        std::cout<<"   Dump: plane = "<<plane<<", num clusters = "
                 <<globalPlaneClusters[plane].size()<<std::endl;
        for(int ic = 0; ic < globalPlaneClusters[plane].size(); ++ic) {
          std::cout<<"      Dump for cluster "<<ic<<": hits = "<<std::endl;;
          const TDCHitWPPtrCollection& hits = globalPlaneClusters[plane][ic].hits();
          for(int ihit = 0; ihit < hits.size(); ++ihit) {
            std::cout<<"         hit wire = "<<hits[ihit]->cell()
                     <<", time = "<<hits[ihit]->time()
                     <<", width = "<<hits[ihit]->width()
                     <<std::endl;
          }
        }
      }
    }
  }
}

//================================================================
