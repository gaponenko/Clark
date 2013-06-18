// Andrei Gaponenko, 2013

#include "HistProtonWindow.h"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "TH1.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
void HistProtonWindow::init(HistogramFactory &hf,
                            const std::string& hdir,
                            const DetectorGeo& geom,
                            const ConfigFile &conf)
{
  hClustersVsPlane_ = hf.DefineTH2D(hdir, "numClustersVsPlane", "num clusters vs plane",
                                    1+geom.numGlobal(), -0.5, geom.numGlobal()+0.5,
                                    10, -0.5, 9.5);

  hClustersVsPlane_->SetOption("colz");

  hMaxClustersVsLastPlane_ = hf.DefineTH2D(hdir, "maxClustersVsLastPlane", "max clusters vs last plane",
                                           1+geom.numGlobal(), -0.5, geom.numGlobal()+0.5,
                                           10, -0.5, 9.5);

  hMaxClustersVsLastPlane_->SetOption("colz");

  hClustersVsRemainingRange_ = hf.DefineTH2D(hdir, "numClustersVsRemainingRange",
                                             "num clusters vs remaining range",
                                             29, -0.5, 28.5,
                                             10, -0.5, 9.5);

  hClustersVsRemainingRange_->SetOption("colz");
}

//================================================================
void HistProtonWindow::fill(const ClustersByPlane& globalPlaneClusters,
                            const EventClass& evt)
{
  const PlaneRange gr = findPlaneRange(globalPlaneClusters);
  // All selection cuts were done upstream, we are not cutting
  // on range gaps or anything here

  int maxClusters(0);
  for(int plane=0; plane<globalPlaneClusters.size(); ++plane) {
    const int numClusters = globalPlaneClusters[plane].size();
    maxClusters=std::max(maxClusters, numClusters);
    hClustersVsPlane_->Fill(plane, numClusters);

    const int remainingRange = gr.max() - plane;
    if(remainingRange >= 0) {
      hClustersVsRemainingRange_->Fill(remainingRange, numClusters);
    }
  }

  hMaxClustersVsLastPlane_->Fill(gr.max(), maxClusters);
}

//================================================================
