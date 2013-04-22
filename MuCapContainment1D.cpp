// Andrei Gaponenko, 2013

#include "MuCapContainment1D.h"

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
void MuCapContainment1D::init(HistogramFactory &hf,
                              const std::string& hdir,
                              const DetectorGeo& geom,
                              const ConfigFile &conf)
{
  geom_ = &geom;
  minPlanedz_ = conf.read<double>("MuCapture/ProtonWindow/Containment1D/minPlanedz");

  hNumClusters_ = hf.DefineTH2D(hdir, "numClusters", "last vs prev clusters",
                                7, -0.5, 6.5, 7, -0.5, 6.5);

  hNumClusters_->SetOption("colz");

  hSlopeAll_ = hf.DefineTH1D(hdir, "slopeAll", "slope, all windows",
                             241, -6.025, 6.025);

  hSlopeMulticluster_ = hf.DefineTH1D(hdir, "slopeMulti", "slope, multicluster windows",
                                      241, -6.025, 6.025);
}

//================================================================
WirePlane::Measurement MuCapContainment1D::limit(int lastPlane,
                                                 double extrapolateToZ,
                                                 const ClustersByPlane& globalClusters)
{
  const WireClusterCollection& clLast = globalClusters[lastPlane];
  if(clLast.empty()) {
    throw std::runtime_error("MuCapContainment1D::limit(): no clusters at input lastPlane");
  }

  const WirePlane& wpLast = geom_->global(lastPlane);
  const double z1 = wpLast.center().z();

  const int step = (z1 < extrapolateToZ) ? -1 : +1;

  // Find the measurement closest to lastPlane of the right type
  int iplane = lastPlane+step;
  while( (iplane > 0)&&(iplane <= geom_->numGlobal())) {
    const WirePlane& wp = geom_->global(iplane);

    if((wp.direction() == wpLast.direction())&&
       (!globalClusters[iplane].empty())&&
       (std::abs(z1 - geom_->global(iplane).center().z()) > minPlanedz_))
      {
        break;
      }

    iplane += step;
  }

  if( (iplane <= 0)||(iplane > geom_->numGlobal())) {
    throw std::runtime_error("MuCapContainment1D::limit(): not enough data to exrapolate");
  }

  const double z0 = geom_->global(iplane).center().z();
  const WireClusterCollection& clPrev = globalClusters[iplane];
  const WirePlane& wpPrev = geom_->global(iplane);
  hNumClusters_->Fill(clPrev.size(), clLast.size());

  const bool multicluster = (clPrev.size()>1)||(clLast.size()>1);
  double absymax(0);
  for(int iprev=0; iprev<clPrev.size(); ++iprev) {
    const double y0 = wpPrev.measurement(clPrev[iprev].centralCell()).coordinate;
    for(int ilast = 0; ilast < clLast.size(); ++ilast) {
      const double y1 = wpLast.measurement(clLast[ilast].centralCell()).coordinate;
      const double slope = (y1 - y0)/(z1 - z0);
      const double y = y1 + slope*(extrapolateToZ - z1);
      absymax = std::max(std::abs(y), absymax);
      hSlopeAll_->Fill(slope);
      if(multicluster) {
        hSlopeMulticluster_->Fill(slope);
      }
    }
  }

  return WirePlane::Measurement(wpLast.direction(), absymax);
}

//================================================================
