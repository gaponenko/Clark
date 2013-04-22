// Andrei Gaponenko, 2013

#include "MuCapContainmentCheck.h"

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
void MuCapContainmentCheck::init(HistogramFactory &hf,
                                 const std::string& hdir,
                                 const DetectorGeo& geom,
                                 const ConfigFile &conf)
{
  geom_ = &geom;

  cu_.init(hf, hdir+"/U", geom, conf);
  cv_.init(hf, hdir+"/V", geom, conf);

  huv_ = hf.DefineTH2D(hdir, "uv", "extrapolated V vs U",
                       250, 0., 25., 250, 0., 25.);

  huv_->SetOption("colz");

  hr_ = hf.DefineTH1D(hdir, "rmax", "R max", 300, 0., 30.);
}

//================================================================
double MuCapContainmentCheck::rmax(int lastPlane,
                                   const ClustersByPlane& globalClusters)
{
  if((lastPlane < 2) || (lastPlane > geom_->numGlobal() - 1)) {
    throw std::runtime_error("MuCapContainmentCheck::rmax(): lastPlane out of range");
  }

  const double toZ = geom_->global(lastPlane+1).center().z();

  const WirePlane& wpLast = geom_->global(lastPlane);
  const WirePlane::Measurement m2 =
    (wpLast.direction() == WirePlane::U ? cu_ : cv_ )
    .limit(lastPlane, toZ, globalClusters);

  // lastPlane-1 may be of the same direction as lastPlane
  // take care of this
  int prevPlane = lastPlane - 1;
  while(0 < prevPlane) {
    if(geom_->global(prevPlane).direction() != wpLast.direction()) {
      break;
    }

    --prevPlane;
  }
  if(0 == prevPlane) {
    throw std::runtime_error("MuCapContainmentCheck::rmax(): no suitable prevPlane");
  }

  const WirePlane& wpPrev = geom_->global(prevPlane);
  const WirePlane::Measurement m1 =
    (wpPrev.direction() == WirePlane::U ? cu_ : cv_ )
    .limit(prevPlane, toZ, globalClusters);

  const ROOT::Math::XYPoint res = WirePlane::uv(m1, m2);
  huv_->Fill(res.x(), res.y());
  hr_->Fill(res.R());

  return res.R();
}

//================================================================
