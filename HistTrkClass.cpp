// Andrei Gaponenko, 2015

#include "HistTrkClass.h"
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include "HitBasedObservables.h"
#include "EventClass.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

//================================================================
void HistTrkClass::init(HistogramFactory& hf,
                        const std::string& hdir,
                        const DetectorGeo& geom,
                        const ConfigFile& conf)
{
  //----------------------------------------------------------------
  hptot_all_  = hf.DefineTH1D(hdir, "ptot_all", "ptot, all", 500, 0., 500.);

  hptot_zcrc_ = hf.DefineTH1D(hdir, "ptot_zcrc", "ptot, zcrc", 500, 0., 500.);
  hptot_zcrn_ = hf.DefineTH1D(hdir, "ptot_zcrn", "ptot, zcrn", 500, 0., 500.);
  hptot_znrc_ = hf.DefineTH1D(hdir, "ptot_znrc", "ptot, znrc", 500, 0., 500.);
  hptot_znrn_ = hf.DefineTH1D(hdir, "ptot_znrn", "ptot, znrn", 500, 0., 500.);

  //----------------------------------------------------------------
  hpcos_all_  = hf.DefineTH2D(hdir, "pcos_all", "pcos, all", 500, 0., 500., 100, -1., +1.);

  hpcos_zcrc_ = hf.DefineTH2D(hdir, "pcos_zcrc", "pcos, zcrc", 500, 0., 500., 100, -1., +1.);
  hpcos_zcrn_ = hf.DefineTH2D(hdir, "pcos_zcrn", "pcos, zcrn", 500, 0., 500., 100, -1., +1.);
  hpcos_znrc_ = hf.DefineTH2D(hdir, "pcos_znrc", "pcos, znrc", 500, 0., 500., 100, -1., +1.);
  hpcos_znrn_ = hf.DefineTH2D(hdir, "pcos_znrn", "pcos, znrn", 500, 0., 500., 100, -1., +1.);
}

//================================================================
void HistTrkClass::fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters) {
  if(itrack != -1) {
    const double ptot = evt.ptot[itrack];
    const double costh = evt.costh[itrack];

    const int cutMaxPlane = 51;
    const double cutMaxRout = 15.;

    const int extrange = MuCapUtilities::findExtendedLastPlane(evt, itrack, protonGlobalClusters);
    const bool zc = (extrange <= cutMaxPlane);

    const double rout =
      sqrt(std::pow(evt.hefit_ucenter[itrack], 2) + std::pow(evt.hefit_vcenter[itrack], 2))
      + evt.radius[itrack];

    const bool rc = (rout <= cutMaxRout);

    //----------------
    hptot_all_->Fill(ptot);
    (zc ? (rc ? hptot_zcrc_ : hptot_zcrn_ ) : (rc ? hptot_znrc_ : hptot_znrn_) )->Fill(ptot);

    hpcos_all_->Fill(ptot, costh);
    (zc ? (rc ? hpcos_zcrc_ : hpcos_zcrn_ ) : (rc ? hpcos_znrc_ : hpcos_znrn_) )->Fill(ptot, costh);

  }
}

//================================================================
