// Andrei Gaponenko, 2013

#include "HistXTPlane.h"

#include <sstream>
#include <iomanip>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "TH1.h"
#include "TProfile2D.h"

#include "DetectorGeo.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
void HistXTPlane::init(HistogramFactory& hf,
                       const std::string& hdir,
                       const DetectorGeo& geom,
                       const ConfigFile& conf,
                       double cutBroadHitWidth)
{
  geom_ = &geom;
  cutNarrowHitWidth_[WirePlane::PC] = 40.;
  cutNarrowHitWidth_[WirePlane::DC] = 50.;
  cutBroadHitWidth_[WirePlane::PC] = cutBroadHitWidth;
  cutBroadHitWidth_[WirePlane::DC] = cutBroadHitWidth;

  // Initalize the array with NULLs
  // Then we'll put non-NULL values for histograms we want to look at.
  hdt_.resize(1+geom_->numGlobal()); // will use MOFIA plane numbering

  //  Look at planes going to the same crate/slot in fbc1_map.00045,
  //  fbc2_map.00045.  List slots for all downstream PCs, than see
  //  what else goes to the same slot.
  //
  // Crate 1:
  //
  // slot  plane
  // 3     PC7
  // 3     PC8
  //
  // 22    PC9 (1-80)
  // 22    DC42
  // 22    DC44
  //
  // 23    PC9 (81-160)
  // 23    PC11
  //
  // Crate 2:
  //
  // 22    PC10
  // 22    DC41
  // 22    DC43
  //
  // 23    PC10
  // 23    PC12

  // look at the above combinations as well as some other pairs
  bookHisto(hf, hdir, 29, 27); // pc7 <- pc5
  bookHisto(hf, hdir, 29, 28); // pc7 <- pc6
  bookHisto(hf, hdir, 29, 30); // pc7 <- pc8
  bookHisto(hf, hdir, 29, 31); // pc7 <- dc23
  bookHisto(hf, hdir, 29, 53); // pc7 <- pc9

  bookHisto(hf, hdir, 30, 28); // pc8 <- pc6
  bookHisto(hf, hdir, 30, 29); // pc8 <- pc7
  bookHisto(hf, hdir, 30, 31); // pc8 <- dc23
  bookHisto(hf, hdir, 30, 32); // pc8 <- dc24
  bookHisto(hf, hdir, 30, 53); // pc8 <- pc9

  bookHisto(hf, hdir, 53, 30); // pc9 <- pc8
  bookHisto(hf, hdir, 53, 50); // pc9 <- dc42
  bookHisto(hf, hdir, 53, 51); // pc9 <- dc43
  bookHisto(hf, hdir, 53, 52); // pc9 <- dc44
  bookHisto(hf, hdir, 53, 54); // pc9 <- pc10
  bookHisto(hf, hdir, 53, 55); // pc9 <- pc11

  bookHisto(hf, hdir, 54, 49); // pc10 <- dc41
  bookHisto(hf, hdir, 54, 51); // pc10 <- dc43
  bookHisto(hf, hdir, 54, 53); // pc10 <- pc9
  bookHisto(hf, hdir, 54, 55); // pc10 <- pc11
  bookHisto(hf, hdir, 54, 56); // pc10 <- pc12

  bookHisto(hf, hdir, 55, 53); // pc11 <- pc9
  bookHisto(hf, hdir, 55, 54); // pc11 <- pc10
  bookHisto(hf, hdir, 55, 56); // pc11 <- pc12

  bookHisto(hf, hdir, 56, 54); // pc12 <- pc10
  bookHisto(hf, hdir, 56, 55); // pc12 <- pc11
}

//================================================================
void HistXTPlane::bookHisto(HistogramFactory& hf, const std::string& hdir, int induced, int inducing) {
  std::ostringstream ost;
  ost<<"t(narrow in plane "<<induced<<") - t(wide in plane "<<inducing<<")";
  std::ostringstream osn;
  osn<<std::setw(2)<<std::setfill('0')<<"dtxt_ "<<induced<<"_"<<inducing;

  if(hdt_[induced].empty()) {
    hdt_[induced].resize(1+geom_->numGlobal());
  }

  hdt_[induced][inducing] = hf.DefineTH1D(hdir, osn.str(), ost.str(), 200, -50., 150.);
}

//================================================================
void HistXTPlane::fill(const EventClass& evt, const ClustersByPlane& globalClusters) {
  for(int induced=1; induced < hdt_.size(); ++induced) {
    if(!hdt_[induced].empty()) {

      TDCHitWPPtrCollection narrowHits;
      for(int ic=0; ic<globalClusters[induced].size(); ++ic) {
        for(int ih=0; ih<globalClusters[induced][ic].hits().size(); ++ih) {
          if(globalClusters[induced][ic].hits()[ih]->width()
             < cutNarrowHitWidth_[geom_->global(induced).planeType()]) {
            narrowHits.push_back(globalClusters[induced][ic].hits()[ih]);
          }
        }
      }

      for(int inducing=1; inducing < hdt_[induced].size(); ++inducing) {
        if(hdt_[induced][inducing]) {

          for(int ic=0; ic<globalClusters[inducing].size(); ++ic) {
            for(int ih=0; ih<globalClusters[inducing][ic].hits().size(); ++ih) {
              if(globalClusters[inducing][ic].hits()[ih]->width()
                 > cutBroadHitWidth_[geom_->global(inducing).planeType()]) {

                const double wideTime = globalClusters[inducing][ic].hits()[ih]->time();
                for(int i=0; i<narrowHits.size(); ++i) {
                  hdt_[induced][inducing]->Fill(narrowHits[i]->time() - wideTime);
                }
              }
            }
          }
        }
      }
    }
  }
}

//================================================================
