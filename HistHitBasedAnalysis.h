// Andrei Gaponenko, 2014

#ifndef HistHitBasedAnalysis_h
#define HistHitBasedAnalysis_h

#include <string>

#include "HistMuCapTruth.h"
#include "HistMuCapFinal.h"
#include "HistHotSpot.h"
#include "TimeWindow.h"
#include "HistTDCBCSWidth.h"
#include "HistXTPlane.h"
#include "HitBasedObservables.h"

#include "TAxis.h"

class TH1;
class TH2;
class TH3;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistHitBasedAnalysis {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_DIOVETO, "DIO veto");
    ax->SetBinLabel(1+CUT_ZVETO, "Z veto");
    ax->SetBinLabel(1+CUT_NOPC7, "No PC7");
    ax->SetBinLabel(1+CUT_PCWIDTH, "PC TDC width");
    ax->SetBinLabel(1+CUTS_ACCEPTED, "Accepted");
  }

public:
  enum CutNumber {
    CUT_DIOVETO,
    CUT_ZVETO,
    CUT_NOPC7,
    CUT_PCWIDTH,
    CUTS_ACCEPTED,
    CUTS_END
  };

  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  bool accepted(const EventClass& evt, const ClustersByPlane& globalPlaneClusters, int iDIOVetoTrack);

  HistHitBasedAnalysis() : geom_(0), doMCTruth_(false) {}

private :
  const DetectorGeo *geom_;
  bool doMCTruth_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH2* lastconPlaneVsCWires_; // last plane in a contiguous range vs sum(largest cluster size)

  TH3* migration_;

  //----------------------------------------------------------------
  // Extra histograms to plot efficiencies etc.  Not essential for the
  // unfolding.

  // same content as the reco hists, but which one is filled depends on MC PID
  TH2* lastconPlaneVsCWires_mcproton_;
  TH2* lastconPlaneVsCWires_mcdeuteron_;
  TH2* lastconPlaneVsCWires_mcdio_;

  HistMuCapFinal noncontiguous_;
  HistHitBasedAmbiguities hambig_;

  HistMuCapTruth hTruth_in_;
  HistMuCapTruth hTruth_accepted_;

  TH1* hOuterVetoNumHitPlanes_;
  TH1* hNumPC7Clusters_;

  //----------------------------------------------------------------
  // Understanding the "hot spot" feature in data

  HistHotSpot hshot_;
  HistHotSpot hscold_;

  HistTDCBCSWidth htdcwidthMaxWires_;
  HistTDCBCSWidth htdcwidthMaxTDCWidth_;

  HistXTPlane hxtplane100_;
  HistXTPlane hxtplane300_;
  //----------------------------------------------------------------

  CutNumber analyzeEvent(const EventClass& evt, const ClustersByPlane& globalPlaneClusters, int iDIOVetoTrack);

};

#endif/*HistHitBasedAnalysis_h*/
