// Andrei Gaponenko, 2014

#ifndef HistMuCapHitbasedChannel_h
#define HistMuCapHitbasedChannel_h

#include <string>

#include "HistMuCapTruth.h"
#include "HistMuCapFinal.h"
#include "HistHotSpot.h"
#include "TimeWindow.h"
#include "HistTDCBCSWidth.h"
#include "HistXTPlane.h"
#include "HitBasedObservables.h"
#include "MuCapPC78Cut.h"

#include "TAxis.h"

class TH1;
class TH2;
class TH3;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistMuCapHitbasedChannel {
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

  bool accepted(const EventClass& evt, const ClustersByPlane& globalPlaneClusters, int iDIOVetoTrack, bool referenceSampleAccepted);

  HistMuCapHitbasedChannel() : geom_(0), doMCTruth_(false), tdcWidthFilterCutPC_() {}

private :
  const DetectorGeo *geom_;
  bool doMCTruth_;
  double tdcWidthFilterCutPC_;
  int maxClusterWiresFilterCutPC_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  MuCapPC78Cut cutpc78_;

  TH2* lastconPlaneVsCWires_; // last plane in a contiguous range vs sum(largest cluster size)

  TH3* migration_;
  TH3* migration_mcproton_;
  TH3* migration_mcdeuteron_;

  TH3* contamination_;
  TH3* contamination_mcproton_;
  TH3* contamination_mcdeuteron_;

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

  TH2 *hFilterEffectPC7_;
  TH2 *hFilterEffectPC8_;

  //----------------------------------------------------------------
  HistTDCBCSWidth htdcwidthInput_;
  HistTDCBCSWidth htdcwidthDoubleFiltered_;

  // Understanding the "hot spot" feature in data
  HistHotSpot hshot_;
  HistHotSpot hscold_;

  TH1 *hshotrun_;

  HistXTPlane hxtplane100_;
  HistXTPlane hxtplane300_;
  //----------------------------------------------------------------

  CutNumber analyzeEvent(const EventClass& evt, const ClustersByPlane& globalPlaneClusters, int iDIOVetoTrack, bool referenceSampleAccepted);

  void filterDnPCNoise(ClustersByPlane *out, const ClustersByPlane& in);
  void filterClusterSize(WireClusterCollection *out, const WireClusterCollection& in);
  void fillFilterEffectHist(TH2* h, const WireClusterCollection& orig, const WireClusterCollection& filtered);
};

#endif/*HistMuCapHitbasedChannel_h*/
