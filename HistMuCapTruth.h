// These are the truth-only plots that can be made without any event
// selection.  That is, the do not use reconstructed quantities.
//
// Andrei Gaponenko, 2013

#ifndef HistMuCapTruth_h
#define HistMuCapTruth_h

#include <string>

#include "HistMCElectrons.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistMuCapTruth {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const ConfigFile &conf);

  void fill(const EventClass& evt);

  HistMuCapTruth()
    : hNumMCCaptureTracks_()
    , hCaptureTime_()
    , hptot_()
    , hptot_proton_()
    , hptot_deuteron_()
    , hptot_triton_()
    , hptot_alpha_()
    , hek_()
    , hbeta_()
    , hphi_()
    , hpcos_()
    , hVUend_()
    , hRZend_()
    , hRendVsPstart_()
    , hZendVsPstart_()
    , hZStart1_()
    , hZStart2_()
    , hZStop1_()
    , hZStop2_()
    , hZStop3_()
    , hZStop4_()
  {}

private :
  HistMCElectrons helectrons_;
  TH1 *hNumMCCaptureTracks_;
  TH1 *hCaptureTime_;
  TH1 *hptot_;
  TH1 *hptot_proton_;
  TH1 *hptot_deuteron_;
  TH1 *hptot_triton_;
  TH1 *hptot_alpha_;
  TH1 *hek_;
  TH1 *hbeta_;
  TH1 *hphi_;
  TH2 *hpcos_;
  TH2 *hVUend_;
  TH2 *hRZend_;
  TH2 *hRendVsPstart_;
  TH2 *hZendVsPstart_;
  TH1 *hZStart1_;
  TH1 *hZStart2_;
  TH1 *hZStop1_;
  TH1 *hZStop2_;
  TH1 *hZStop3_;
  TH1 *hZStop4_;
};

#endif/*HistMuCapTruth_h*/
