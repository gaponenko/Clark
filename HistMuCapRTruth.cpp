// Andrei Gaponenko, 2013

#include "HistMuCapRTruth.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"

//================================================================
void HistMuCapRTruth::init(HistogramFactory &hf,
                           const std::string& hdir,
                           const ConfigFile &conf)
{
  hLastPlaneVsMCPstart_ = hf.DefineTH2D(hdir, "lastPlaneVsMCPstart", "Last plane vs MC pstart",
                                        200, 0., 200., 56, 0.5, 56.5);

  hLastPlaneVsMCPstart_->SetOption("colz");

  hRmaxContained_ = hf.DefineTH1D(hdir, "rmaxContained", "R max for MC Rend<16 cm", 300, 0., 30.);
  hRmaxUncontained_ = hf.DefineTH1D(hdir, "rmaxUncontained", "R max for MC Rend>16 cm", 300, 0., 30.);
}

//================================================================
void HistMuCapRTruth::fill(const EventClass& evt, int lastPlane, double extrapolatedRmax) {
  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  const unsigned imcvtxEnd = evt.iCaptureMcVtxEnd;
  const unsigned imctrk = evt.iCaptureMcTrk;

  hLastPlaneVsMCPstart_->Fill(evt.mcvertex_ptot[imcvtxStart], lastPlane);

  const double mcRend = std::sqrt(std::pow(evt.mcvertex_vu[imcvtxEnd], 2) +
                                  std::pow(evt.mcvertex_vv[imcvtxEnd], 2));

  if(std::abs(evt.mcvertex_vz[imcvtxEnd]) < 60.) {
    (mcRend < 16. ? hRmaxContained_ : hRmaxUncontained_)->Fill(extrapolatedRmax);
  }
}

//================================================================
