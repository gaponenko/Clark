#include "MuCapPACTQuadrant.h"

#include <algorithm>
#include <stdexcept>

#include "TH2.h"

#include "HistogramFactory.h"
#include "DetectorGeo.h"
#include "ConfigFile.h"
#include "WireCluster.h"
#include "EventClass.h"

//================================================================
MuCapPACTQuadrant::MuCapPACTQuadrant(HistogramFactory &hf, const DetectorGeo& geom, const ConfigFile& conf,
                                     const std::string& hdir, const std::string& suffix)
  : slopea_(conf.read<double>("MuCapture/PACT/slopea"+suffix))
  , intercepta_(conf.read<double>("MuCapture/PACT/intercepta"+suffix))
  , slopeb_(conf.read<double>("MuCapture/PACT/slopeb"+suffix))
  , interceptb_(conf.read<double>("MuCapture/PACT/interceptb"+suffix))
  , hpc6vs5widthAll_(hf.DefineTH2D(hdir, "pc6vs5widthAll"+suffix, "PC6 vs 5 TDC width, all"+suffix, 500, 0., 500.,  500, 0., 500.))
  , hpc6vs5widthQ1_(hf.DefineTH2D(hdir, "pc6vs5widthQ1"+suffix, "PC6 vs 5 TDC width, quadrant 1"+suffix, 500, 0., 500.,  500, 0., 500.))
  , dia_vs_dib_all_(hf.DefineTH2D(hdir, "dia_vs_dib_all", "Delta(ia) vs Delta(ib)", 1200, -600., 600.,  800, -300., 500.))
  , dia_vs_dib_q1_(hf.DefineTH2D(hdir, "dia_vs_dib_q1", "Delta(ia) vs Delta(ib)", 1200, -600., 600.,  800, -300., 500.))

  , doMCTruth_(conf.read<bool>("TruthBank/Do"))
  , targetCenterZ_(geom.zTargetCenter())
  , targetThickness_(geom.targetThickness())
  , pc6CenterZ_(geom.pc(6).center().z())
  , pc6wireRadius_(geom.wireRadiusPC())

  , mctruthTargetStops_()
  , mctruthWireStops_()
  , mctruthOtherStops_()
{
  hpc6vs5widthAll_->SetOption("colz");
  hpc6vs5widthQ1_->SetOption("colz");
  dia_vs_dib_all_->SetOption("colz");
  dia_vs_dib_q1_->SetOption("colz");


  if(doMCTruth_) {
    mctruthTargetStops_   = hf.DefineTH2D(hdir, "mcTargetStops", "PC6 vs 5 TDC width, MC target stops", 500, 0., 500.,  500, 0., 500.);
    mctruthTargetStops_->SetOption("colz");

    mctruthWireStops_   = hf.DefineTH2D(hdir, "mcWireStops", "PC6 vs 5 TDC width, MC PC6 wire stops", 500, 0., 500.,  500, 0., 500.);
    mctruthWireStops_->SetOption("colz");

    mctruthOtherStops_   = hf.DefineTH2D(hdir, "mcOtherStops", "PC6 vs 5 TDC width, MC other stops", 500, 0., 500.,  500, 0., 500.);
    mctruthOtherStops_->SetOption("colz");
  }
}

//================================================================
int MuCapPACTQuadrant::quadrant(const WireCluster& pc5cluster,
                                const WireCluster& pc6cluster,
                                const EventClass& evt)
{
  const double pc5width = pc5cluster.totalTDCWidth();
  const double pc6width = pc6cluster.totalTDCWidth();


  hpc6vs5widthAll_->Fill(pc5width, pc6width);

  const double linea = intercepta_ + slopea_ * pc5width;
  const double lineb = interceptb_ + slopeb_ * pc5width;

  int res = (pc6width < linea) ?
    ((pc6width < lineb)? 4 : 3):
    ((pc6width < lineb)? 1 : 2);

  if(res == 1) {
    hpc6vs5widthQ1_->Fill(pc5width, pc6width);
  }

  const double dia = pc6width - slopea_ * pc5width - intercepta_;
  const double dib = pc6width - slopeb_ * pc5width - interceptb_;
  dia_vs_dib_all_->Fill(dib, dia);
  if(res==1) {
    dia_vs_dib_q1_->Fill(dib, dia);
  }


  if(doMCTruth_) {
    switch(muStopKind(evt)) {
    case MuStopRegion::TARGET: mctruthTargetStops_->Fill(pc5width, pc6width); break;
    case MuStopRegion::PC6WIRE: mctruthWireStops_->Fill(pc5width, pc6width); break;
    case MuStopRegion::OTHER: mctruthOtherStops_->Fill(pc5width, pc6width); break;
    default:
      throw std::runtime_error("MuCapPACTQuadrant: internal errror interpreting muStopKind() return value");
    }
  }

  return res;
}

//================================================================
MuCapPACTQuadrant::MuStopRegion
MuCapPACTQuadrant::muStopKind(const EventClass& evt) const {
  MuStopRegion res = MuStopRegion::OTHER;

  if(evt.iMuStopMcVtxEnd != -1) {
    const double zstop = evt.mcvertex_vz[evt.iMuStopMcVtxEnd];
    if(std::abs(zstop - targetCenterZ_) <= targetThickness_/2.) {
      res = MuStopRegion::TARGET;
    }
    else if(std::abs(zstop - pc6CenterZ_) <= pc6wireRadius_) {
      res = MuStopRegion::PC6WIRE;
    }
  }

  return res;
}

//================================================================
