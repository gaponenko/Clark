// Andrei Gaponenko, 2013

#include "HistMCElectrons.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "TH2.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMCElectrons::init(HistogramFactory &hf,
                           const std::string& hdir,
                           const ConfigFile &conf)
{
  hNumAllElectrons_ = hf.DefineTH1D(hdir, "numAllMCElectrons", "num MC electrons, all", 100, -0.5, 99.5);
  hNumAllElectrons_->StatOverflows();

  hNumElectrons_ = hf.DefineTH1D(hdir, "numMCElectrons", "num MC electrons, non-upstream production", 100, -0.5, 99.5);
  hNumElectrons_->StatOverflows();

  hEkAll_ = hf.DefineTH1D(hdir, "ekall", "MC electron kinetic energy, all", 500, 0., 50.);
  hEk_ = hf.DefineTH1D(hdir, "ek", "MC electron kinetic energy, non-upstream production", 500, 0., 50.);

  hZStartAll_ = hf.DefineTH1D(hdir, "zstartall", "Electron track start Z, all", 1200, -60., 60.);
  hZStart_ = hf.DefineTH1D(hdir, "zstart", "Electron track start Z, non-upstream production", 1200, -60., 60.);
}

//================================================================
void HistMCElectrons::fill(const EventClass& evt) {
  if(evt.mctype==EventClass::G3) {
    int numAllElectrons = 0;
    int numElectrons = 0;

    for(int imctrk=0; imctrk<evt.nmctr; ++imctrk) {
      if(evt.mctrack_pid[imctrk] == MuCapUtilities::PID_G3_EMINUS) {
        ++numAllElectrons;

        const int imcvtx = evt.getFirstMCVertexIndexForTrack(imctrk);
        if(evt.mcvertex_istop[imcvtx] != 0) {
          throw std::runtime_error("HistMCElectrons: mcvertex_istop != 0 for the initial vertex");
        }

        hEkAll_->Fill(MuCapUtilities::kineticEnergy(evt.mctrack_pid[imctrk], evt.mcvertex_ptot[imcvtx], evt));
        hZStartAll_->Fill(evt.mcvertex_vz[imcvtx]);

        if(evt.mcvertex_vz[imcvtx] > -0.018 /*Include the target*/) {
          ++numElectrons;
          hEk_->Fill(MuCapUtilities::kineticEnergy(evt.mctrack_pid[imctrk], evt.mcvertex_ptot[imcvtx], evt));
          hZStart_->Fill(evt.mcvertex_vz[imcvtx]);
        }
      }
    }

    hNumAllElectrons_->Fill(numAllElectrons);
    hNumElectrons_->Fill(numElectrons);
  }
}

//================================================================
