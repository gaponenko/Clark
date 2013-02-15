// Andrei Gaponenko, 2013

#include "MuCapture.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TH1.h"
#include "TH2.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)


bool MuCapture::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
  Log     = TmpLog;
  Log->info( "Register MuCapture module");
  //    -------- Name of the cut ---------     //
  Name            = "mu capture";

  //       --------- Special list of cut in this class for the global histograms ---------                //
  H.AddCut(Name);

  if (not Conf.read<bool>("MuCapture/Do")) {
    Log->info( "MuCapture code will not be run");
    return false ;
  }

  //       --------- Histograms initialization ---------          //
  h_dt_pc_ = H.DefineTH1D( "MuCapture", "dt_pc",   "dt to next PC hit", 1000, 0., 5000.);
  h_dt_dc_ = H.DefineTH1D( "MuCapture", "dt_dc",   "dt to next DC hit", 1000, 0., 5000.);
  h_dt_any_ = H.DefineTH1D( "MuCapture", "dt_any",   "dt to next wire hit", 1000, 0., 5000.);
  h_dt_any2_ = H.DefineTH2D( "MuCapture", "dt_any2",   "dt to next wire hit vs dt to previous hit", 1000, 0., 5000., 1000, 0., 5000.);

  h_m12_t0_ = H.DefineTH1D( "MuCapture", "m12t0",   "m12 time closest to 0", 201, -50.5, 50.5);

  //       --------- Parameters initialization ---------          //
  doDefaultTWIST = Conf.read<bool>("MuCapture/doDefaultTWIST");
  return true;
}

bool MuCapture::Process(EventClass &E, HistogramFactory &H)
{

  AGDEBUG("pc_nhits = "<<E.pc_nhits);
  std::vector<double> pc_times(E.pc_time, E.pc_time+E.pc_nhits);
  std::sort(pc_times.begin(), pc_times.end());

  std::vector<double> dc_times(E.dc_time, E.dc_time+E.dc_nhits);
  std::sort(dc_times.begin(), dc_times.end());

  AGDEBUG("here");

  if(pc_times.size() > 1) {
    for(int ihit=0; ihit < pc_times.size()-1; ++ihit) {
      AGDEBUG("here: h_dt_pc_ = "<<h_dt_pc_<<", ihit = "<<ihit);
      h_dt_pc_->Fill(pc_times[ihit+1] - pc_times[ihit]);
    }
  }

  AGDEBUG("here");

  if(dc_times.size() > 1) {
    for(int ihit=0; ihit < dc_times.size()-1; ++ihit) {
      AGDEBUG("here: h_dt_dc_ = "<<h_dt_dc_<<", ihit = "<<ihit);
      h_dt_dc_->Fill(dc_times[ihit+1] - dc_times[ihit]);
    }
  }

  std::vector<double> all_times(pc_times);
  all_times.insert(all_times.end(), dc_times.begin(), dc_times.end());
  std::sort(all_times.begin(), all_times.end());
  if(all_times.size() > 1) {
    for(int ihit=0; ihit < all_times.size()-1; ++ihit) {
      AGDEBUG("here: h_dt_any_ = "<<h_dt_any_<<", ihit = "<<ihit);
      h_dt_any_->Fill(all_times[ihit+1] - all_times[ihit]);
    }
  }

  if(all_times.size() > 2) {
    for(int ihit=1; ihit < all_times.size()-1; ++ihit) {
      AGDEBUG("here: h_dt_any_ = "<<h_dt_any_<<", ihit = "<<ihit);
      h_dt_any2_->Fill(all_times[ihit] - all_times[ihit-1], all_times[ihit+1] - all_times[ihit]);
    }
  }

  //----------------
  if(E.sc_nhits > 0) {
    double t0 = E.sc_time[0];
    for(int ihit=1; ihit < E.sc_nhits; ++ihit) {
      if(std::abs(E.sc_time[ihit]) < std::abs(t0)) {
        t0 = E.sc_time[ihit];
      }
    }
    h_m12_t0_->Fill(t0);
  }

  //----------------

  return doDefaultTWIST;
}
