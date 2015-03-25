// Andrei Gaponenko, 2014

#include "HistMuCapAnalysisChannels.h"

#include <stdexcept>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapAnalysisChannels::init(HistogramFactory& hf,
                                     const std::string& htopdir,
                                     const std::string& channelsetname,
                                     const DetectorGeo& geom,
                                     const ConfigFile& conf)
{
  const std::string hdir = htopdir + "/" + channelsetname;

  doMCTruth_ = conf.read<bool>("TruthBank/Do");
  fillExtras_TDC_ = conf.read<bool>(hdir+"/fillExtras_TDC", false);

  if(doMCTruth_) {
    refsample_muminus_multiplicity_ = hf.DefineTH1D(hdir+"/refsample", "muminus_multiplicity", "Mu- multiplicity, input", 10, -0.5, 9.5);
    refsample_endvtx_time_ = hf.DefineTH1D(hdir+"/refsample", "endvtx_time", "Candidate vtx time", 400, -500., 3500.);
    refsample_num_stops_ = hf.DefineTH1D(hdir+"/refsample", "num_stops", "Num selected stops per event", 10, -0.5, 9.5);

    refsample_in_zstop_ = hf.DefineTH1D(hdir+"/refsample", "in_zstop", "Z stop position, input events", 120, -0.0200, -0.0080);
    refsample_accepted_count_ = hf.DefineTH1D(hdir+"/refsample", "acc_count", "Count of accepted ref sample events", 1, -0.5, 0.5);
  }

  //----------------------------------------------------------------
  const std::string containedVars = conf.read<std::string>(hdir+"/contained/vars");
  if(containedVars == "RangeCosVsP") {
    // leak the memory - no elegant solution here without C++11
    contained_.init(hf, htopdir, channelsetname, geom, conf, *new MuCapContainedVars::RangeCosVsP());
  }
  else if(containedVars == "RangeVsP") {
    contained_.init(hf, htopdir, channelsetname, geom, conf, *new MuCapContainedVars::RangeVsP());
  }
  else if(containedVars == "RangeVsPz") {
    contained_.init(hf, htopdir, channelsetname, geom, conf, *new MuCapContainedVars::RangeVsPz());
  }
  else {
    throw std::runtime_error("HistMuCapAnalysisChannels::init(): unknow set of contained vars \""+containedVars+"\"");
  }

  uncontained_.init(hf, htopdir, channelsetname, geom, conf);

  hitbased_.init(hf, htopdir, channelsetname, geom, conf);

  if(doMCTruth_) {
    const int gen1nbins = conf.read<int>(hdir+"/numGeneratorBins");
    const double gen1pmin = conf.read<double>(hdir+"/genpmin");
    const double gen1pmax = conf.read<double>(hdir+"/genpmax");
    // truth level binning must be consistent for all channels
    mcin_proton_ptot_ = hf.DefineTH1D(hdir, "mcin_proton_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcin_deuteron_ptot_ = hf.DefineTH1D(hdir, "mcin_deuteron_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcin_dio_count_ = hf.DefineTH1D(hdir, "mcin_dio_count", "noncapture count, input", 1, -0.5, 0.5);
  }

  //----------------------------------------------------------------
  // Extra distributions

  if(fillExtras_TDC_) {
    hTDCWidthContained_.init(hf, hdir+"/contained/tdcwidth", geom, conf);
    hTDCWidthUncontained_.init(hf, hdir+"/uncontained/tdcwidth", geom, conf);
    hTDCWidthHitbased_.init(hf, hdir+"/hitbased/tdcwidth", geom, conf);
    hTDCWidthNone_.init(hf, hdir+"/nochannel/tdcwidth", geom, conf);
  }

  if(doMCTruth_) {
    hTruthContained_.init(hf, hdir+"/contained/MCTruth", conf);
    hTruthUncontained_.init(hf, hdir+"/uncontained/MCTruth", conf);
    hTruthHitbased_.init(hf, hdir+"/hitbased/MCTruth", conf);
    hTruthNone_.init(hf, hdir+"/nochannel/MCTruth", conf);

    hResolutionContained_.init(hf, hdir+"/contained/resolution", conf);
    hResolutionUncontained_.init(hf, hdir+"/uncontained/resolution", conf);
  }

}

//================================================================
namespace {
  int getFirstMCVertexIndexForTrack(const EventClass& evt, int imctrk) {
    int res = 0;
    for(unsigned i = 0; i < imctrk; ++i) {
      res += evt.mctrack_nv[i];
    }
    return res;
  }
}

void HistMuCapAnalysisChannels::fillReferenceSample(const EventClass& evt) {
  referenceSampleAccepted_ = false;
  if(doMCTruth_) {
    referenceSample_nrun_ = evt.nrun;
    referenceSample_nevt_ = evt.nevt;

    // Look for the trigger muon stop
    int numInputMuons = 0;
    int numMuStopCandidates = 0;
    int iMuStopTrack = -1;
    int iMuStopVtxEnd = -1;
    int iMuStopVtxStart = -1;

    for(unsigned i=0; i<evt.nmctr; ++i) {
      if(evt.mctrack_pid[i] == MuCapUtilities::PID_G3_MUMINUS) {
        ++numInputMuons;

        // Look at the end vertex of the muon track
        const int itmpvtxstart = getFirstMCVertexIndexForTrack(evt, i);
        const int itmpvtxend = itmpvtxstart + evt.mctrack_nv[i] - 1;
        const double stoptime = evt.mcvertex_time[itmpvtxend];
        refsample_endvtx_time_->Fill(stoptime);
        if(std::abs(stoptime) < 100.) {
          ++numMuStopCandidates;
          if(iMuStopTrack == -1) {
            iMuStopTrack = i;
          }
        }
      }
    }
    if(iMuStopTrack != -1) {
      // Set vertex indexes
      iMuStopVtxStart = getFirstMCVertexIndexForTrack(evt, iMuStopTrack);
      iMuStopVtxEnd = iMuStopVtxStart + evt.mctrack_nv[iMuStopTrack] - 1;
    }

    refsample_muminus_multiplicity_->Fill(numInputMuons);
    refsample_num_stops_->Fill(numMuStopCandidates);

    // Apply reference sample cuts
    if(iMuStopVtxEnd != -1) {
      const double zstop = evt.mcvertex_vz[iMuStopVtxEnd];
      refsample_in_zstop_->Fill(zstop);
      // From dt_geo.00066:
      // tgt mylar foil at -0.00904 cm, thickness 0.0025 cm
      // tgt Al  0.0071 cm thick
      static const double zmax = -0.00904 - 0.0025/2;
      static const double zmin = zmax - 0.0071;
      referenceSampleAccepted_ = (zmin<=zstop) && (zstop<=zmax);

      if(referenceSampleAccepted_) {
        refsample_accepted_count_->Fill(0.);

        // Truth momentum with the binning used in the unfolding
        // Keep protons and deuterons separately to compare hadd-ed
        // pseudodata truth to unfolding results.

        const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
        // Simulated DIO have no easily accessible MC truth.  We'll tread PID=zero as DIO down in this code.
        const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

        switch(mcParticle) {
        case MuCapUtilities::PID_G3_PROTON:
          mcin_proton_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
          break;
        case MuCapUtilities::PID_G3_DEUTERON:
          mcin_deuteron_ptot_->Fill(evt.mcvertex_ptot[imcvtxStart]);
          break;
        case 0:
          mcin_dio_count_->Fill(0.);
          break;
        }
      }
    }
    else {
      std::ostringstream os;
      os<<"HistMuCapAnalysisChannels::fillReferenceSample(): no mu stop vtx in run "<<evt.nrun<<", event "<<evt.nevt;
      //throw std::runtime_error(os.str());
      std::cout<<os.str()<<std::endl;
    }
  }
}

//================================================================
void HistMuCapAnalysisChannels::fill(const EventClass& evt,
                                     int iPosTrack,
                                     int iNegTrack,
                                     const ClustersByPlane& globalPlaneClusters )
{
  if(doMCTruth_ && ((referenceSample_nrun_ != evt.nrun)||(referenceSample_nevt_ != evt.nevt))) {
    throw std::runtime_error("Error: HistMuCapAnalysisChannels::fill() is called before fillReferenceSample() on that event.");
  }

  //----------------------------------------------------------------
  // Figure out an exclusive analysis channel for this event

  enum Channel { CONTAINED, UNCONTAINED, HITBASED, NONE};
  const Channel ch =
    contained_.accepted(evt, referenceSampleAccepted_, iPosTrack, iNegTrack, globalPlaneClusters) ?
    CONTAINED
    : ( uncontained_.accepted(evt, referenceSampleAccepted_, iPosTrack, iNegTrack) ?
        UNCONTAINED
        : ( hitbased_.accepted(evt, globalPlaneClusters, iNegTrack, referenceSampleAccepted_) ?
            HITBASED
            : NONE
            )
        );

  // Fill extra distributions
  switch(ch) {
  case CONTAINED: if(fillExtras_TDC_) { hTDCWidthContained_.fill(evt, globalPlaneClusters); }
    break;

  case UNCONTAINED: if(fillExtras_TDC_) { hTDCWidthUncontained_.fill(evt, globalPlaneClusters); }
    break;

  case HITBASED: if(fillExtras_TDC_) { hTDCWidthHitbased_.fill(evt, globalPlaneClusters); }
    break;

  case NONE: if(fillExtras_TDC_) { hTDCWidthNone_.fill(evt, globalPlaneClusters); }
    break;
  }

  if(doMCTruth_) {
    switch(ch) {
    case CONTAINED:
      hTruthContained_.fill(evt);
      hResolutionContained_.fill(evt, iPosTrack);
      break;

    case UNCONTAINED:
      hTruthUncontained_.fill(evt);
      hResolutionUncontained_.fill(evt, iPosTrack);
      break;

    case HITBASED:
      hTruthHitbased_.fill(evt);
      break;

    case NONE:
      hTruthNone_.fill(evt);
      break;
    }
  }
}

//================================================================
