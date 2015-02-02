// Andrei Gaponenko, 2014

#include "HistMuCapAnalysisChannels.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapAnalysisChannels::init(HistogramFactory& hf,
                                     const std::string& hdir,
                                     const DetectorGeo& geom,
                                     const ConfigFile& conf)
{
  doMCTruth_ = conf.read<bool>("TruthBank/Do");
  if(doMCTruth_) {
    refsample_muminus_multiplicity_ = hf.DefineTH1D(hdir+"/refsample", "muminus_multiplicity", "Mu- multiplicity, input", 10, -0.5, 9.5);
    refsample_endvtx_time_ = hf.DefineTH1D(hdir+"/refsample", "endvtx_time", "Candidate vtx time", 400, -500., 3500.);
    refsample_num_stops_ = hf.DefineTH1D(hdir+"/refsample", "num_stops", "Num selected stops per event", 10, -0.5, 9.5);

    refsample_in_zstop_ = hf.DefineTH1D(hdir+"/refsample", "in_zstop", "Z stop position, input events", 120, -0.0200, -0.0080);
    refsample_accepted_count_ = hf.DefineTH1D(hdir+"/refsample", "acc_count", "Count of accepted ref sample events", 1, -0.5, 0.5);
  }

  //----------------------------------------------------------------
  // "channel" analysis histograms
  const int recoContPNbins = 88;
  const double recoContPMin = 30.;
  const double recoContPMax = 250.;
  const int recoContRNbins = 6;
  const double recoContRMin = 8.;
  const double recoContRMax = 32.;

  const int recoUncPNbins = 88;
  const double recoUncPMin = 30.;
  const double recoUncPMax = 250.;

  //----------------------------------------------------------------
  contained_prange_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP", "Last plane hit/|cos(theta)| vs p, contained",
                                    recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
  contained_prange_->SetOption("colz");
  contained_prange_->GetXaxis()->SetTitle("p [MeV/c]");
  contained_prange_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

  if(doMCTruth_) {
    contained_prange_mcproton_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP_mcproton", "Last plane hit/|cos(theta)| vs p, contained mcproton",
                                               recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
    contained_prange_mcproton_->SetOption("colz");
    contained_prange_mcproton_->GetXaxis()->SetTitle("p [MeV/c]");
    contained_prange_mcproton_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

    contained_prange_mcdeuteron_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP_mcdeuteron", "Last plane hit/|cos(theta)| vs p, contained mcdeuteron",
                                                 recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
    contained_prange_mcdeuteron_->SetOption("colz");
    contained_prange_mcdeuteron_->GetXaxis()->SetTitle("p [MeV/c]");
    contained_prange_mcdeuteron_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");

    contained_prange_mcdio_ = hf.DefineTH2D(hdir+"/contained", "rangecosVsP_mcdio", "Last plane hit/|cos(theta)| vs p, contained mcdio",
                                            recoContPNbins, recoContPMin, recoContPMax, recoContRNbins, recoContRMin, recoContRMax);
    contained_prange_mcdio_->SetOption("colz");
    contained_prange_mcdio_->GetXaxis()->SetTitle("p [MeV/c]");
    contained_prange_mcdio_->GetYaxis()->SetTitle("(plane-28)/|cos(theta)|");
  }

  //----------------------------------------------------------------
  uncontained_p_ = hf.DefineTH1D(hdir+"/uncontained", "ptot", "ptot, non contained", recoUncPNbins, recoUncPMin, recoUncPMax);
  uncontained_p_->GetXaxis()->SetTitle("p, MeV/c");

  if(doMCTruth_) {
    uncontained_p_mcproton_ = hf.DefineTH1D(hdir+"/uncontained", "ptot_mcproton", "ptot, non contained mcproton", recoUncPNbins, recoUncPMin, recoUncPMax);
    uncontained_p_mcproton_->GetXaxis()->SetTitle("p, MeV/c");

    uncontained_p_mcdeuteron_ = hf.DefineTH1D(hdir+"/uncontained", "ptot_mcdeuteron", "ptot, non contained mcdeuteron", recoUncPNbins, recoUncPMin, recoUncPMax);
    uncontained_p_mcdeuteron_->GetXaxis()->SetTitle("p, MeV/c");

    uncontained_p_mcdio_ = hf.DefineTH1D(hdir+"/uncontained", "ptot_mcdio", "ptot, non contained mcdio", recoUncPNbins, recoUncPMin, recoUncPMax);
    uncontained_p_mcdio_->GetXaxis()->SetTitle("p, MeV/c");
  }

  //----------------------------------------------------------------
  hitbased_.init(hf, hdir+"/hitbased", geom, conf);

  hTDCWidthContained_.init(hf, hdir+"/contained/tdcwidth", geom, conf);
  hTDCWidthUncontained_.init(hf, hdir+"/uncontained/tdcwidth", geom, conf);

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning
    const int gen1nbins = 400; // 2.5 MeV/c bins
    const double gen1pmin = 0.;
    const double gen1pmax = 400.;

    // True distribution of all events used for the unfolding
    mcin_proton_ptot_ = hf.DefineTH1D(hdir, "mcin_proton_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcin_deuteron_ptot_ = hf.DefineTH1D(hdir, "mcin_deuteron_ptot", "mcptot, input", gen1nbins, gen1pmin, gen1pmax);
    mcin_dio_count_ = hf.DefineTH1D(hdir, "mcin_dio_count", "noncapture count, input", 1, -0.5, 0.5);

    // Migration matrices for the contained channel
    containedMigration_ = hf.DefineTH3D(hdir+"/contained", "migration",
                                        "Contained channel migration",
                                        gen1nbins, gen1pmin, gen1pmax,
                                        recoContPNbins, recoContPMin, recoContPMax,
                                        recoContRNbins, recoContRMin, recoContRMax);

    containedMigration_->GetXaxis()->SetTitle("p true, MeV/c");
    containedMigration_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedMigration_->GetZaxis()->SetTitle("range");

    containedMigration_mcproton_ = hf.DefineTH3D(hdir+"/contained", "migration_mcproton",
                                                 "Contained channel migration, proton",
                                                 gen1nbins, gen1pmin, gen1pmax,
                                                 recoContPNbins, recoContPMin, recoContPMax,
                                                 recoContRNbins, recoContRMin, recoContRMax);

    containedMigration_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    containedMigration_mcproton_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedMigration_mcproton_->GetZaxis()->SetTitle("range");

    containedMigration_mcdeuteron_ = hf.DefineTH3D(hdir+"/contained", "migration_mcdeuteron",
                                                 "Contained channel migration, deuteron",
                                                 gen1nbins, gen1pmin, gen1pmax,
                                                 recoContPNbins, recoContPMin, recoContPMax,
                                                 recoContRNbins, recoContRMin, recoContRMax);

    containedMigration_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    containedMigration_mcdeuteron_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedMigration_mcdeuteron_->GetZaxis()->SetTitle("range");

    // Migration matrices for the non contained channel
    uncontainedMigration_ = hf.DefineTH2D(hdir+"/uncontained","migration",
                                          "Non-contained channel migration",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedMigration_->SetOption("colz");
    uncontainedMigration_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedMigration_->GetYaxis()->SetTitle("p reco, MeV/c");

    uncontainedMigration_mcproton_ = hf.DefineTH2D(hdir+"/uncontained","migration_mcproton",
                                                   "Non-contained channel migration, proton",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedMigration_mcproton_->SetOption("colz");
    uncontainedMigration_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedMigration_mcproton_->GetYaxis()->SetTitle("p reco, MeV/c");

    uncontainedMigration_mcdeuteron_ = hf.DefineTH2D(hdir+"/uncontained","migration_mcdeuteron",
                                                   "Non-contained channel migration, deuteron",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedMigration_mcdeuteron_->SetOption("colz");
    uncontainedMigration_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedMigration_mcdeuteron_->GetYaxis()->SetTitle("p reco, MeV/c");

    // Contamination matrices for the contained channel
    containedContamination_ = hf.DefineTH3D(hdir+"/contained", "contamination",
                                        "Contained channel contamination",
                                        gen1nbins, gen1pmin, gen1pmax,
                                        recoContPNbins, recoContPMin, recoContPMax,
                                        recoContRNbins, recoContRMin, recoContRMax);

    containedContamination_->GetXaxis()->SetTitle("p true, MeV/c");
    containedContamination_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedContamination_->GetZaxis()->SetTitle("range");

    containedContamination_mcproton_ = hf.DefineTH3D(hdir+"/contained", "contamination_mcproton",
                                                 "Contained channel contamination, proton",
                                                 gen1nbins, gen1pmin, gen1pmax,
                                                 recoContPNbins, recoContPMin, recoContPMax,
                                                 recoContRNbins, recoContRMin, recoContRMax);

    containedContamination_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    containedContamination_mcproton_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedContamination_mcproton_->GetZaxis()->SetTitle("range");

    containedContamination_mcdeuteron_ = hf.DefineTH3D(hdir+"/contained", "contamination_mcdeuteron",
                                                 "Contained channel contamination, deuteron",
                                                 gen1nbins, gen1pmin, gen1pmax,
                                                 recoContPNbins, recoContPMin, recoContPMax,
                                                 recoContRNbins, recoContRMin, recoContRMax);

    containedContamination_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    containedContamination_mcdeuteron_->GetYaxis()->SetTitle("p reco, MeV/c");
    containedContamination_mcdeuteron_->GetZaxis()->SetTitle("range");

    // Contamination matrices for the non contained channel
    uncontainedContamination_ = hf.DefineTH2D(hdir+"/uncontained","contamination",
                                          "Non-contained channel contamination",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedContamination_->SetOption("colz");
    uncontainedContamination_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedContamination_->GetYaxis()->SetTitle("p reco, MeV/c");

    uncontainedContamination_mcproton_ = hf.DefineTH2D(hdir+"/uncontained","contamination_mcproton",
                                                   "Non-contained channel contamination, proton",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedContamination_mcproton_->SetOption("colz");
    uncontainedContamination_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedContamination_mcproton_->GetYaxis()->SetTitle("p reco, MeV/c");

    uncontainedContamination_mcdeuteron_ = hf.DefineTH2D(hdir+"/uncontained","contamination_mcdeuteron",
                                                   "Non-contained channel contamination, deuteron",
                                                   gen1nbins, gen1pmin, gen1pmax,
                                                   recoUncPNbins, recoUncPMin, recoUncPMax);

    uncontainedContamination_mcdeuteron_->SetOption("colz");
    uncontainedContamination_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    uncontainedContamination_mcdeuteron_->GetYaxis()->SetTitle("p reco, MeV/c");

    //----------------------------------------------------------------
    hTruthTrkContained_.init(hf, hdir+"/contained/MCTruthTrk", conf);
    hTruthTrkUncontained_.init(hf, hdir+"/uncontained/MCTruthTrk", conf);

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
      throw std::runtime_error(os.str());
    }
  }
}

//================================================================
void HistMuCapAnalysisChannels::fill(const EventClass& evt,
                                     int iPosTrack,
                                     int iNegTrack,
                                     bool isPosTrackContained,
                                     double rangePIDVar,
                                     const ClustersByPlane& globalPlaneClusters
                                     )
{
  if(doMCTruth_ && ((referenceSample_nrun_ != evt.nrun)||(referenceSample_nevt_ != evt.nevt))) {
    throw std::runtime_error("Error: HistMuCapAnalysisChannels::fill() is called before fillReferenceSample() on that event.");
  }

  //----------------------------------------------------------------
  // Figure out an exclusive analysis channel for this event

  bool eventUsedInAChannel2 = false;

  const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
  // Simulated DIO have no easily accessible MC truth.  We'll tread PID=zero as DIO down in this code.
  const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

  if(iNegTrack == -1) { // Veto DIO events
    if(iPosTrack != -1) { // Got a reconstructed capture track

      const double prec = evt.ptot[iPosTrack];

      if(isPosTrackContained) {
        // The "contained tracks" analysis channel
        eventUsedInAChannel2 = true;
        contained_prange_->Fill(prec, rangePIDVar);
        hTDCWidthContained_.fill(evt, globalPlaneClusters);

        if(doMCTruth_) {

          if(imcvtxStart  != -1) {
            (referenceSampleAccepted_ ? containedMigration_ : containedContamination_)
              ->Fill(evt.mcvertex_ptot[imcvtxStart], prec, rangePIDVar);
          }

          switch(mcParticle) {
          case MuCapUtilities::PID_G3_PROTON:
            contained_prange_mcproton_->Fill(prec, rangePIDVar);
            (referenceSampleAccepted_ ? containedMigration_mcproton_ : containedContamination_mcproton_)
              ->Fill(evt.mcvertex_ptot[imcvtxStart], prec, rangePIDVar);
            break;
          case MuCapUtilities::PID_G3_DEUTERON:
            contained_prange_mcdeuteron_->Fill(prec, rangePIDVar);
            (referenceSampleAccepted_ ? containedMigration_mcdeuteron_ : containedContamination_mcdeuteron_)
              ->Fill(evt.mcvertex_ptot[imcvtxStart], prec, rangePIDVar);
            break;
          case 0:
            contained_prange_mcdio_->Fill(prec, rangePIDVar);
            break;
          }

          hTruthTrkContained_.fill(evt);
          hResolutionContained_.fill(evt, iPosTrack);
        }
      }
      else { // The non-contained tracks channel
        eventUsedInAChannel2 = true;
        uncontained_p_->Fill(prec);

        hTDCWidthUncontained_.fill(evt, globalPlaneClusters);

        if(doMCTruth_) {

          if(imcvtxStart  != -1) {
            (referenceSampleAccepted_ ? uncontainedMigration_ : uncontainedContamination_)
              ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
          }

          switch(mcParticle) {
          case MuCapUtilities::PID_G3_PROTON:
            uncontained_p_mcproton_->Fill(prec);
            (referenceSampleAccepted_ ? uncontainedMigration_mcproton_ : uncontainedContamination_mcproton_)
              ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
            break;
          case MuCapUtilities::PID_G3_DEUTERON:
            uncontained_p_mcdeuteron_->Fill(prec);
            (referenceSampleAccepted_ ? uncontainedMigration_mcdeuteron_ : uncontainedContamination_mcdeuteron_)
              ->Fill(evt.mcvertex_ptot[imcvtxStart], prec);
            break;
          case 0:
            uncontained_p_mcdio_->Fill(prec);
            break;
          }

          hTruthTrkUncontained_.fill(evt);
          hResolutionUncontained_.fill(evt, iPosTrack);
        }
      }
    }
  }

  if(!eventUsedInAChannel2) {
    // No good positive tracks. Do a hit-based analysis here.
    hitbased_.accepted(evt, globalPlaneClusters, iNegTrack, referenceSampleAccepted_);
  }
}

//================================================================
