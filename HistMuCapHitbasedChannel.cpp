// Andrei Gaponenko, 2014

#include "HistMuCapHitbasedChannel.h"
#include "HitBasedObservables.h"

#include <fstream>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "PlaneRange.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "MuCapUtilities.h"

//================================================================
void HistMuCapHitbasedChannel::init(HistogramFactory& hf,
                                    const std::string& hgrandtopdir,
                                    const std::string& channelsetname,
                                    const DetectorGeo& geom,
                                    const ConfigFile& conf)
{
  geom_ = &geom;

  doMCTruth_ = conf.read<bool>("TruthBank/Do");

  const std::string htopdir = hgrandtopdir+"/"+channelsetname;
  const std::string hdir = htopdir+"/hitbased";

  tdcWidthFilterCutPC_ = conf.read<double>("HitBasedAnalysis/tdcWidthFilterCutPC");
  // Art's suggestion: discard clusters that are too large
  maxClusterWiresFilterCutPC_ = conf.read<int>("HitBasedAnalysis/maxClusterWiresFilterCutPC");

  // Binning for the principal analysis channel histogram
  const int recoCWiresNBins = conf.read<int>(hdir+"/cwiresnbins");
  const double recoCWiresXMin = conf.read<double>(hdir+"/cwiresmin");
  const double recoCWiresXMax = conf.read<double>(hdir+"/cwiresmax");

  const int recoCPlanesNBins = conf.read<int>(hdir+"/cplanesnbins");
  const double recoCPlanesXMin = conf.read<double>(hdir+"/cplanesmin"); 
  const double recoCPlanesXMax = conf.read<double>(hdir+"/cplanesmax");

  //----------------------------------------------------------------
  h_cuts_r = hf.DefineTH1D(hdir, "cuts_r", "Events rejected by cut", CUTS_END, -0.5, CUTS_END-0.5);
  h_cuts_r->SetStats(kFALSE);
  set_cut_bin_labels(h_cuts_r->GetXaxis());
  h_cuts_r->SetOption("hist text");

  h_cuts_p = hf.DefineTH1D(hdir, "cuts_p", "Events before cut", CUTS_END, -0.5, CUTS_END-0.5);
  set_cut_bin_labels(h_cuts_p->GetXaxis());
  h_cuts_p->SetStats(kFALSE);
  h_cuts_p->SetOption("hist text");

  //----------------------------------------------------------------
  cutpc78_.init(hf, hdir+"/pc78cut", geom, conf);

  //----------------------------------------------------------------
  lastconPlaneVsCWires_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires", "Last contiguous plane vs sum(largest cluster size)",
                                        recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

  lastconPlaneVsCWires_->SetOption("colz");
  lastconPlaneVsCWires_->GetXaxis()->SetTitle("cluster wires");
  lastconPlaneVsCWires_->GetYaxis()->SetTitle("cplane");

  if(doMCTruth_) {
    lastconPlaneVsCWires_mcproton_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires_mcproton", "Last contiguous plane vs sum(largest cluster size), mcproton",
                                                   recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);
    lastconPlaneVsCWires_mcproton_->SetOption("colz");
    lastconPlaneVsCWires_mcproton_->GetXaxis()->SetTitle("cluster wires");
    lastconPlaneVsCWires_mcproton_->GetYaxis()->SetTitle("cplane");

    lastconPlaneVsCWires_mcdeuteron_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires_mcdeuteron", "Last contiguous plane vs sum(largest cluster size), mcdeuteron",
                                                     recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);
    lastconPlaneVsCWires_mcdeuteron_->SetOption("colz");
    lastconPlaneVsCWires_mcdeuteron_->GetXaxis()->SetTitle("cluster wires");
    lastconPlaneVsCWires_mcdeuteron_->GetYaxis()->SetTitle("cplane");

    lastconPlaneVsCWires_mcdio_ = hf.DefineTH2D(hdir, "cplanes_vs_cwires_mcdio", "Last contiguous plane vs sum(largest cluster size), mcdio",
                                                recoCWiresNBins, recoCWiresXMin, recoCWiresXMax, recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);
    lastconPlaneVsCWires_mcdio_->SetOption("colz");
    lastconPlaneVsCWires_mcdio_->GetXaxis()->SetTitle("cluster wires");
    lastconPlaneVsCWires_mcdio_->GetYaxis()->SetTitle("cplane");
  }

  hambig_.init(hf, hdir+"/clusterAmbiguities", geom, conf);

  //----------------------------------------------------------------
  if(doMCTruth_) {
    // truth level binning must be consistent for all channels
    const int gen1nbins = conf.read<int>(htopdir+"/numGeneratorBins");
    const double gen1pmin = conf.read<double>(htopdir+"/genpmin");
    const double gen1pmax = conf.read<double>(htopdir+"/genpmax");

    // Migration matrices for the contained channel, two generator binnings
    migration_ = hf.DefineTH3D(hdir, "migration", "Hit based channel migration",
                               gen1nbins, gen1pmin, gen1pmax,
                               recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                               recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    migration_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_->GetYaxis()->SetTitle("cluster wires");
    migration_->GetZaxis()->SetTitle("cplane");

    migration_mcproton_ = hf.DefineTH3D(hdir, "migration_mcproton", "Hit based channel migration, proton",
                                        gen1nbins, gen1pmin, gen1pmax,
                                        recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                                        recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    migration_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_mcproton_->GetYaxis()->SetTitle("cluster wires");
    migration_mcproton_->GetZaxis()->SetTitle("cplane");

    migration_mcdeuteron_ = hf.DefineTH3D(hdir, "migration_mcdeuteron", "Hit based channel migration, deuteron",
                                          gen1nbins, gen1pmin, gen1pmax,
                                          recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                                          recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    migration_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    migration_mcdeuteron_->GetYaxis()->SetTitle("cluster wires");
    migration_mcdeuteron_->GetZaxis()->SetTitle("cplane");

    // Contamination matrices for the contained channel, two generator binnings
    contamination_ = hf.DefineTH3D(hdir, "contamination", "Hit based channel contamination",
                                   gen1nbins, gen1pmin, gen1pmax,
                                   recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                                   recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    contamination_->GetXaxis()->SetTitle("p true, MeV/c");
    contamination_->GetYaxis()->SetTitle("cluster wires");
    contamination_->GetZaxis()->SetTitle("cplane");

    contamination_mcproton_ = hf.DefineTH3D(hdir, "contamination_mcproton", "Hit based channel contamination, proton",
                                            gen1nbins, gen1pmin, gen1pmax,
                                            recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                                            recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    contamination_mcproton_->GetXaxis()->SetTitle("p true, MeV/c");
    contamination_mcproton_->GetYaxis()->SetTitle("cluster wires");
    contamination_mcproton_->GetZaxis()->SetTitle("cplane");

    contamination_mcdeuteron_ = hf.DefineTH3D(hdir, "contamination_mcdeuteron", "Hit based channel contamination, deuteron",
                                              gen1nbins, gen1pmin, gen1pmax,
                                              recoCWiresNBins, recoCWiresXMin, recoCWiresXMax,
                                              recoCPlanesNBins, recoCPlanesXMin, recoCPlanesXMax);

    contamination_mcdeuteron_->GetXaxis()->SetTitle("p true, MeV/c");
    contamination_mcdeuteron_->GetYaxis()->SetTitle("cluster wires");
    contamination_mcdeuteron_->GetZaxis()->SetTitle("cplane");

    //----------------
    hTruth_in_.init(hf, hdir+"/MCTruth_in", conf);
    hTruth_accepted_.init(hf, hdir+"/MCTruth_accepted", conf);
  }

  hOuterVetoNumHitPlanes_ = hf.DefineTH1D(hdir, "outerVetoNumHitPlanes", "outerVetoNumHitPlanes", 6, -0.5, 5.5);
  hNumPC7Clusters_ = hf.DefineTH1D(hdir, "numPC7Clusters", "numPC7Clusters", 20, -0.5, 19.5);

  hFilterEffectPC7_ = hf.DefineTH2D(hdir, "hitFilterPC7", "Hit filter effect on PC7", 51, -25.5, +25.5, 50, -49.5, 0.5);
  hFilterEffectPC7_->SetOption("colz");
  hFilterEffectPC7_->GetXaxis()->SetTitle("nclusters filtered - orig");
  hFilterEffectPC7_->GetYaxis()->SetTitle("nwires filtered - orig");

  hFilterEffectPC8_ = hf.DefineTH2D(hdir, "hitFilterPC8", "Hit filter effect on PC8", 51, -25.5, +25.5, 50, -49.5, 0.5);
  hFilterEffectPC8_->SetOption("colz");
  hFilterEffectPC8_->GetXaxis()->SetTitle("nclusters filtered - orig");
  hFilterEffectPC8_->GetYaxis()->SetTitle("nwires filtered - orig");

  //----------------------------------------------------------------
  htdcwidthInput_.init(hf, hdir+"/tdcwidthInput", geom, conf);
//mem:  htdcwidthDoubleFiltered_.init(hf, hdir+"/tdcwidthDoubleFiltered", geom, conf);

//mem:  hshot_.init(hf, hdir+"/hot", geom, conf);
//mem:  hscold_.init(hf, hdir+"/cold", geom, conf);

  // This histogram is only iteresting for real data.  Here is our run number range:
  const int firstrun = 47607;
  const int lastrun =  47798;
  hshotrun_ = hf.DefineTH1D(hdir, "hshotrun", "Run number for hot spot events", 1+lastrun-firstrun, firstrun-0.5, lastrun+0.5);

//mem:  hxtplane100_.init(hf, hdir+"/xtplane100", geom, conf, 100.);
//mem:  hxtplane300_.init(hf, hdir+"/xtplane300", geom, conf, 300.);
  //----------------------------------------------------------------
}

//================================================================
bool HistMuCapHitbasedChannel::accepted(const EventClass& evt, const ClustersByPlane& protonGlobalClusters, int iDIOVetoTrack, bool referenceSampleAccepted) {
  CutNumber c = analyzeEvent(evt, protonGlobalClusters, iDIOVetoTrack, referenceSampleAccepted);
  h_cuts_r->Fill(c);
  for(int cut=0; cut<=c; cut++) {
    h_cuts_p->Fill(cut);
  }
  return c == CUTS_ACCEPTED;
}

//================================================================
HistMuCapHitbasedChannel::CutNumber HistMuCapHitbasedChannel::analyzeEvent(const EventClass& evt, const ClustersByPlane& inputClusters, int iDIOVetoTrack, bool referenceSampleAccepted) {
  if(doMCTruth_) {
    hTruth_in_.fill(evt);
  }

  if(iDIOVetoTrack != -1) {
    return CUT_DIOVETO;
  }

  // Filter out noise hits in downstream PCs.
  // We accept the loss of some real electron hits at this stage.
  ClustersByPlane doubleFilteredClusters;
  filterDnPCNoise(&doubleFilteredClusters, inputClusters);

  // Z containment cut
  int numOuterVetoHitPlanes(0);
  for(int i=0; i<4; ++i) {
    if(!doubleFilteredClusters[56-i].empty()) {
      ++numOuterVetoHitPlanes;
    }
  }

  hOuterVetoNumHitPlanes_->Fill(numOuterVetoHitPlanes);
  if(numOuterVetoHitPlanes > 1) {
    return CUT_ZVETO;
  }

  const unsigned numPC7Clusters = doubleFilteredClusters.at(29).size();
  hNumPC7Clusters_->Fill(numPC7Clusters);
  if(!numPC7Clusters) {
    return CUT_NOPC7;
  }

  if(!cutpc78_.accepted(evt, doubleFilteredClusters)) {
    return CUT_PCWIDTH;
  }

  //----------------------------------------------------------------
  HitBasedObservablesMaxWidth obs(doubleFilteredClusters, &hambig_);

  lastconPlaneVsCWires_->Fill(obs.dnCWires(), obs.dnCPlanes());

  if(doMCTruth_) {
    const unsigned imcvtxStart = evt.iCaptureMcVtxStart;
    // Simulated DIO have no easily accessible MC truth.  We'll treat PID=zero as DIO down in this code.
    const int mcParticle = (imcvtxStart != -1) ? evt.mctrack_pid[evt.iCaptureMcTrk] : 0;

    if(imcvtxStart  != -1) {
      (referenceSampleAccepted ? migration_ : contamination_)->Fill(evt.mcvertex_ptot[imcvtxStart], obs.dnCWires(), obs.dnCPlanes());
    }

    hTruth_accepted_.fill(evt);

    switch(mcParticle) {
    case MuCapUtilities::PID_G3_PROTON:
      lastconPlaneVsCWires_mcproton_->Fill(obs.dnCWires(), obs.dnCPlanes());
      (referenceSampleAccepted ? migration_mcproton_ : contamination_mcproton_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], obs.dnCWires(), obs.dnCPlanes());
      break;
    case MuCapUtilities::PID_G3_DEUTERON:
      lastconPlaneVsCWires_mcdeuteron_->Fill(obs.dnCWires(), obs.dnCPlanes());
      (referenceSampleAccepted ? migration_mcdeuteron_ : contamination_mcdeuteron_)
        ->Fill(evt.mcvertex_ptot[imcvtxStart], obs.dnCWires(), obs.dnCPlanes());
      break;
    case 0:
      lastconPlaneVsCWires_mcdio_->Fill(obs.dnCWires(), obs.dnCPlanes());
      break;
    }
  }

  //----------------------------------------------------------------
  htdcwidthInput_.fill(evt, inputClusters);
  //mem:htdcwidthDoubleFiltered_.fill(evt, doubleFilteredClusters);

  //----------------
  if((obs.dnCPlanes() == 2)&&(obs.dnCWires()>50)) {
    //mem:hshot_.fill(evt, inputClusters, obs);
    hshotrun_->Fill(evt.nrun);
  }
  else {
    //mem:hscold_.fill(evt, inputClusters, obs);
  }

//mem:  hxtplane100_.fill(evt, inputClusters);
//mem:  hxtplane300_.fill(evt, inputClusters);

  //----------------------------------------------------------------

  return CUTS_ACCEPTED;
}

//================================================================
void HistMuCapHitbasedChannel::filterDnPCNoise(ClustersByPlane *out, const ClustersByPlane& in) {

  out->resize(in.size());

  TDCHitWPPtrCollection pchits;

  for(int iplane=1; iplane<in.size(); ++iplane) {
    // Filter hist from donwstream PCs into a new collection.
    // Copy over clusters for other planes;
    if( (iplane > 28) && (geom_->global(iplane).planeType() == WirePlane::PC)) {
      const WireClusterCollection& clusters = in[iplane];
      for(WireClusterCollection::const_iterator ic = clusters.begin(); ic != clusters.end(); ++ic) {
        for(TDCHitWPPtrCollection::const_iterator ih = ic->hits().begin(); ih != ic->hits().end(); ++ih) {
          if((*ih)->width() > tdcWidthFilterCutPC_) {
            pchits.push_back(*ih);
          }
        }
      }
    }
    else {
      (*out)[iplane] = in[iplane];
    }
  }

  // Combine filtered hits into clusters
  ClustersByPlane dnPCClusters = constructPlaneClusters(geom_->numPCs(), pchits);

  fillFilterEffectHist(hFilterEffectPC7_, in[29], dnPCClusters[7]);
  fillFilterEffectHist(hFilterEffectPC8_, in[30], dnPCClusters[8]);

  // Record the new clusters
  filterClusterSize( &(*out)[29], dnPCClusters[7]); // pc7
  filterClusterSize( &(*out)[30], dnPCClusters[8]); // pc8

  filterClusterSize( &(*out)[53], dnPCClusters[9]); // pc9
  filterClusterSize( &(*out)[54], dnPCClusters[10]);
  filterClusterSize( &(*out)[55], dnPCClusters[11]);
  filterClusterSize( &(*out)[56], dnPCClusters[12]);
}

//================================================================
void HistMuCapHitbasedChannel::filterClusterSize(WireClusterCollection *out, const WireClusterCollection& in) {
  out->clear();
  for(WireClusterCollection::const_iterator i=in.begin(); i!=in.end(); ++i) {
    if(i->numCells() <= maxClusterWiresFilterCutPC_) {
      out->push_back(*i);
    }
  }
}

//================================================================
void HistMuCapHitbasedChannel::fillFilterEffectHist(TH2* hh, const WireClusterCollection& orig, const WireClusterCollection& filtered) {
  int dcluster = filtered.size() - orig.size();

  int nworig=0;
  for(WireClusterCollection::const_iterator i=orig.begin(); i!=orig.end(); ++i) {
    nworig += i->numCells();
  }

  int nwfiltered=0;
  for(WireClusterCollection::const_iterator i=filtered.begin(); i!=filtered.end(); ++i) {
    nwfiltered += i->numCells();
  }

  int dnw = nwfiltered - nworig;

  hh->Fill(dcluster, dnw);
}

//================================================================
