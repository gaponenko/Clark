// Quantify the amount of residual cross talk in data.
//
// We use a sample of "good" positive (proton) tracks.
// We impose a 40 degree cut on track angle.
// Tracks under 45 degrees can not hit more than 2 wires.
// Also, a single track should make only a single wire cluster per plane.
// We also combine the two criteria (size and multiplicity)
// into a single boolean: "single track" vs "extra stuff".
//
// We look at cluster multiplicity and cluster size, and their
// correlation between DC23 and DC24.  Assuming that the  cross
// talk in different planes is not correlated,  its contribution to the
// (extra stuff, exra stuff)  bin is predictable from
// (as expected,extra stuff) and (extra stuff, as expected)
// relative to (as expected, as expected).
//
// On the other hand extra particles can make correlated contribution
// to the (extra stuff, extra stuff) bin.
//
//
// Andrei Gaponenko, 2017

#ifndef HistXT5_h
#define HistXT5_h

#include <string>
#include <vector>

#include "TDCHitWP.h"
#include "WireCluster.h"

#include "MuCapContainedVars.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;
class EventClass;

//================================================================
// a helper class
class HistXT5Ana {
public:
  void init(const std::string& hdir, const std::string& suffix, HistogramFactory &hf, const ConfigFile &conf);
  void fill(const TDCHitWPPtrCollection& planeHits);
private :
  TH2 *clustermult_;
  TH2 *clustersize_;

  TH2 *logic_cm_;
  TH2 *logic_cs_;
  TH2 *logic_cmcs_;
};

//================================================================
class HistXT5 {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const EventClass& evt, int iTrack);

  HistXT5() {}

private :
  TH2 *hpcos_all_;
  TH2 *hpcos_passed_;

  TH1 *trackAngle_all_;
  TH1 *trackAngle_passed_;

  HistXT5Ana xtanaAll_;
  HistXT5Ana xtanaNXT_;
};

#endif/*HistXT5_h*/
