// Look at events with an accepted downstream track.
// They guaranteed to have a real track hit in DC23.
// In XT4 we use events with wire multiplicity greater
// than what is kinenatically allowed for the track to
// identify potential cross talk hits.
//
// Andrei Gaponenko, 2017

#ifndef HistXT4_h
#define HistXT4_h

#include <string>
#include <vector>

#include "TDCHitWP.h"
#include "WireCluster.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class ConfigFile;
class EventClass;


//================================================================
// a helper class
class HistXT4WireGap {
public:
  void init(const std::string& hdir, unsigned cutWireGap, const std::string& suffix, HistogramFactory &hf, const ConfigFile &conf);
  void fill(const TDCHitWPPtrCollection& planeHits);
private :
  unsigned cutWireGap_;

  TH1 *outwide_dt_all_;
  TH2 *outwide_ww_all_;
  TH2 *outwide_wt_all_;

  TH1 *outwide_dt_xt_;
  TH2 *outwide_ww_xt_;
  TH2 *outwide_wt_xt_;

  TH1 *outwide_dt_nxt_;
  TH2 *outwide_ww_nxt_;
  TH2 *outwide_wt_nxt_;
};

//================================================================
class HistXT4 {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const EventClass& evt, int iTrack);

  HistXT4() {}

private :
  TH1 *trackCosth_;
  TH1 *trackAngle_;

  HistXT4WireGap wg1_;
  HistXT4WireGap wg2_;
};

#endif/*HistXT4_h*/
