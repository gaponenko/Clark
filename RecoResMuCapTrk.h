// A reconstruction summary for the track-based analysis,
// used to fill the RooUnfold matrix.
//
// Andrei Gaponenko, 2013

#ifndef RecoResMuCapTrk_h
#define RecoResMuCapTrk_h

//================================================================
struct RecoResMuCapTrk {
  bool accepted;
  double momentum;
  RecoResMuCapTrk() : accepted(false), momentum(-1.) {}
};

#endif/*RecoResMuCapTrk_h*/
