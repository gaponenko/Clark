// Wire plane TDC hit.
// Andrei Gaponenko, 2013

#ifndef TDCHitWP_h
#define TDCHitWP_h

#include <vector>

struct TDCHitWP {
  float time;
  float width;
  int plane;
  int cell;

  TDCHitWP(float t, float w, float p, float c)
    : time(t), width(w), plane(p), cell(c)
  {}
};

typedef std::vector<TDCHitWP> TDCHitWPCollection;

//=====
// A comparator, for use with std::sort

struct TDCHitWPCmpTime {
  bool operator()(const TDCHitWP& a, const TDCHitWP& b) {
    return a.time < b.time;
  }
};

//================================================================
// A non-owning pointer-like class for presenting different views
// on a hit collection

class TDCHitWPPtr {
  const TDCHitWPCollection *coll_;
  TDCHitWPCollection::size_type idx_;
public:
  TDCHitWPPtr() : coll_(0), idx_() {}
  TDCHitWPPtr(const TDCHitWPCollection& c, TDCHitWPCollection::size_type i)
    : coll_(&c), idx_(i)
  {}

  const TDCHitWP* operator->() const { return &(*coll_)[idx_]; }
  const TDCHitWP& operator*() const { return (*coll_)[idx_]; }
};

// Comparator
struct TDCHitWPCmpGeom {
  bool operator()(const TDCHitWP& a, const TDCHitWP& b) {
    return (a.plane < b.plane) ||
      (a.plane == b.plane) && (a.cell < b.cell);
  }

  bool operator()(const TDCHitWPPtr& a, const TDCHitWPPtr& b) {
    return (a->plane < b->plane) ||
      (a->plane == b->plane) && (a->cell < b->cell);
  }
};

#endif/*TDCHitWP_h*/
