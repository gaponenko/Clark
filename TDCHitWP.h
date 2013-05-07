// Wire plane TDC hit.
// Andrei Gaponenko, 2013

#ifndef TDCHitWP_h
#define TDCHitWP_h

#include <vector>
#include <ostream>

//================================================================
struct WireCellId {
  int plane;
  int cell;
  WireCellId(int p, int c) : plane(p), cell(c) {}

  bool operator<(const WireCellId& b) const {
    return (plane < b.plane) ||
      (plane == b.plane) && (cell < b.cell);
  }
};

//================================================================
struct TDCHitWP {
  WireCellId cid;
  float time;
  float width;
  int plane() const { return cid.plane; }
  int cell() const { return cid.cell; }

  TDCHitWP(float t, float w, int p, int c)
    : cid(p,c), time(t), width(w)
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

typedef std::vector<TDCHitWPPtr> TDCHitWPPtrCollection;

//----------------------------------------------------------------
// Comparator
struct TDCHitWPCmpGeom {
  bool operator()(const TDCHitWP& a, const TDCHitWP& b) const {
    return (a.cid < b.cid);
  }

  bool operator()(const TDCHitWPPtr& a, const TDCHitWPPtr& b) const {
    return (a->cid < b->cid);
  }
};

//----------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const WireCellId& cid);
std::ostream& operator<<(std::ostream& os, const TDCHitWP& hit);
std::ostream& operator<<(std::ostream& os, const TDCHitWPCollection& hits);
std::ostream& operator<<(std::ostream& os, const TDCHitWPPtrCollection& hits);

#endif/*TDCHitWP_h*/
