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

  bool operator==(const WireCellId& b) const {
    return (plane == b.plane) && (cell == b.cell);
  }

  bool operator!=(const WireCellId& b) const {
    return !(*this == b);
  }
};

//================================================================
class TDCHitWP {
  WireCellId cid_;
  float time_;
  float width_;
  bool xtalk_;
public:
  int plane() const { return cid_.plane; }
  int cell() const { return cid_.cell; }
  float time() const { return time_; }
  float width() const { return width_; }
  bool xtalk() const { return xtalk_; }

  const WireCellId& cid() const { return cid_; }

  TDCHitWP(float t, float w, int p, int c, int xt)
    : cid_(p,c), time_(t), width_(w), xtalk_(xt)
  {}
};

typedef std::vector<TDCHitWP> TDCHitWPCollection;

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
// Comparators for use with std::sort

struct TDCHitWPCmpTime {
  bool operator()(const TDCHitWP& a, const TDCHitWP& b) {
    return a.time() < b.time();
  }

  bool operator()(const TDCHitWPPtr& a, const TDCHitWPPtr& b) const {
    return (a->time() < b->time());
  }
};


struct TDCHitWPCmpGeom {
  bool operator()(const TDCHitWP& a, const TDCHitWP& b) const {
    return (a.cid() < b.cid());
  }

  bool operator()(const TDCHitWPPtr& a, const TDCHitWPPtr& b) const {
    return (a->cid() < b->cid());
  }
};

//----------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const WireCellId& cid);
std::ostream& operator<<(std::ostream& os, const TDCHitWP& hit);
std::ostream& operator<<(std::ostream& os, const TDCHitWPCollection& hits);
std::ostream& operator<<(std::ostream& os, const TDCHitWPPtrCollection& hits);

#endif/*TDCHitWP_h*/
