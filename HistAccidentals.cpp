// Andrei Gaponenko, 2013

#include "HistAccidentals.h"
#include <sstream>
#include <vector>

#include "TH1.h"

#include "HistogramFactory.h"
#include "ConfigFile.h"
#include "TimeWindow.h"

//================================================================
namespace {
  // essentially a histogram
  class Subdivision {
    double xmin_;
    double xmax_;
    std::vector<int> bins_;
  public:
    Subdivision(int n, double xmin, double xmax) : xmin_(xmin), xmax_(xmax), bins_(1+n) {}
    unsigned numIntervals() const { return bins_.size(); }
    int numEntries(unsigned i) const { return bins_[i]; }
    void fill(double x) {
      if((xmin_ <= x)  && (x < xmax_)) {
        const int i = int(bins_.size() * (x - xmin_)/(xmax_ - xmin_));
        //std::cout<<"fill():  x = "<<x<<", i = "<<i<<" / "<<bins_.size()<<std::endl;
        if(i < bins_.size()) {
          ++bins_[i];
        }
      }
    }
  };
}

//================================================================
void HistAccidentals::init(const std::string& hdir,
                           HistogramFactory &hf,
                           const ConfigFile& conf)
{
  cutPreTrigTimeMax_ = conf.read<double>("MuCapture/Accidentals/tmax");
  const double cycleLength = conf.read<double>("MuCapture/Accidentals/cycleLength");
  const double numCycles = conf.read<double>("MuCapture/Accidentals/numCycles");
  cutPreTrigTimeMin_ = cutPreTrigTimeMax_ - numCycles * cycleLength;

  htstartAll_ = hf.DefineTH1D(hdir, "tstartAll", "t, all windows",
                              8*numCycles, cutPreTrigTimeMin_, cutPreTrigTimeMax_);

  htstartDn_ = hf.DefineTH1D(hdir, "tstartDn", "t, Downstream windows",
                             8*numCycles, cutPreTrigTimeMin_, cutPreTrigTimeMax_);

  hnumwinAll_ = hf.DefineTH1D(hdir, "numWinAll", "number of all windows in range",11, -0.5, 10.5);
  hnumwinDn_ = hf.DefineTH1D(hdir, "numWinDn", "number of downstream windows in range",11, -0.5, 10.5);

  const unsigned maxSubdivisions = conf.read<unsigned>("MuCapture/Accidentals/maxSubdivisions");
  hnumwinAllSubDiv_.resize(1+maxSubdivisions);
  for(unsigned i=0; i<1+maxSubdivisions; ++i) {
    std::ostringstream os;
    os<<"numWinAll_s"<<i;
    hnumwinAllSubDiv_[i] = hf.DefineTH1D(hdir, os.str(), os.str(),11, -0.5, 10.5);
  }

  hnumhitsUp_ = hf.DefineTH2D(hdir, "numHitsUp", "number of PC vs DC hits in upstream windows in range", 90, -0.5, 89.5, 30, -0.5, 29.5);
  hnumhitsUp_->SetOption("colz");
  hnumhitsMixed_ = hf.DefineTH2D(hdir, "numHitsMixed", "number of PC vs DC hits in mixed windows in range", 60, -0.5, 59.5, 30, -0.5, 29.5);
  hnumhitsMixed_->SetOption("colz");
  hnumhitsDn_ = hf.DefineTH2D(hdir, "numHitsDn", "number of PC vs DC hits in downstream windows in range", 60, -0.5, 59.5, 30, -0.5, 29.5);
  hnumhitsDn_->SetOption("colz");

  hOccupancyPCNoDC_.init(hdir, "hitMapPCNoDC", 12, 160, hf, conf);
}

//================================================================
void HistAccidentals::fill(const TimeWindowingResults& wres) {
  int numAll(0), numDn(0);
  std::vector<Subdivision> subdivisions;
  for(unsigned s=0; s < hnumwinAllSubDiv_.size(); ++s) {
    subdivisions.push_back(Subdivision(s, cutPreTrigTimeMin_, cutPreTrigTimeMax_));
  }

  for(unsigned iwin = 0; iwin < wres.iTrigWin; ++iwin) {
    htstartAll_->Fill(wres.windows[iwin].tstart);
    if(wres.windows[iwin].stream == TimeWindow::DOWNSTREAM) {
      htstartDn_->Fill(wres.windows[iwin].tstart);
    }

    if((cutPreTrigTimeMin_ < wres.windows[iwin].tstart ) && (wres.windows[iwin].tstart <= cutPreTrigTimeMax_)) {

      ++numAll;
      if(wres.windows[iwin].stream == TimeWindow::DOWNSTREAM) {
        ++numDn;
      }

      for(unsigned s=0; s<subdivisions.size(); ++s) {
        subdivisions[s].fill(wres.windows[iwin].tstart);
      }

      switch(wres.windows[iwin].stream) {
      case TimeWindow::UPSTREAM:
        hnumhitsUp_->Fill(wres.windows[iwin].dcHits.size(), wres.windows[iwin].pcHits.size());
        break;

      case TimeWindow::MIXED:
        hnumhitsMixed_->Fill(wres.windows[iwin].dcHits.size(), wres.windows[iwin].pcHits.size());
        break;

      case TimeWindow::DOWNSTREAM:
        hnumhitsDn_->Fill(wres.windows[iwin].dcHits.size(), wres.windows[iwin].pcHits.size());
        break;
      } // switch

      if(wres.windows[iwin].dcHits.empty()) {
        hOccupancyPCNoDC_.fill(wres.windows[iwin].pcHits);
      }

    } // if(in range)
  } // for(window)

  hnumwinAll_->Fill(numAll);
  hnumwinDn_->Fill(numDn);

  for(unsigned s=0; s<hnumwinAllSubDiv_.size(); ++s) {
    for(unsigned interval=0; interval < subdivisions[s].numIntervals(); ++interval) {
      hnumwinAllSubDiv_[s]->Fill(subdivisions[s].numEntries(interval));
    }
  }
}

//================================================================
