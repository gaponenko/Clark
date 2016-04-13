#ifndef RooUnfoldDummy_h
#define RooUnfoldDummy_h

class TH2D;

class RooUnfoldDummy {
public:
  void Setup(int, double, double) {}
  void Setup(int, double, double, int, double, double) {}
  void Setup(TH2D*, TH2D*) {}

  void Fill(double, double) {}
  void Miss(double) {}

  void Fill(int, double, int, double) {}
  void Miss(int, double) {}

};

#endif/*RooUnfoldDummy_h*/
