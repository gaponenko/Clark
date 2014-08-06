// Andrei Gaponenko, 2013

#ifndef TDCPIDCalib_h
#define TDCPIDCalib_h

#include <string>
#include <vector>

class ConfigFile;

struct TDCPIDCalibParameters {
  double x0;
  double y0;
  double slope;
  TDCPIDCalibParameters() : x0(), y0(), slope() {}
};

//================================================================
class TDCPIDCalib {
public:
  void init(const std::string& configPrefix, const ConfigFile &conf);

  std::pair<bool,double> yCalibrated(int multiplicityBin, double x, double yraw) const;

private :
  typedef std::vector<TDCPIDCalibParameters> CalibMap;
  CalibMap db_;
  std::string configPrefix_; // keep this only to identify ourself in runtime exceptions
};

#endif/*TDCPIDCalib_h*/
