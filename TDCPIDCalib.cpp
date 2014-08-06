#include "TDCPIDCalib.h"

#include <limits>
#include <sstream>
#include <stdexcept>

#include "ConfigFile.h"

void TDCPIDCalib::init(const std::string& configPrefix, const ConfigFile &conf) {
  configPrefix_ = configPrefix;

  db_.resize(1); // dummy entry at[0], we start with 1-sized clusters

  for(int mult=1; ; ++mult) {
    std::ostringstream os;
    os<<configPrefix<<"/"<<mult;
    if(!conf.keyExists(os.str())) {
      break;
    }

    std::cout<<"Initializing calib: "<<os.str()<<std::endl;
    std::istringstream is(conf.read<std::string>(os.str()));
    TDCPIDCalibParameters par;
    if(!(is>>par.x0)) throw std::runtime_error("Error reading x0 for "+os.str());
    if(!(is>>par.y0)) throw std::runtime_error("Error reading y0 for "+os.str());
    if(!(is>>par.slope)) throw std::runtime_error("Error reading slope for "+os.str());
    std::string dummy;
    if(is>>dummy) throw std::runtime_error("Error: trailing garbage for "+os.str()+": "+dummy);

    db_.push_back(par);
  }
}

std::pair<bool,double> TDCPIDCalib::yCalibrated(int multiplicityBin, double x, double yraw) const {
  std::pair<bool,double> res(false, std::numeric_limits<double>::quiet_NaN());
  if(multiplicityBin < db_.size()) {
    const TDCPIDCalibParameters& par = db_[multiplicityBin];
    res.first = true;
    res.second = yraw - (par.y0 + par.slope*(x - par.x0));
  }
  return res;
}
