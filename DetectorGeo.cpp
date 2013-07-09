#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>

#include <boost/regex.hpp>

#include "log4cpp/Category.hh"

#include "DetectorGeo.h"
#include "ConfigFile.h"

//================================================================
WirePlane::WirePlane(WirePlane::DetType d,
                     double z,
                     double rotDegrees,
                     double alignmentShift)
  : det_(d)
  , wireSpacing_(d == PC ? 0.2 : 0.4)
  , centralWire_(d == PC ? 161./2. : 81./2.)
{
  static const double pi = 4.*atan(1.);
  const double rot = rotDegrees * pi/180.;

  // This is equivalent to mvector computation in chambers_mod.f90
  // but we just subtract 45 degrees instead of doing xyz2uvz3
  mvector_ = ROOT::Math::XYZVector(cos(rot-pi/4.),
                                   sin(rot-pi/4.),
                                   0.);

  dir_ = (std::abs(mvector_.x()) > std::abs(mvector_.y())) ? U : V;

  center_ = ROOT::Math::XYZPoint((dir_ == U)? alignmentShift : 0.,
                                 (dir_ == V)? alignmentShift : 0.,
                                 z);
}

//================================================================
std::string WirePlane::detName(WirePlane::DetType d) {
  switch(d) {
  case PC: return "PC";
  case DC: return "DC";
  default: throw std::runtime_error("WirePlane::detName(DetType d): invalid input");
  }
}

//================================================================
WirePlane::Measurement WirePlane::measurement(double wireNumber) const {
  ROOT::Math::XYZPoint wireCenter = center_ +
    mvector_ * wireSpacing_ * (wireNumber - centralWire_);

  return Measurement(dir_,
                     (dir_ == U) ? wireCenter.x() : wireCenter.y()
                     );
}

//================================================================
ROOT::Math::XYPoint WirePlane::uv(Measurement m1, Measurement m2) {
  if(m1.dir == m2.dir) {
    throw std::runtime_error("ERROR: WirePlane::uv(): inputs measure the same direction.");
  }
  const double u = (m1.dir == U) ? m1.coordinate : m2.coordinate;
  const double v = (m1.dir == V) ? m1.coordinate : m2.coordinate;
  return ROOT::Math::XYPoint(u, v);
}

//================================================================
namespace {
  struct WirePlaneZSorter {
    bool operator()(const WirePlane& a, const WirePlane& b) const {
      return a.center().z() < b.center().z();
    }
  };
}

DetectorGeo::DetectorGeo(const ConfigFile& conf, log4cpp::Category& logger) {
  if(!ReadGeometry(conf, &logger)) {
    throw std::runtime_error("DetectorGeo: error reading geometry.");
  }

  for(int i=0; i<npplanes; ++i) {
    pcplanes_.push_back(WirePlane(WirePlane::PC, zpplane[i], prot[i], pshift[i]));
  }

  for(int i=0; i<ndplanes; ++i) {
    dcplanes_.push_back(WirePlane(WirePlane::DC, zdplane[i], drot[i], dshift[i]));
  }

  std::copy(pcplanes_.begin(), pcplanes_.end(), std::back_inserter(globalplanes_));
  std::copy(dcplanes_.begin(), dcplanes_.end(), std::back_inserter(globalplanes_));
  std::sort(globalplanes_.begin(), globalplanes_.end(), WirePlaneZSorter());
}

//================================================================
const std::vector<WirePlane>& DetectorGeo::planes(WirePlane::DetType d) const {
  switch(d) {
  case WirePlane::PC: return pcplanes_;
  case WirePlane::DC: return dcplanes_;
  default: throw std::runtime_error("DetectorGeo::planes(DetType d): invalid input");
  }
}

//================================================================
unsigned DetectorGeo::maxCellNumber(WirePlane::DetType d) const {
  switch(d) {
  case WirePlane::PC: return 160;
  case WirePlane::DC: return 80;
  default: throw std::runtime_error("DetectorGeo::maxCellNumber(DetType d): invalid input");
  }
}

//================================================================
bool DetectorGeo::ReadGeometry(const ConfigFile& conf, log4cpp::Category *L)
{
  int GeoNum = conf.read<int>("Detector/GeometryFile", 0);
  if(!GeoNum) {
    L->warn("No geometry file given. The geometry of the detector is not set. Do not use any geo->* variable from the EventClass");
    return false;
  }

  if(GeoNum < 0) {
    L->warn("Using the default geometry file number %d",abs(GeoNum));
    GeoNum = std::abs(GeoNum);
  }

  char TmpFile[256];
  sprintf(TmpFile,"%s/dt_geo.%05d",getenv("CAL_DB"),GeoNum);
  L->info("Reading geometry file %s",TmpFile);
  const string geofile(TmpFile);

	std::ifstream file(geofile.c_str());
	if(!file)
	{
		L->error("Failed to open geometry file %s", geofile.c_str());
		return false;
	}

	char  ReadIn[280];

	do 
	{
		file.getline(ReadIn,280,'\n');

		// (1) Drift Chamber

		if(!strncmp(ReadIn,"DRFT", 4)) 
		{
			ReadDRFT(file);
		}
		
		// (2) PC variables

		if(!strncmp(ReadIn,"PROP", 4)) 
		{
			ReadPROP(file);
		}

		// (3) Target variables

		if( ( !strncmp(ReadIn,"TAR2", 4))		// "old" target
			or ( !strncmp(ReadIn,"TAR3", 4)))	// "new" target
		{
                  ReadTARX(file, ReadIn);
		}

		// (4) Scint variables

		if( ( !strncmp(ReadIn,"SCI2", 4))
			or ( !strncmp(ReadIn,"SCI3", 4)))
		{
                  ReadSCIX(file, ReadIn);
		}

	} while (! file.eof() );   // while not eof
	file.close();

        applyCorrections(conf, *L);

	return true;
}

// =============== (1) Drift Chamber ================ //
void DetectorGeo::ReadDRFT(std::ifstream &file)
{
	char  ReadIn[280];
	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');    

	//Radius, thickness of mylar foil in DC
	file >> rd_myl >> td_myl ; 
	// std::cout << "rd_myl, td_myl\n";
	// std::cout << rd_myl << ", " << td_myl << "\n";

	//Total number of foils in drift chamber
	file >>ndfoils;
	// std::cout<<"ndfoils " <<  ndfoils <<std::endl;
	// std::cout << "foils\n";
	for(int idfoil=0;idfoil<ndfoils;idfoil++) 
	{
		//Array of z coords for drift foils 
		file >>zdfoil[idfoil]; 
		// std::cout << idfoil << " "<< zdfoil[idfoil] << "\n";
	}

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// Radius, length of DC sense wires ,Spacing between wires in DC and,
	// Number of physical wires/DC plane

	file >>dw_rad>>dw_len>>dw_space>> ndwires_physical;
	// std::cout <<" dw_rad, dw_len, dw_space, ndwires_physical \n";
	// std::cout <<dw_rad<<", " <<dw_len<<", " <<dw_space<<", " 
	//	<< ndwires_physical << "\n";

	//Total number of drift planes and Thickness of drift plane

	file >>ndplanes>>td_plane;
	// std::cout<<"ndplanes " << ndplanes <<std::endl;
	// std::cout<<"td_plane " << td_plane <<std::endl;

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// drift plane positions
	int plane;
	for(int iplane=0;iplane<ndplanes;iplane++) 
	{
		file>> plane >>zdplane[iplane]>>dshift[iplane]
			>>drot[iplane]>>ndwires[iplane];
		// std::cout<<iplane<< ", " <<zdplane[iplane]<< ", " 
		//     <<dshift[iplane]<< ", " <<drot[iplane]
		//     << ", " <<ndwires[iplane]<< "\n";
	} 
	// std::cout << "rd_myl, td_myl\n";
	// std::cout << rd_myl << ", " << td_myl << "\n";
}


// ================== (2) PC variables ==================== //
void DetectorGeo::ReadPROP(std::ifstream &file)
{
	char  ReadIn[280];
	// skip to and over the + sign
	do
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// Radius, thickness of mylar foil in PC
	file>>rp_myl>>tp_myl;
	// std::cout << "Mylar radius, thickness " <<  rp_myl << " "<<tp_myl <<"\n";

	//Total number of foils in prop. chamber
	file >> npfoils;
	// std::cout<<"npfoils " << npfoils <<std::endl;

	//position the Mylar foils defining PC wire planes

	for(int ipfoil=0; ipfoil<npfoils; ipfoil++)
	{
		file >> zpfoil[ipfoil];
		// std::cout << ipfoil << " " << zpfoil[ipfoil] << "\n";
	}

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// Radius, length of PC sense wires, spacing between wires in PC and,
	// Number of physical wires/PC plane

	file >> pw_rad >> pw_len >> pw_space >> npwires_physical;
	// std::cout << "pw_rad, pw_len, pw_space, npwires_physical \n";
	// std::cout << pw_rad << " " << pw_len << " " << pw_space << " " 
	//	<< npwires_physical << "\n";

	//Total number of prop planes and Thickness of gas volume

	file >>npplanes>>tp_plane;
	// std::cout<<"npplanes, tp_plane " << npplanes << tp_plane<<std::endl;

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// plane positions in PC
	int plane;
	for(int ipplane=0; ipplane<npplanes; ipplane++) 
	{
		file >> plane >> zpplane[ipplane] >> pshift[ipplane] >> prot[ipplane] 
			>> npwires[ipplane];
		// std::cout << ipplane << " " <<zpplane[ipplane]<< " " <<pshift[ipplane]
		//	    << " " <<prot[ipplane]<< " " <<npwires[ipplane]<< "\n";
	} 
}


// ================ (3) Target variables ================== //
void DetectorGeo::ReadTARX(std::ifstream &file, char* ReadIn)
{
	// Read in target info
	if(!strncmp(ReadIn,"TAR2", 4)) // "old" target
	{
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');

		// Radius, thickness and medium of target foil
		// thickness of conduction layer
		// TARGET MEDIUM IS HARD CODED IN GEANT 4!! Maher. Feb 6, 2003.
		file>> targ_rad >> targ_thick >> itarg_mat >> cond_thick; 
		// std::cout<< targ_rad << targ_thick << itarg_mat << cond_thick << "\n"; 
	}  //end readin of Targ

	else if(!strncmp(ReadIn,"TAR3", 4)) // "new" target
	{
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');
		file >> tfoil_r1 >> tfoil_r2 >> tfoil_thick >> 
			targ_rad >> targ_thick >> itarg_mat >>
			cfoil_r1 >> cfoil_r2 >> cfoil_thick;
		// std::cout << tfoil_r1 << tfoil_r2 << tfoil_thick << 
		//  targ_rad << targ_thick << itarg_mat <<
		//  cfoil_r1 << cfoil_r2 << cfoil_thick;
	}
}


// ================ (4) Scint variables =================== //
void DetectorGeo::ReadSCIX(std::ifstream &file, char *ReadIn)
{
	// Read in scintillator info
	if(!strncmp(ReadIn,"SCI2", 4)) 
	{
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');

		//Total number of mu scints, t0 scints
		file>>nmuscints >> nt0scints;
		/*
		   std::cout<<"nmuscints " << nmuscints <<std::endl;
		   std::cout<<"nt0scints " << nt0scints <<std::endl;
		   */
		for(int isc=0;isc<nmuscints;isc++)
		{
			file>> rmuscint[isc] >> tmuscint[isc];
			file>> zmuscint[isc];
			rmuscint[isc] = rmuscint[isc]; 
			tmuscint[isc]=tmuscint[isc];
			zmuscint[isc]=zmuscint[isc];
			/*
			   std::cout << rmuscint[isc] << " " 
			   <<tmuscint[isc] << " " 
			   << zmuscint[isc] << std::endl;
			   */
		}

		//Inner/outer radius and thickness of t0 scintillator
		for(int isc=0;isc<nt0scints;isc++)
		{
			file>> rt0scint[0][isc] >> rt0scint[1][isc] >> tt0scint[isc];
			file>> zt0scint[isc];
			rt0scint[0][isc]=rt0scint[0][isc];
			rt0scint[1][isc]=rt0scint[1][isc];
			tt0scint[isc]=tt0scint[isc];
			zt0scint[isc]=zt0scint[isc];
			/*
			   std::cout << rt0scint[0][isc] << " "
			   << rt0scint[1][isc] << " " << tt0scint[isc] << " " 
			   << zt0scint[isc] <<std::endl;
			   */
		}
	}  
	else if(!strncmp(ReadIn,"SCI3", 4)) 
	{
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');

		//Total number of mu scints, t0 scints
		file >> nmuscints;
		//	 std::cout<<"nmuscints " << nmuscints <<std::endl;

		for(int isc=0;isc<nmuscints;isc++)
		{
		file >> imuplane[isc] >> nmuwires[isc] >> 
				rmuscint[isc] >> tmuscint[isc] >> zmuscint[isc];
			rmuscint[isc]=rmuscint[isc]; 
			tmuscint[isc]=tmuscint[isc]; 
			zmuscint[isc]=zmuscint[isc];
			/*
			   std::cout << imuplane[isc] << " " << nmuwires[isc] << " "
			   << rmuscint[isc] << " " << tmuscint[isc] << " " 
			   << zmuscint[isc] << std::endl;
			   */
		}

		file >> nt0scints;
		//	 std::cout<<"nt0scints " << nt0scints <<std::endl;
		//Inner/outer radius and thickness of t0 scintillator
		for(int isc=0;isc<nt0scints;isc++)
		{
			file >> it0plane[isc] >> nt0wires[isc] >> rt0scint[0][isc] >>
				rt0scint[1][isc] >> tt0scint[isc] >> zt0scint[isc];
			rt0scint[0][isc]=rt0scint[0][isc];
			rt0scint[1][isc]=rt0scint[1][isc]; 
			tt0scint[isc]=tt0scint[isc]; 
			zt0scint[isc]=zt0scint[isc];
			/*
			   std::cout << it0plane[isc] << " " << nt0wires[isc] << " "
			   << rt0scint[0][isc] << " " << rt0scint[1][isc] << " "
			   << tt0scint[isc] << " " << zt0scint[isc] << std::endl;
			   */
		}
	}  //end readin of SCIN
}

//================================================================
void DetectorGeo::applyCorrections(const ConfigFile& conf, log4cpp::Category& logger) {
  const int ipc_ppc = conf.read<int>("Detector/Geometry/pc_ppc", 0);
  if(ipc_ppc) {
    std::ostringstream os;
    os<<getenv("CAL_DB")<<"/pc_ppc."<<std::setfill('0')<<std::setw(5)<<ipc_ppc;
    logger.infoStream()<<"Reading PC UV corrections file "<<os.str();
    std::ifstream file(os.str().c_str());
    if(!file.is_open()) {
      throw std::runtime_error("Could not open calibratio file "+os.str());
    }

    applyCorrectionsPPC(os.str(), logger, pshift, MAX_PLANES_P);
  }
  else {
    logger.warn("NO PC UV corrections file Detector/Geometry/pc_ppc is specified");
  }

  //----------------
  const int idc_ppc = conf.read<int>("Detector/Geometry/dc_ppc", 0);
  if(idc_ppc) {
    std::ostringstream os;
    os<<getenv("CAL_DB")<<"/dc_ppc."<<std::setfill('0')<<std::setw(5)<<idc_ppc;
    logger.infoStream()<<"Reading DC UV corrections file "<<os.str();
    std::ifstream file(os.str().c_str());
    if(!file.is_open()) {
      throw std::runtime_error("Could not open calibratio file "+os.str());
    }

    applyCorrectionsPPC(os.str(), logger, dshift, MAX_PLANES_D);
  }
  else {
    logger.warn("NO DC UV corrections file Detector/Geometry/dc_ppc is specified");
  }
}

//================================================================
void DetectorGeo::applyCorrectionsPPC(const std::string& calibFileName, log4cpp::Category& logger, double *posArray, int arrsize) {
  std::ifstream file(calibFileName.c_str());
  if(!file.is_open()) {
    throw std::runtime_error("Could not open calibration file "+calibFileName);
  }

  std::string line;
  const boost::regex COMMENT("^C");
  while(std::getline(file, line)) {
    if(!boost::regex_search(line, COMMENT)) {
      // std::cout<<"data line: "<<line<<std::endl;
      unsigned iplane(0);
      double correction(0.);
      std::istringstream is(line);
      if(! (is>>iplane>>correction) ) {
        throw std::runtime_error("DetectorGeo::applyCorrections(): error parsing the line \""+line+"\"");
      }

      // DetectorGeo uses 0-based plane numbering, but rest of the code is 1-based.
      if( (iplane < 1) || (iplane > arrsize)) {
        throw std::runtime_error("DetectorGeo::applyCorrectionsPPC(): iplane is out of range");
      }
      posArray[iplane-1] += correction;
    }
  }
}

//================================================================
