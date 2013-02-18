//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef DetectorGeo_h
#define DetectorGeo_h
// Include de C++
using namespace std;

// DCs
#define MAX_FOILS_D 60
#define MAX_PLANES_D 44
#define MAX_WIRES_D 80
// PCs
#define MAX_FOILS_P 15
#define MAX_PLANES_P 12
#define MAX_WIRES_P 160
// Scintillator
#define MAXM_SCINTS 2
#define MAXT_SCINTS 4
#define MAXO_SCINTS 128

#include <iostream>
#include <fstream>
#include <string>

// Include Program
#include "FuncLib.h"
#include "log4cpp/Category.hh"

class DetectorGeo {
	public :
		DetectorGeo()	{};
		~DetectorGeo()	{};

		bool ReadGeometry( string geofile, log4cpp::Category *L);

		// (1) DC variables
		int		ndfoils;                 // Total number of foils in drift chambers
		int		ndplanes;                // Total number of drift planes
		int		ndwires_physical;        // Number of physical wires/DC plane
		int		ndwires[MAX_PLANES_D];   // Number of instrumented wires/DC plane
		double	zdfoil[MAX_FOILS_D];   // Array of z coords for drift foils 
		double	zdplane[MAX_PLANES_D]; // Array of z coords for drift planes
		double	dshift[MAX_PLANES_D];  // Array of nominal drift plane shifts
		double	drot[MAX_PLANES_D];    // Array of nominal drift plane rotations
		double	rd_myl;                // Radius of DC nylar foil
		double	td_myl;                // Thickness of DC mylar foil
		double	dw_rad;                // Radius of DC sense wires
		double	dw_len;                // Length of DC sense wires
		double	dw_space;              // Spacing between wires in DC
		double	td_plane;              // Thickness of drift plane

		// (2) PC variables    
		int		npfoils;                  // Total number of foils in prop. chambers
		int		npplanes;                 // Total number of prop. planes
		int		npwires_physical;         // Number of physical wires/PC plane
		int		npwires[MAX_PLANES_P];    // Number of instrumented wires/PC plane
		double	zpfoil[MAX_FOILS_P];    // Array of z coords for prop. foils 
		double	zpplane[MAX_PLANES_P];  // Array of z coords for prop. planes
		double	pshift[MAX_PLANES_P];   // Array of nominal prop. plane shifts
		double	prot[MAX_PLANES_P];     // Array of nominal prop. plane rotations
		double	rp_myl,tp_myl;          // Radius, thickness of mylar foil in PC
		double	pw_rad,pw_len;          // Radius, length of PC sense wires
		double	pw_space;               // Spacing between wires in PC
		double	tp_plane;               // Thickness of proportional plane
  
		// (3) Target variables

		// TAR2
		double	targ_rad;           // Radius of stopping target
		double	targ_thick;         // Thickness of stopping target
		int		itarg_mat;          // Medium of target(=material)
		double	cond_thick;         // Thickness of target conduction layer
		// The following doesn't appear anywhere
		// int		icond_med;          // Medium of target conduction layer

		// Additional variables for TAR3
		double tfoil_r1, tfoil_r2, tfoil_thick;
		double cfoil_r1, cfoil_r2, cfoil_thick;

		// (4) Scint variables

		int		nmuscints;                      // Total number of beam scintillators
		int		imuplane[MAXM_SCINTS];
		int		nmuwires[MAXM_SCINTS];
		double	zmuscint[MAXM_SCINTS];       // z-coords of mu scintillators
		double	rmuscint[MAXM_SCINTS];       // Radius of mu scintillators
		double	tmuscint[MAXM_SCINTS];       // Thickness of mu scintillators

		int		nt0scints;                      // Total number of t0 scintillators
		int		it0plane[MAXT_SCINTS];
		int		nt0wires[MAXT_SCINTS];
		double	zt0scint[MAXT_SCINTS];       // z-coords of t0 scintillators
		double	rt0scint[2][MAXT_SCINTS];    // Inner/Outer Radius of t0 scintillators
		double	tt0scint[MAXT_SCINTS];       // Thickness of t0 scintillators

	private :
		void	ReadDRFT( ifstream &file);
		void	ReadPROP( ifstream &file);
		void	ReadTARX( ifstream &file);
		void	ReadSCIX( ifstream &file);

};

#endif
