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

class DetectorGeo {
	public :
		DetectorGeo()	{};
		~DetectorGeo()	{};

		bool ReadGeometry( string geofile);

		// (1) DC variables
		int		geo_ndfoils;                 // Total number of foils in drift chambers
		int		geo_ndplanes;                // Total number of drift planes
		int		geo_ndwires_physical;        // Number of physical wires/DC plane
		int		geo_ndwires[MAX_PLANES_D];   // Number of instrumented wires/DC plane
		double	geo_zdfoil[MAX_FOILS_D];   // Array of z coords for drift foils 
		double	geo_zdplane[MAX_PLANES_D]; // Array of z coords for drift planes
		double	geo_dshift[MAX_PLANES_D];  // Array of nominal drift plane shifts
		double	geo_drot[MAX_PLANES_D];    // Array of nominal drift plane rotations
		double	geo_rd_myl;                // Radius of DC nylar foil
		double	geo_td_myl;                // Thickness of DC mylar foil
		double	geo_dw_rad;                // Radius of DC sense wires
		double	geo_dw_len;                // Length of DC sense wires
		double	geo_dw_space;              // Spacing between wires in DC
		double	geo_td_plane;              // Thickness of drift plane

		// (2) PC variables    
		int		geo_npfoils;                  // Total number of foils in prop. chambers
		int		geo_npplanes;                 // Total number of prop. planes
		int		geo_npwires_physical;         // Number of physical wires/PC plane
		int		geo_npwires[MAX_PLANES_P];    // Number of instrumented wires/PC plane
		double	geo_zpfoil[MAX_FOILS_P];    // Array of z coords for prop. foils 
		double	geo_zpplane[MAX_PLANES_P];  // Array of z coords for prop. planes
		double	geo_pshift[MAX_PLANES_P];   // Array of nominal prop. plane shifts
		double	geo_prot[MAX_PLANES_P];     // Array of nominal prop. plane rotations
		double	geo_rp_myl,geo_tp_myl;          // Radius, thickness of mylar foil in PC
		double	geo_pw_rad,geo_pw_len;          // Radius, length of PC sense wires
		double	geo_pw_space;               // Spacing between wires in PC
		double	geo_tp_plane;               // Thickness of proportional plane
  
		// (3) Target variables

		// TAR2
		double	geo_targ_rad;           // Radius of stopping target
		double	geo_targ_thick;         // Thickness of stopping target
		int		geo_itarg_mat;          // Medium of target(=material)
		double	geo_cond_thick;         // Thickness of target conduction layer
		// The following doesn't appear anywhere
		// int		geo_icond_med;          // Medium of target conduction layer

		// Additional variables for TAR3
		double geo_tfoil_r1, geo_tfoil_r2, geo_tfoil_thick;
		double geo_cfoil_r1, geo_cfoil_r2, geo_cfoil_thick;

		// (4) Scint variables

		int		geo_nmuscints;                      // Total number of beam scintillators
		int		geo_imuplane[MAXM_SCINTS];
		int		geo_nmuwires[MAXM_SCINTS];
		double	geo_zmuscint[MAXM_SCINTS];       // z-coords of mu scintillators
		double	geo_rmuscint[MAXM_SCINTS];       // Radius of mu scintillators
		double	geo_tmuscint[MAXM_SCINTS];       // Thickness of mu scintillators

		int		geo_nt0scints;                      // Total number of t0 scintillators
		int		geo_it0plane[MAXT_SCINTS];
		int		geo_nt0wires[MAXT_SCINTS];
		double	geo_zt0scint[MAXT_SCINTS];       // z-coords of t0 scintillators
		double	geo_rt0scint[2][MAXT_SCINTS];    // Inner/Outer Radius of t0 scintillators
		double	geo_tt0scint[MAXT_SCINTS];       // Thickness of t0 scintillators

	private :
		bool	ReadDRFT( ifstream &file);
		bool	ReadPROP( ifstream &file);
		bool	ReadTARX( ifstream &file);
		bool	ReadSCIX( ifstream &file);
};

#endif
