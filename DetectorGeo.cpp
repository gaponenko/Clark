#include "DetectorGeo.h"

bool DetectorGeo::ReadGeometry( string geofile)
{
	
	std::ifstream file(geofile.c_str());
	if(!file)
	{
		std::cerr << "Failed to open geometry file\n";
		//     return(EXIT_FAILURE);
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
			ReadTARX(file);
		}

		// (4) Scint variables

		if( ( !strncmp(ReadIn,"SCI2", 4))
			or ( !strncmp(ReadIn,"SCI3", 4)))
		{
			ReadSCIX(file);
		}

	} while (! file.eof() );   // while not eof
	file.close();
}

// =============== (1) Drift Chamber ================ //
bool DetectorGeo::ReadDRFT( ifstream &file)
{
	char  ReadIn[280];
	std::cout<<" Found DRFT"<<std::endl;
	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');    

	//Radius, thickness of mylar foil in DC
	file >> geo_rd_myl >> geo_td_myl ; 
	// std::cout << "rd_myl, td_myl\n";
	// std::cout << rd_myl << ", " << td_myl << "\n";

	//Total number of foils in drift chamber
	file >>geo_ndfoils;
	// std::cout<<"ndfoils " <<  ndfoils <<std::endl;
	// std::cout << "foils\n";
	for(int idfoil=0;idfoil<geo_ndfoils;idfoil++) 
	{
		//Array of z coords for drift foils 
		file >>geo_zdfoil[idfoil]; 
		// std::cout << idfoil << " "<< zdfoil[idfoil] << "\n";
	}

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// Radius, length of DC sense wires ,Spacing between wires in DC and,
	// Number of physical wires/DC plane

	file >>geo_dw_rad>>geo_dw_len>>geo_dw_space>> geo_ndwires_physical;
	// std::cout <<" dw_rad, dw_len, dw_space, ndwires_physical \n";
	// std::cout <<dw_rad<<", " <<dw_len<<", " <<dw_space<<", " 
	//	<< ndwires_physical << "\n";

	//Total number of drift planes and Thickness of drift plane

	file >>geo_ndplanes>>geo_td_plane;
	// std::cout<<"geo_ndplanes " << geo_ndplanes <<std::endl;
	// std::cout<<"td_plane " << td_plane <<std::endl;

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// drift plane positions
	int plane;
	for(int iplane=0;iplane<geo_ndplanes;iplane++) 
	{
		file>> plane >>geo_zdplane[iplane]>>geo_dshift[iplane]
			>>geo_drot[iplane]>>geo_ndwires[iplane];
		// std::cout<<iplane<< ", " <<zdplane[iplane]<< ", " 
		//	   <<dshift[iplane]<< ", " <<drot[iplane]
		//	   << ", " <<ndwires[iplane]<< "\n";
	} 
	// std::cout << "rd_myl, td_myl\n";
	// std::cout << rd_myl << ", " << td_myl << "\n";
}


// ================== (2) PC variables ==================== //
bool DetectorGeo::ReadPROP( ifstream &file)
{
	char  ReadIn[280];
	// skip to and over the + sign
	do
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// Radius, thickness of mylar foil in PC
	file>>geo_rp_myl>>geo_tp_myl;
	// std::cout << "Mylar radius, thickness " <<  geo_rp_myl << " "<<geo_tp_myl <<"\n";

	//Total number of foils in prop. chamber
	file >> geo_npfoils;
	// std::cout<<"geo_npfoils " << geo_npfoils <<std::endl;

	//position the Mylar foils defining PC wire planes

	for(int ipfoil=0; ipfoil<geo_npfoils; ipfoil++)
	{
		file >> geo_zpfoil[ipfoil];
		// std::cout << ipfoil << " " << geo_zpfoil[ipfoil] << "\n";
	}

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// Radius, length of PC sense wires, spacing between wires in PC and,
	// Number of physical wires/PC plane

	file >> geo_pw_rad >> geo_pw_len >> geo_pw_space >> geo_npwires_physical;
	// std::cout << "geo_pw_rad, geo_pw_len, geo_pw_space, geo_npwires_physical \n";
	// std::cout << geo_pw_rad << " " << geo_pw_len << " " << geo_pw_space << " " 
	//	<< geo_npwires_physical << "\n";

	//Total number of prop planes and Thickness of gas volume

	file >>geo_npplanes>>geo_tp_plane;
	// std::cout<<"geo_npplanes, geo_tp_plane " << geo_npplanes << geo_tp_plane<<std::endl;

	// skip to and over the + sign
	do 
	{
		file.getline(ReadIn,280,'\n');
	} while(ReadIn[0] != '+');

	// plane positions in PC
	int plane;
	for(int ipplane=0; ipplane<geo_npplanes; ipplane++) 
	{
		file >> plane >> geo_zpplane[ipplane] >> geo_pshift[ipplane] >> geo_prot[ipplane] 
			>> geo_npwires[ipplane];
		// std::cout << ipplane << " " <<geo_zpplane[ipplane]<< " " <<geo_pshift[ipplane]
		//	    << " " <<geo_prot[ipplane]<< " " <<geo_npwires[ipplane]<< "\n";
	} 
}


// ================ (3) Target variables ================== //
bool DetectorGeo::ReadTARX( ifstream &file)
{
	char  ReadIn[280];
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
		file>> geo_targ_rad >> geo_targ_thick >> geo_itarg_mat >> geo_cond_thick; 
		// std::cout<< geo_targ_rad << geo_targ_thick << geo_itarg_mat << geo_cond_thick << "\n"; 
	}  //end readin of Targ

	else if(!strncmp(ReadIn,"TAR3", 4)) // "new" target
	{
		std::cout<<" Found TAR3"<<std::endl;
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');
		file >> geo_tfoil_r1 >> geo_tfoil_r2 >> geo_tfoil_thick >> 
			geo_targ_rad >> geo_targ_thick >> geo_itarg_mat >>
			geo_cfoil_r1 >> geo_cfoil_r2 >> geo_cfoil_thick;
		// std::cout << geo_tfoil_r1 << geo_tfoil_r2 << geo_tfoil_thick << 
		//  geo_targ_rad << geo_targ_thick << geo_itarg_mat <<
		//  geo_cfoil_r1 << geo_cfoil_r2 << geo_cfoil_thick;
	}
}


// ================ (4) Scint variables =================== //
bool DetectorGeo::ReadSCIX( ifstream &file)
{
	char  ReadIn[280];
	// Read in scintillator info
	if(!strncmp(ReadIn,"SCI2", 4)) 
	{
		std::cout<<" Found SCI2"<<std::endl;
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');

		//Total number of mu scints, t0 scints
		file>>geo_nmuscints >> geo_nt0scints;
		/*
		   std::cout<<"geo_nmuscints " << geo_nmuscints <<std::endl;
		   std::cout<<"geo_nt0scints " << geo_nt0scints <<std::endl;
		   */
		for(int isc=0;isc<geo_nmuscints;isc++)
		{
			file>> geo_rmuscint[isc] >> geo_tmuscint[isc];
			file>> geo_zmuscint[isc];
			geo_rmuscint[isc] = geo_rmuscint[isc]; 
			geo_tmuscint[isc]=geo_tmuscint[isc];
			geo_zmuscint[isc]=geo_zmuscint[isc];
			/*
			   std::cout << geo_rmuscint[isc] << " " 
			   <<geo_tmuscint[isc] << " " 
			   << geo_zmuscint[isc] << std::endl;
			   */
		}

		//Inner/outer radius and thickness of t0 scintillator
		for(int isc=0;isc<geo_nt0scints;isc++)
		{
			file>> geo_rt0scint[0][isc] >> geo_rt0scint[1][isc] >> geo_tt0scint[isc];
			file>> geo_zt0scint[isc];
			geo_rt0scint[0][isc]=geo_rt0scint[0][isc];
			geo_rt0scint[1][isc]=geo_rt0scint[1][isc];
			geo_tt0scint[isc]=geo_tt0scint[isc];
			geo_zt0scint[isc]=geo_zt0scint[isc];
			/*
			   std::cout << geo_rt0scint[0][isc] << " "
			   << geo_rt0scint[1][isc] << " " << geo_tt0scint[isc] << " " 
			   << geo_zt0scint[isc] <<std::endl;
			   */
		}
	}  
	else if(!strncmp(ReadIn,"SCI3", 4)) 
	{
		std::cout<<" Found SCI3"<<std::endl;
		// skip to and over the + sign
		do 
		{
			file.getline(ReadIn,280,'\n');
		} while(ReadIn[0] != '+');

		//Total number of mu scints, t0 scints
		file >> geo_nmuscints;
		//	 std::cout<<"geo_nmuscints " << geo_nmuscints <<std::endl;

		for(int isc=0;isc<geo_nmuscints;isc++)
		{
			file >> geo_imuplane[isc] >> geo_nmuwires[isc] >> 
				geo_rmuscint[isc] >> geo_tmuscint[isc] >> geo_zmuscint[isc];
			geo_rmuscint[isc]=geo_rmuscint[isc]; 
			geo_tmuscint[isc]=geo_tmuscint[isc]; 
			geo_zmuscint[isc]=geo_zmuscint[isc];
			/*
			   std::cout << geo_imuplane[isc] << " " << geo_nmuwires[isc] << " "
			   << geo_rmuscint[isc] << " " << geo_tmuscint[isc] << " " 
			   << geo_zmuscint[isc] << std::endl;
			   */
		}

		file >> geo_nt0scints;
		//	 std::cout<<"geo_nt0scints " << geo_nt0scints <<std::endl;
		//Inner/outer radius and thickness of t0 scintillator
		for(int isc=0;isc<geo_nt0scints;isc++)
		{
			file >> geo_it0plane[isc] >> geo_nt0wires[isc] >> geo_rt0scint[0][isc] >>
				geo_rt0scint[1][isc] >> geo_tt0scint[isc] >> geo_zt0scint[isc];
			geo_rt0scint[0][isc]=geo_rt0scint[0][isc];
			geo_rt0scint[1][isc]=geo_rt0scint[1][isc]; 
			geo_tt0scint[isc]=geo_tt0scint[isc]; 
			geo_zt0scint[isc]=geo_zt0scint[isc];
			/*
			   std::cout << geo_it0plane[isc] << " " << geo_nt0wires[isc] << " "
			   << geo_rt0scint[0][isc] << " " << geo_rt0scint[1][isc] << " "
			   << geo_tt0scint[isc] << " " << geo_zt0scint[isc] << std::endl;
			   */
		}
	}  //end readin of SCIN
}


