//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EventClass_h
#define EventClass_h
// Include de C++
using namespace std;

#include <iostream>
#include <string>
#include <map>
#include <math.h>

// Include de ROOT
// #include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
// #include "TMath.h"
#include "TObjArray.h"
#include "TVectorT.h"
#include "TObjString.h"

// Include Program
#include "log4cpp/Category.hh"
#include "ConfigFile.h"
#include "FuncLib.h"
#include "DetectorGeo.h"

#define MAXWIN			21
#define MAXTECHITS		384
#define MAXHIT			24		// TEC max hit per track (one hit per wire in each module)
#define MAXPLANE		2		// TEC number of modules ... yes 2
#define MAXTRK			150		// Max number of helix fit tracks
#define MAXMCTRK		1000	// Max number of MC tracks
#define MAX_WIRES_P		160		// Max number of hits in PsPACT planes (?)
#define MAXMICHELDDKS	32		// Max number of Michel decay from micheld server.
#define MAXCDA			190		// Max number of CDA assuming a max of 20 good tracks.

class EventClass{
	public :
		EventClass()	{};
		~EventClass()	{};

		void Init( ConfigFile &C, log4cpp::Category *L );
		bool InitGeometry(ConfigFile &C);
		void InitVar( TTree* T);
		void GetVar( TTree* T, const char* Branch, const char* Leaf, Int_t* V);
		void GetVar( TTree* T, const char* Branch, const char* Leaf, Float_t* V);
		bool Load( );
		bool GetEcalib(string EcalibFile, string EcalibArray);

		bool Exists( const char* Name)	{return ExistList[Name];}

		log4cpp::Category *Log;
		// Geometry class of the TWIST detector
		DetectorGeo *geo;

		//__________________________ Config _________________________//
		bool DoEcalib;
		double BField;
		int NumDC;
		int NumPC;
		
		//__________________________ Parameters _________________________//
		float	pi;

		int MuonWinType;
		vector<int> UpDkWinType;
		vector<int> DownDkWinType;

		map<int, int> GlobToDC;
		map<int, int> GlobToPC;
		vector<string> PlaneType;
		// vector<float> DCzposition;
		
		// For the truth bank
		double MuUVMaxRadius;
		double TriggerTimeJitter;
		double MuZAroundTarget;
		double MueVertexEpsilon;
		double MinDkTimeAfterMu;
		double MaxDkTimeAfterMu;

		
		//________________________ Energy calibration ___________________//

		double	Ecal_au;		// upstream slope from energy calibration
		double	Ecal_ad;		// downstream slope from energy calibration
		double	Ecal_bu;		// upstream intercept from energy calibration
		double	Ecal_bd;		// downstream intercept from energy calibration

		//__________________________ EVID branch _________________________//

		Int_t		nrun;
		Int_t		nevt;


		//__________________________ Event branch ________________________//

		Int_t		treeversion;
		Int_t		timestamp;
		Int_t		type;
		Float_t		m12width;
		Float_t		cptime[3];
		Float_t		rftime[3];
		Int_t		nwin;
		Int_t		ntr;
		Int_t		pienuitr;


		//_____________________________ Windows __________________________//

		Int_t		win_type[MAXWIN];			//[nwin]
		Int_t		win_flag[MAXWIN];			//[nwin]
		Float_t		win_time[MAXWIN];			//[nwin]
		Float_t		win_ufirst[MAXWIN];			//[nwin]
		Float_t		win_vfirst[MAXWIN];			//[nwin]
		Int_t		win_pfirst[MAXWIN];			//[nwin]
		Float_t		win_ulast[MAXWIN];			//[nwin]
		Float_t		win_vlast[MAXWIN];			//[nwin]
		Int_t		win_pulast[MAXWIN];			//[nwin]
		Int_t		win_pvlast[MAXWIN];			//[nwin]
		Float_t		win_pctof[MAXWIN];			//[nwin]
		Float_t		win_pcwidthavg[MAXWIN];		//[nwin]


		//__________________________ Helix Fit ___________________________//

		Int_t		hefit_winidx[MAXTRK]; 		//[ntr]
		Float_t		hefit_pu[MAXTRK];     		//[ntr]
		Float_t		hefit_pv[MAXTRK];     		//[ntr]
		Float_t		hefit_pz[MAXTRK];     		//[ntr]
		Int_t		hefit_q[MAXTRK];      		//[ntr]
		Float_t		hefit_u[MAXTRK];      		//[ntr]
		Float_t		hefit_v[MAXTRK];      		//[ntr]
		Float_t		hefit_u0[MAXTRK];     		//[ntr]
		Float_t		hefit_v0[MAXTRK];     		//[ntr]
		Float_t		hefit_z[MAXTRK];      		//[ntr]
		Int_t		hefit_pstart[MAXTRK]; 		//[ntr]
		Int_t		hefit_pstop[MAXTRK];  		//[ntr]
		Float_t		hefit_time[MAXTRK];   		//[ntr]
		Float_t		hefit_chi2[MAXTRK];   		//[ntr]
		Int_t		hefit_ndof[MAXTRK];   		//[ntr]
		Int_t		hefit_ierror[MAXTRK]; 		//[ntr]


		//__________________________ CDA ___________________________//

		Int_t		dkwin_ncda;
		Float_t		dkwin_cda[MAXCDA];      		//[ncda]
		Float_t		dkwin_cdaz[MAXCDA];      		//[ncda]
		Float_t		dkwin_cdadefl[MAXCDA];     		//[ncda]


		//_________________________ TEC __________________________//

		Float_t		xy0[MAXPLANE];
		Float_t		tanthxy[MAXPLANE];
		Float_t		sigma[MAXPLANE];
		Float_t		zspan[MAXPLANE];
		Float_t		zrms[MAXPLANE];
		Int_t		ifail[MAXPLANE];
		Int_t		tec_fitxn;
		Int_t		tec_fityn;
		Float_t		tec_fitxt[MAXHIT];			//[tec_fitxn]
		Float_t		tec_fitxz[MAXHIT];			//[tec_fitxn]
		Float_t		tec_xres[MAXHIT]; 			//[tec_fitxn]
		Float_t		tec_fityt[MAXHIT];			//[tec_fityn]
		Float_t		tec_fityz[MAXHIT];			//[tec_fityn]
		Float_t		tec_yres[MAXHIT]; 			//[tec_fityn]
		Int_t		TEC_nhits;
		Float_t		tec_time[MAXTECHITS];		//[TEC_nhits]
		Float_t		tec_width[MAXTECHITS];		//[TEC_nhits]
		Float_t		tec_z[MAXTECHITS];   		//[TEC_nhits]
		Int_t		nhits_in_scin;
		Int_t		failed_cptof_m12;
		Int_t		global_wire[MAXTECHITS];	//[TEC_nhits]
		Float_t		TC_xy_TWIST[MAXTECHITS];	//[TEC_nhits]


		//_________________________ FirstLastDC  __________________________//

		Int_t		win_pudcfirst[MAXWIN];		//[nwin]
		Int_t		win_pudclast[MAXWIN];		//[nwin]
		Int_t		win_pvdcfirst[MAXWIN];		//[nwin]
		Int_t		win_pvdclast[MAXWIN];		//[nwin]
		Float_t		win_udcfirst[MAXWIN];		//[nwin]
		Float_t		win_udclast[MAXWIN];		//[nwin]
		Float_t		win_vdcfirst[MAXWIN];		//[nwin]
		Float_t		win_vdclast[MAXWIN];		//[nwin]


		//_________________________ FGOutPut  __________________________//

		Int_t		nfgtr;
		Int_t		fgfit_winidx[MAXTRK];		//[nfgtr]
		Float_t		fgfit_pu[MAXTRK];    		//[nfgtr]
		Float_t		fgfit_pv[MAXTRK];    		//[nfgtr]
		Float_t		fgfit_pz[MAXTRK];    		//[nfgtr]
		Float_t		fgfit_u[MAXTRK];     		//[nfgtr]
		Float_t		fgfit_v[MAXTRK];     		//[nfgtr]
		Float_t		fgfit_z[MAXTRK];     		//[nfgtr]
		Int_t		fgfit_pstart[MAXTRK];		//[nfgtr]
		Int_t		fgfit_pstop[MAXTRK]; 		//[nfgtr]
		Float_t		fgfit_time[MAXTRK];  		//[nfgtr]
		Float_t		fgfit_chi2[MAXTRK];  		//[nfgtr]
		Int_t		fgfit_ndof[MAXTRK];  		//[nfgtr]
		Int_t		fgfit_ierror[MAXTRK];		//[nfgtr]
		Int_t		fgfit_type[MAXTRK];  		//[nfgtr]


		//_________________________ MHitsOutput  __________________________//

		Int_t		nmhits;
		Float_t		time[15];


		//_________________________ PactOutput  __________________________//

		Float_t		pact_elost[4];


		//_________________________ PseudoPact  __________________________//

		Int_t		npc5;
		Int_t		npc6;
		Int_t		pc5wire[MAX_WIRES_P];		//[npc5]
		Float_t		pc5time[MAX_WIRES_P];		//[npc5]
		Float_t		pc5width[MAX_WIRES_P];		//[npc5]
		Int_t		pc6wire[MAX_WIRES_P];		//[npc6]
		Float_t		pc6time[MAX_WIRES_P];		//[npc6]
		Float_t		pc6width[MAX_WIRES_P];		//[npc6]


		//_________________________ PseudoPactTest  __________________________//

		Int_t		npc7;
		Int_t		npc8;
		Int_t		pc7wire[MAX_WIRES_P];  	//[npc7]
		Float_t		pc7time[MAX_WIRES_P];  	//[npc7]
		Float_t		pc7width[MAX_WIRES_P]; 	//[npc7]
		Int_t		pc8wire[MAX_WIRES_P];  	//[npc8]
		Float_t		pc8time[MAX_WIRES_P];  	//[npc8]
		Float_t		pc8width[MAX_WIRES_P]; 	//[npc8]


		//_________________________ WinStatOutput  __________________________//

		Int_t		win_numdc[MAXWIN];          		//[nwin]
		Int_t		win_numpc[MAXWIN];          		//[nwin]
		Float_t		win_maxuv[MAXWIN];          		//[nwin]
		Float_t		win_uavg[MAXWIN];           		//[nwin]
		Float_t		win_vavg[MAXWIN];           		//[nwin]
		Float_t		win_usig[MAXWIN];           		//[nwin]
		Float_t		win_vsig[MAXWIN];           		//[nwin]
		Float_t		win_hitspp[MAXWIN];         		//[nwin]
		Float_t		win_clareaavg[MAXWIN];      		//[nwin]
		Float_t		win_dcwidthavg[MAXWIN];     		//[nwin]
		Float_t		win_pctsig[MAXWIN];         		//[nwin]
		Int_t		win_up_dense_start[MAXWIN]; 		//[nwin]
		Int_t		win_dn_dense_start[MAXWIN]; 		//[nwin]


		//_________________________ HeFitNHitsOutput  __________________________//

		Int_t		hefit_numu[MAXTRK];					//[ntr]
		Int_t		hefit_numv[MAXTRK];					//[ntr]
		Int_t		hefit_nunused[MAXTRK];				//[ntr]
		Int_t		hefit_nmissingpl[MAXTRK];			//[ntr]


		//_________________________ FgFitNHitsOutput  __________________________//

		Int_t		fgfit_numu[MAXTRK];					//[nfgtr]
		Int_t		fgfit_numv[MAXTRK];					//[nfgtr]
		Int_t		fgfit_nunused[MAXTRK];				//[nfgtr]
		Int_t		fgfit_nmissingpl[MAXTRK];			//[nfgtr]


		//_________________________ MCBankOutput  __________________________//

		Int_t		spectrum;
		Int_t		sample;
		Int_t		ndecays;
		Float_t		spin3;
		Int_t		nmctr;
		Int_t		mctrack_itrack[MAXMCTRK];			//[nmctr]
		Int_t		mctrack_pid[MAXMCTRK];   			//[nmctr]
		Int_t		mctrack_voff[MAXMCTRK];  			//[nmctr]
		Int_t		mctrack_nv[MAXMCTRK];    			//[nmctr]
		Int_t		nmcvtx;
		Float_t		mcvertex_ptot[MAXMCTRK];			//[nmcvtx]
		Float_t		mcvertex_costh[MAXMCTRK];			//[nmcvtx]
		Float_t		mcvertex_phimuv[MAXMCTRK];			//[nmcvtx]
		Float_t		mcvertex_time[MAXMCTRK];			//[nmcvtx]
		Float_t		mcvertex_vu[MAXMCTRK];				//[nmcvtx]
		Float_t		mcvertex_vv[MAXMCTRK];				//[nmcvtx]
		Float_t		mcvertex_vz[MAXMCTRK];				//[nmcvtx]
		Int_t		mcvertex_istop[MAXMCTRK];			//[nmcvtx]
		Int_t		micheld_itrack[MAXMICHELDDKS];  	//[ndecays]
		Int_t		micheld_accflag[MAXMICHELDDKS]; 	//[ndecays]


		//_________________________ MCBankOutput  __________________________//

		Float_t		event_user[10];
		Float_t		win_user1[MAXWIN]; 					//[nwin]
		Float_t		win_user2[MAXWIN]; 					//[nwin]
		Float_t		win_user3[MAXWIN]; 					//[nwin]
		Float_t		win_user4[MAXWIN]; 					//[nwin]
		Float_t		track_user1[MAXTRK];				//[ntr]
		Float_t		track_user2[MAXTRK];				//[ntr]
		Float_t		track_user3[MAXTRK];				//[ntr]
		Float_t		track_user4[MAXTRK];				//[ntr]
		Float_t		mctrack_user1[MAXMCTRK];			//[nmctr]
		Float_t		mctrack_user2[MAXMCTRK];			//[nmctr]
		Float_t		mctrack_user3[MAXMCTRK];			//[nmctr]
		Float_t		mctrack_user4[MAXMCTRK];			//[nmctr]


		//______________ New variables (not tree variables) ______________//

		Float_t		ptot[MAXTRK];		// Total momentum (calibrated if energy calibration)
		Float_t		hefit_ptot[MAXTRK];	// Total momentum (NOT calibrated)
		Float_t		pu[MAXTRK];			// Track pu (calibrated if energy calibration)
		Float_t		pv[MAXTRK];			// Track pv (calibrated if energy calibration)
		Float_t		pz[MAXTRK];			// Track pz (calibrated if energy calibration)
		Float_t		costh[MAXTRK];		// Track cosinus(theta)
		Float_t		sinth[MAXTRK];		// Track sinus(theta)
		Float_t		hefit_phi[MAXTRK];	// Track phi angle (NOT calibrated)
		Float_t		hefit_pt[MAXTRK];	// Transverse momentum (NOT calibrated)
		Float_t		pt[MAXTRK];			// Transverse momentum (calibrated if energy calibration)
		Float_t		radius[MAXTRK];		// Radius of the helix
		Float_t		wavelen[MAXTRK];	// Wavelength of the helix
		Float_t		dcmin[MAXTRK];		// DC min of the track
		Float_t		dcmax[MAXTRK];		// DC max of the track


		int			imuonwin;			// Muon window index
		int			iewin;				// Decay window index
		bool		is_upstreamdk;		// Direction of decay from classification
		int			muon_plast;			// Muon last plane
		float		muon_ulast;			// Muon last U position
		float		muon_vlast;			// Muon last V position
		float		muon_radius;		// Muon last position radius

		vector<int>	seltrack;			// Selected track at a given point in the analysis
		vector<int>	dkwintrack;			// All good tracks (ierror==0) in the decay window

		
		// Truth bank
		int			tb_nmu;				// Number of muons in the event
		int			tb_mu_trki;			// Muon track index
		int			tb_mu_firstvtx;		// Muon first vertex index
		int			tb_mu_lastvtx;		// Muon last vertex index

		int			tb_e_trki;			// Positron track index
		int			tb_e_firstvtx;		// Positron first vertex index
		int			tb_e_lastvtx;		// Positron last vertex index
		int			tb_e_firstdcvtx;	// Positron vertex index at the first DC
		int			tb_e_lastdcvtx;		// Positron vertex index at the last DC



	private :
		map <const char*, bool>ExistList;

		bool CheckBranchLeaf( TTree* T, const char* Branch, const char* Leaf);
		bool CheckBranchLeaf( TTree* T, const char* Leaf);

};


#endif
