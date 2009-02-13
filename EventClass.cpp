#include "EventClass.h"

void EventClass::Init( ConfigFile &C, log4cpp::Category *L )
{
	Log	= L;
	//__________ Extract the energy calibration parameters ___________//
	DoEcalib	= false;
	if (C.read<string>("CommandLine/EcalFile") != "")
	{
		if (GetEcalib(C.read<string>("CommandLine/EcalFile"), C.read<string>("CommandLine/EcalArray")))
			DoEcalib	= true;
		else
		{
			Log->error("The energy calibration cannot be performed. Exit now.");
			exit(1);
		}
	}


	MuonWinType	=	C.read<int>( "Parameters/MuonWinType");

	// TODO: code the parsing of the ranges for the window types
	string Tmp;
	Tmp				= C.read<string>( "Parameters/UpDkWinType");
	UpDkWinType		= StrToIntVect(Tmp);
	Tmp				= C.read<string>( "Parameters/DownDkWinType");
	DownDkWinType	= StrToIntVect(Tmp);

	// For Truth Bank
	MuUVMaxRadius		=	C.read<double>( "TruthBank/MuUVMaxRadius");
	TriggerTimeJitter	=	C.read<double>( "TruthBank/TriggerTimeJitter");
	MuZAroundTarget		=	C.read<double>( "TruthBank/MuZAroundTarget");
	MueVertexEpsilon	=	C.read<double>( "TruthBank/MueVertexEpsilon");
	MinDkTimeAfterMu	=	C.read<double>( "TruthBank/MinDkTimeAfterMu");
	MaxDkTimeAfterMu	=	C.read<double>( "TruthBank/MaxDkTimeAfterMu");

	// It is safer to do it at least once
	tb_nmu			= 0;
	tb_mu_trki		= -1;
	tb_mu_firstvtx	= -1;
	tb_mu_lastvtx	= -1;

	tb_e_trki		= -1;
	tb_e_firstvtx	= -1;
	tb_e_lastvtx	= -1;

	cumulative_accflag = 0;
}

bool EventClass::InitGeometry(ConfigFile &C)
{
	BField		=	C.read<double>( "Detector/BField");

	vector<int> TmpVect;
	string Tmp;
	Tmp			=	C.read<string>( "Parameters/DCplanes");
	TmpVect		= StrToIntVect(Tmp);
	for (int i = 0; i < TmpVect.size(); i++)
		GlobToDC[TmpVect[i]]	= i;
	Tmp			=	C.read<string>( "Parameters/PCplanes");
	TmpVect		= StrToIntVect(Tmp);
	for (int i = 0; i < Tmp.size(); i++)
		GlobToPC[Tmp[i]]	= i;

	NumDC	= GlobToDC.size();
	NumPC	= GlobToPC.size();

	Tmp			=	C.read<string>( "Parameters/PlaneType");
	PlaneType	= StrToStrVect(Tmp);

	// Tmp			=	C.read<string>( "Parameters/DCzposition");
	// DCzposition	= StrToFloatVect(Tmp);

	char TmpFile[256];
	if (C.read<int>( "Detector/GeometryFile") != 0)
	{
		int GeoNum = C.read<int>( "Detector/GeometryFile");
		if ( GeoNum < 0)
			Log->warn("Using the default geometry file number %d",abs(GeoNum));
		GeoNum = abs(GeoNum);
		sprintf(TmpFile,"%s/dt_geo.%05d",getenv("CAL_DB"),GeoNum);
		geo = new DetectorGeo();
		Log->info("Reading geometry file %s",TmpFile);
		return geo->ReadGeometry( TmpFile, Log);
	}
	else
	{
		Log->warn("No geometry file given. The geometry of the detector is not set. Do not use any geo->* variable from the EventClass");
		return true;
	}
}

bool EventClass::Load( )
{
	//	#################################### Windowing #################################### //

	imuonwin	= -1;	// Default value. This will crash the prog. That's the goal.
	iewin		= -1;	// Default value. This will crash the prog. That's the goal.
	is_upstreamdk= true;

	if ( not nwin > 0 )
		return false;

	// Look for the muon and the decay windows.
	for (int i = 0; i < nwin; i++)
	{
		if ( win_type[i] == MuonWinType)
		{
			imuonwin	= i;
			// Therefore not a decay.
			continue;
		}
		for(vector<int>::iterator  U=UpDkWinType.begin(); U != UpDkWinType.end(); U++)
			if (win_type[i] == *U)
			{
				iewin			= i;
				is_upstreamdk	= true;
				break;
			}
		for(vector<int>::iterator  D=DownDkWinType.begin(); D != DownDkWinType.end(); D++)
			if (win_type[i] == *D)
			{
				iewin			= i;
				is_upstreamdk	= false;
				break;
			}
	}

	//	#################################### Muon Variables #################################### //

	if (imuonwin != -1)
	{
		if (win_pulast[imuonwin] > win_pvlast[imuonwin])
			muon_plast	= win_pulast[imuonwin];
		else
			muon_plast	= win_pvlast[imuonwin];

        
		muon_ulast	= win_ulast[imuonwin];
		muon_vlast	= win_vlast[imuonwin];
        
		muon_radius = sqrt((muon_ulast * muon_ulast) + (muon_vlast * muon_vlast));
	}
	else
	{
		// default in case no muon window is found.
		muon_plast	= -1;
		muon_ulast	= -1;
		muon_vlast	= -1;
		muon_radius = -1;
	}

	//	####################################### Tracking  ##################################### //

	for( int t = 0; t < ntr; t++)
	{
		// First calculate uncalibrated using the tree variables
		hefit_ptot[t]	= sqrt((hefit_pu[t]*hefit_pu[t])+(hefit_pv[t]*hefit_pv[t])+(hefit_pz[t]*hefit_pz[t]));
		if (not hefit_ptot[t] == 0)
		{
			costh[t]	= hefit_pz[t] / hefit_ptot[t];
			sinth[t]	= sqrt((hefit_pu[t]*hefit_pu[t])+(hefit_pv[t]*hefit_pv[t])) / hefit_ptot[t];
		}
		else
		{
			costh[t]	= 1.0;
			sinth[t]	= 0.0;
		}

		hefit_phi[t]	= atan2(-1*hefit_q[t]*hefit_pv[t] , -1*hefit_q[t]*hefit_pu[t]) - M_PI /2.0;
		hefit_pt[t]	= sqrt((hefit_pu[t]*hefit_pu[t])+(hefit_pv[t]*hefit_pv[t]));

		// Apply the energy calibration if needed
		if (DoEcalib)
		{
			// Upstream
			if (costh[t] < 0.0) 
				ptot[t] = hefit_ptot[t] + ( Ecal_au / fabs(costh[t]) ) - Ecal_bu;

			// Downstream
			if (costh[t] > 0.0)
				ptot[t] = hefit_ptot[t] + ( Ecal_ad / fabs(costh[t]) ) - Ecal_bd;

			pz[t]	= costh[t] * ptot[t];

			pu[t]	= -1 * (sinth[t] * ptot[t] * cos(hefit_phi[t]+M_PI/2.0) );
			pv[t]	= -1 * (sinth[t] * ptot[t] * sin(hefit_phi[t]+M_PI/2.0) );

			// Transverse momentum
			pt[t]		= sqrt((pu[t]*pu[t])+(pv[t]*pv[t]));
			// From that point on, don't use the tree variables even if there is
			// no energy calibration
			// phi[t]		= atan2(-1*hefit_q[t]*pv[t] , -1*hefit_q[t]*pu[t]) - M_PI /2.0
		}
		else
		{
			// By default equal to the tree variables.
			// They will be modified with the energy calibration if needed.
			ptot[t]		= hefit_ptot[t];
			pu[t]		= hefit_pu[t];
			pv[t]		= hefit_pv[t];
			pz[t]		= hefit_pz[t];
			pt[t]		= hefit_pt[t];
			// phi[t]	= hefit_phi[t]
		}
		radius[t]	= 8.4 * ( 2.0 / BField ) * ( 1/50.0 ) * hefit_pt[t];	// From tta, TrackParConverter.h
		////// radius[t]	= ( 0.3 / BField ) * hefit_pt[t];	// From textbooks
		wavelen[t]	= 2.0 * M_PI * radius[t] * (-1*hefit_q[t]*hefit_pz[t] / hefit_pt[t]);
		// cout<<" track "<<t<<" => "<<hefit_pt[t]<<"   "<<radius[t]<<"    "<<wavelen[t]<<endl;
		
		// By default put some impossible values.
		// It seems that some trees from 2004 and before put 0 for the pstart and pstop
		// in the case of ierror=10 (at least). So I have to check that or it crashes.
		// Sometimes the start or stop planes are not DCs.
		// In that case the ierror is not 0 so the track should be thrown away anyway.
		dcmin[t]	= 0;
		dcmax[t]	= 57;
		if ( hefit_pstart[t] != 0 || hefit_pstop[t] != 0 )
		{
			if ( PlaneType[hefit_pstart[t]-1] != "PC" && PlaneType[hefit_pstop[t]-1] != "PC")
			{
				dcmin[t]	= GlobToDC[min(hefit_pstart[t], hefit_pstop[t])] + 1;
				dcmax[t]	= GlobToDC[max(hefit_pstart[t], hefit_pstop[t])] + 1;
			}
		}
	}


	// List of the indices of the selected tracks.
	// To select a track, just add the corresponding index in the list
	// The list of "good" tracks (i.e. ierror=0) is produced here as well.
	seltrack.clear();
	dkwintrack.clear();
	if ( iewin != -1 )
		for (int t = 0; t < ntr; t++)
			if (hefit_winidx[t] == iewin)
			{
				seltrack.push_back(t);
				if (hefit_ierror[t] == 0)		// ZERO IS HARDCODED HERE ... It could probably use the config
					dkwintrack.push_back(t);
			}
	// if ( dkwintrack.size() > 1)
	// 	cout<<" ************** Event nevt = "<<nevt<<"  dkwintrack size = "<<dkwintrack.size()<<endl;

	//	####################################### Truth bank  ##################################### //
	if ( Exists("nmctr") && Exists("nmcvtx") )
	{
		tb_nmu			= 0;
		tb_mu_trki		= -1;
		tb_mu_firstvtx	= -1;
		tb_mu_lastvtx	= -1;

		tb_e_trki		= -1;
		tb_e_firstvtx	= -1;
		tb_e_lastvtx	= -1;

		// First of all, how many muons were simulated in this event ?
		for ( int mctrk = 0; mctrk < nmctr; mctrk++ )
			if(mctrack_pid[mctrk] == 5 || mctrack_pid[mctrk] == 65)
			{
				tb_nmu++;
			}

		for ( int mctrk = 0; mctrk < nmctr; mctrk++ )
		{
			if( ( mctrack_pid[mctrk] == 5 || mctrack_pid[mctrk] == 65 ) && tb_mu_trki < 0)
			{
				tb_mu_firstvtx	= mctrack_voff[mctrk];
				tb_mu_lastvtx	= mctrack_voff[mctrk] + mctrack_nv[mctrk] - 1;
				// If this is the only muon this is easy, just check that it decays.
				if( tb_nmu == 1)
				{
					if( mcvertex_istop[tb_mu_lastvtx] == 1 ) 
						tb_mu_trki = mctrk;
				}
				// More than one muon. We must select only one
				else
				{
					if( mcvertex_istop[tb_mu_lastvtx] == 1 &&							/* Did the muon decay ? */
							fabs(mcvertex_time[tb_mu_lastvtx]) < TriggerTimeJitter &&	/* In time of trigger */
							fabs(mcvertex_vz[tb_mu_lastvtx]) < MuZAroundTarget &&		/* Stop in the target module */
							sqrt(pow(mcvertex_vu[tb_mu_lastvtx],2) + 
							pow(mcvertex_vv[tb_mu_lastvtx],2)) < MuUVMaxRadius)		/* Stop in Beamspot */
						tb_mu_trki = mctrk;
				}
			}

			if( ( mctrack_pid[mctrk] == 2 || mctrack_pid[mctrk] == 3 ) && tb_mu_trki > -1 && tb_e_trki < 0 )
			{
				tb_e_firstvtx	= mctrack_voff[mctrk];
				tb_e_lastvtx	= mctrack_voff[mctrk] + mctrack_nv[mctrk] - 1;
				if( fabs(mcvertex_vu[tb_e_firstvtx] - mcvertex_vu[tb_mu_lastvtx]) < MueVertexEpsilon &&
						fabs(mcvertex_vv[tb_e_firstvtx] - mcvertex_vv[tb_mu_lastvtx]) < MueVertexEpsilon &&
						fabs(mcvertex_vz[tb_e_firstvtx] - mcvertex_vz[tb_mu_lastvtx]) < MueVertexEpsilon &&
						/* positron originates from position of muon stop */
						MinDkTimeAfterMu < mcvertex_time[tb_e_firstvtx] - mcvertex_time[tb_mu_lastvtx] &&
						mcvertex_time[tb_e_firstvtx] - mcvertex_time[tb_mu_lastvtx] < MaxDkTimeAfterMu
						/* Decay time within reconstructable constraints */)
				{
					tb_e_trki	= mctrk;
					for ( int V = tb_e_firstvtx; V < tb_e_lastvtx+1; V++ )
					{
						if ( mcvertex_istop[V] == 10 )
							tb_e_firstdcvtx = V;
						if ( mcvertex_istop[V] == 11 )
							tb_e_lastdcvtx = V;
					}
				}
			}

		}

	}

	return true;
}

void EventClass::InitVar( TTree* T)
{
	//__________________________ EVID branch _________________________//

	GetVar(T, "EVID", "nrun",		&nrun);
	GetVar(T, "EVID", "nevt",		&nevt);


	//__________________________ Event branch ________________________//

	GetVar(T, "Event", "treeversion",	&treeversion);
	GetVar(T, "Event", "timestamp",		&timestamp);
	GetVar(T, "Event", "type",			&type);
	GetVar(T, "Event", "m12width",		&m12width);
	GetVar(T, "Event", "cptime",		cptime);
	GetVar(T, "Event", "rftime",		rftime);
	GetVar(T, "Event", "nwin",			&nwin);
	GetVar(T, "Event", "ntr",			&ntr);
	GetVar(T, "Event", "pienuitr",		&pienuitr);


	//_____________________________ Windows __________________________//

	GetVar(T, "", "win_type",		win_type);
	GetVar(T, "", "win_flag",		win_flag);
	GetVar(T, "", "win_time",		win_time);
	GetVar(T, "", "win_ufirst",		win_ufirst);
	GetVar(T, "", "win_vfirst",		win_vfirst);
	GetVar(T, "", "win_pfirst",		win_pfirst);
	GetVar(T, "", "win_ulast",		win_ulast);
	GetVar(T, "", "win_vlast",		win_vlast);
	GetVar(T, "", "win_pulast",		win_pulast);
	GetVar(T, "", "win_pvlast",		win_pvlast);
	GetVar(T, "", "win_pctof",		win_pctof);
	GetVar(T, "", "win_pcwidthavg",	win_pcwidthavg);


	//__________________________ Helix Fit ___________________________//

	GetVar(T, "", "hefit_winidx",	hefit_winidx);
	GetVar(T, "", "hefit_pu",		hefit_pu);
	GetVar(T, "", "hefit_pv",		hefit_pv);
	GetVar(T, "", "hefit_pz",		hefit_pz);
	GetVar(T, "", "hefit_q",		hefit_q);
	GetVar(T, "", "hefit_u",		hefit_u);
	GetVar(T, "", "hefit_v",		hefit_v);
	GetVar(T, "", "hefit_u0",		hefit_u0);
	GetVar(T, "", "hefit_v0",		hefit_v0);
	GetVar(T, "", "hefit_z",		hefit_z);
	GetVar(T, "", "hefit_pstart",	hefit_pstart);
	GetVar(T, "", "hefit_pstop",	hefit_pstop);
	GetVar(T, "", "hefit_time",		hefit_time);
	GetVar(T, "", "hefit_chi2",		hefit_chi2);
	GetVar(T, "", "hefit_ndof",		hefit_ndof);
	GetVar(T, "", "hefit_ierror",	hefit_ierror);

	GetVar(T, "", "hefit_lastpu",	hefit_lastpu);
	GetVar(T, "", "hefit_lastpv",	hefit_lastpv);
	GetVar(T, "", "hefit_lastpz",	hefit_lastpz);
	GetVar(T, "", "hefit_lastu",	hefit_lastu);
	GetVar(T, "", "hefit_lastv",	hefit_lastv);
	GetVar(T, "", "hefit_lastz",	hefit_lastz);


	//__________________________ CDA ___________________________//

	GetVar(T, "", "dkwin_ncda",		&dkwin_ncda);
	GetVar(T, "", "dkwin_cda",		dkwin_cda);
	GetVar(T, "", "dkwin_cdaz",		dkwin_cdaz);
	GetVar(T, "", "dkwin_cdadefl",	dkwin_cdadefl);


	//_________________________ TEC __________________________//

	GetVar(T, "", "TEC_nhits",			&TEC_nhits);
	GetVar(T, "", "tec_time",			tec_time);
	GetVar(T, "", "tec_width",			tec_width);
	GetVar(T, "", "tec_z",				tec_z);
	GetVar(T, "TEC", "xy0",				xy0);
	GetVar(T, "TEC", "tanthxy",			tanthxy);
	GetVar(T, "TEC", "sigma",			sigma);
	GetVar(T, "TEC", "zspan",			zspan);
	GetVar(T, "TEC", "zrms",			zrms);
	GetVar(T, "TEC", "ifail",			ifail);
	GetVar(T, "", "tec_fitxn",			&tec_fitxn);
	GetVar(T, "", "tec_fityn",			&tec_fityn);
	GetVar(T, "", "tec_fitxt",			tec_fitxt);
	GetVar(T, "", "tec_fitxz",			tec_fitxz);
	GetVar(T, "", "tec_xres",			tec_xres);
	GetVar(T, "", "tec_fityt",			tec_fityt);
	GetVar(T, "", "tec_fityz",			tec_fityz);
	GetVar(T, "", "tec_yres",			tec_yres);
	GetVar(T, "", "nhits_in_scin",		&nhits_in_scin);
	GetVar(T, "", "failed_cptof_m12",	&failed_cptof_m12);
	GetVar(T, "", "global_wire",		global_wire);
	GetVar(T, "", "TC_xy_TWIST",		TC_xy_TWIST);


	//_________________________ FirstLastDC  __________________________//

	GetVar(T, "", "win_pudcfirst",		win_pudcfirst);
	GetVar(T, "", "win_pudclast",		win_pudclast);
	GetVar(T, "", "win_pvdcfirst",		win_pvdcfirst);
	GetVar(T, "", "win_pvdclast",		win_pvdclast);
	GetVar(T, "", "win_udcfirst",		win_udcfirst);
	GetVar(T, "", "win_udclast",		win_udclast);
	GetVar(T, "", "win_vdcfirst",		win_vdcfirst);
	GetVar(T, "", "win_vdclast",		win_vdclast);


	//_________________________ FGOutPut  __________________________//

	GetVar(T, "", "nfgtr",			&nfgtr);
	GetVar(T, "", "fgfit_winidx",	fgfit_winidx);
	GetVar(T, "", "fgfit_pu",		fgfit_pu);
	GetVar(T, "", "fgfit_pv",		fgfit_pv);
	GetVar(T, "", "fgfit_pz",		fgfit_pz);
	GetVar(T, "", "fgfit_u",		fgfit_u);
	GetVar(T, "", "fgfit_v",		fgfit_v);
	GetVar(T, "", "fgfit_z",		fgfit_z);
	GetVar(T, "", "fgfit_pstart",	fgfit_pstart);
	GetVar(T, "", "fgfit_pstop",	fgfit_pstop);
	GetVar(T, "", "fgfit_time",		fgfit_time);
	GetVar(T, "", "fgfit_chi2",		fgfit_chi2);
	GetVar(T, "", "fgfit_ndof",		fgfit_ndof);
	GetVar(T, "", "fgfit_ierror",	fgfit_ierror);
	GetVar(T, "", "fgfit_type",		fgfit_type);


	//_________________________ MHitsOutput  __________________________//

	GetVar(T, "", "nmhits",		&nmhits);
	GetVar(T, "", "time",		time);


	//_________________________ PactOutput  __________________________//

	GetVar(T, "", "pact_elost",		pact_elost);


	//_________________________ PseudoPact  __________________________//

	GetVar(T, "", "npc5",				&npc5);
	GetVar(T, "", "npc6",				&npc6);
	GetVar(T, "", "pspact_pc5wire",		pc5wire);
	GetVar(T, "", "pspact_pc5time",		pc5time);
	GetVar(T, "", "pspact_pc5width",	pc5width);
	GetVar(T, "", "pspact_pc6wire",		pc6wire);
	GetVar(T, "", "pspact_pc6time",		pc6time);
	GetVar(T, "", "pspact_pc6width",	pc6width);


	//_________________________ PseudoPactTest  __________________________//

	GetVar(T, "", "npc7",					&npc7);
	GetVar(T, "", "npc8",					&npc8);
	GetVar(T, "", "pspacttest_pc7wire",		pc7wire);
	GetVar(T, "", "pspacttest_pc7time",		pc7time);
	GetVar(T, "", "pspacttest_pc7width",	pc7width);
	GetVar(T, "", "pspacttest_pc8wire",		pc8wire);
	GetVar(T, "", "pspacttest_pc8time",		pc8time);
	GetVar(T, "", "pspacttest_pc8width",	pc8width);


	//_________________________ WinStatOutput  __________________________//

	GetVar(T, "", "win_numdc",			win_numdc);
	GetVar(T, "", "win_numpc",			win_numpc);
	GetVar(T, "", "win_maxuv",			win_maxuv);
	GetVar(T, "", "win_uavg",			win_uavg);
	GetVar(T, "", "win_vavg",			win_vavg);
	GetVar(T, "", "win_usig",			win_usig);
	GetVar(T, "", "win_vsig",			win_vsig);
	GetVar(T, "", "win_hitspp",			win_hitspp);
	GetVar(T, "", "win_clareaavg",		win_clareaavg);
	GetVar(T, "", "win_dcwidthavg",		win_dcwidthavg);
	GetVar(T, "", "win_pctsig",			win_pctsig);
	GetVar(T, "", "win_up_dense_start",	win_up_dense_start);
	GetVar(T, "", "win_dn_dense_start",	win_dn_dense_start);


	//_________________________ HeFitNHitsOutput  __________________________//

	GetVar(T, "", "hefit_numu",			hefit_numu);
	GetVar(T, "", "hefit_numv",			hefit_numv);
	GetVar(T, "", "hefit_nunused",		hefit_nunused);
	GetVar(T, "", "hefit_nmissingpl",	hefit_nmissingpl);


	//_________________________ FgFitNHitsOutput  __________________________//

	GetVar(T, "", "fgfit_numu",			fgfit_numu);
	GetVar(T, "", "fgfit_numv",			fgfit_numv);
	GetVar(T, "", "fgfit_nunused",		fgfit_nunused);
	GetVar(T, "", "fgfit_nmissingpl",	fgfit_nmissingpl);


	//_________________________ MCBankOutput  __________________________//

	GetVar(T, "MCSet", "spectrum",		&spectrum);
	GetVar(T, "MCSet", "sample",		&sample);
	GetVar(T, "MCSet", "ndecays",		&ndecays);
	GetVar(T, "MCSet", "spin3",			&spin3);
	GetVar(T, "", "nmctr",				&nmctr);
	GetVar(T, "", "mctrack_itrack",		mctrack_itrack);
	GetVar(T, "", "mctrack_pid",		mctrack_pid);
	GetVar(T, "", "mctrack_voff",		mctrack_voff);
	GetVar(T, "", "mctrack_nv",			mctrack_nv);
	GetVar(T, "", "nmcvtx",				&nmcvtx);
	GetVar(T, "", "mcvertex_ptot",		mcvertex_ptot);
	GetVar(T, "", "mcvertex_costh",		mcvertex_costh);
	GetVar(T, "", "mcvertex_phimuv",	mcvertex_phimuv);
	GetVar(T, "", "mcvertex_time",		mcvertex_time);
	GetVar(T, "", "mcvertex_vu",		mcvertex_vu);
	GetVar(T, "", "mcvertex_vv",		mcvertex_vv);
	GetVar(T, "", "mcvertex_vz",		mcvertex_vz);
	GetVar(T, "", "mcvertex_istop",		mcvertex_istop);
	GetVar(T, "", "micheld_itrack",		micheld_itrack);
	GetVar(T, "", "micheld_accflag",	micheld_accflag);


	//_________________________ User  __________________________//

	GetVar(T, "", "event_user",		event_user);
	GetVar(T, "", "win_user1",		win_user1);
	GetVar(T, "", "win_user2",		win_user2);
	GetVar(T, "", "win_user3",		win_user3);
	GetVar(T, "", "win_user4",		win_user4);
	GetVar(T, "", "track_user1",	track_user1);
	GetVar(T, "", "track_user2",	track_user2);
	GetVar(T, "", "track_user3",	track_user3);
	GetVar(T, "", "track_user4",	track_user4);
	GetVar(T, "", "mctrack_user1",	mctrack_user1);
	GetVar(T, "", "mctrack_user2",	mctrack_user2);
	GetVar(T, "", "mctrack_user3",	mctrack_user3);
	GetVar(T, "", "mctrack_user4",	mctrack_user4);
}


void EventClass::GetVar( TTree* T, const char* Branch, const char* Leaf, Int_t* V)
{
	if ( Branch == "")
		if (CheckBranchLeaf( T, Leaf))
		{
			ExistList[Leaf]	= true;
			T->GetLeaf(Leaf)->SetAddress(V);
		}
		else
			ExistList[Leaf]	= false;
	else
		if (CheckBranchLeaf( T, Branch, Leaf))
		{
			ExistList[Leaf]	= true;
			T->GetBranch(Branch)->GetLeaf(Leaf)->SetAddress(V);
		}
		else
			ExistList[Leaf]	= false;
}

void EventClass::GetVar( TTree* T, const char* Branch, const char* Leaf, Float_t* V)
{
	if ( Branch == "")
		if (CheckBranchLeaf( T, Leaf))
		{
			ExistList[Leaf]	= true;
			T->GetLeaf(Leaf)->SetAddress(V);
		}
		else
			ExistList[Leaf]	= false;
	else
		if (CheckBranchLeaf( T, Branch, Leaf))
		{
			ExistList[Leaf]	= true;
			T->GetBranch(Branch)->GetLeaf(Leaf)->SetAddress(V);
		}
		else
			ExistList[Leaf]	= false;
}

bool EventClass::CheckBranchLeaf( TTree* T, const char* Branch, const char* Leaf)
{
	TBranch *tmpBranch;
	if (Branch != "")
		if ( NULL == T->GetBranch(Branch))
		{
			return false;
		}
		else
		{
			tmpBranch	= T->GetBranch(Branch);
		}

	
	if ( NULL == tmpBranch->GetLeaf(Leaf))
	{
		return false;
	}

	return true;
}

bool EventClass::CheckBranchLeaf( TTree* T, const char* Leaf)
{
	if ( NULL == T->GetLeaf(Leaf))
	{
		return false;
	}

	return true;
}


bool EventClass::GetEcalib(string EcalibFile, string EcalibArray)
{
	TFile *File	= new TFile(EcalibFile.c_str());
	if (File->IsZombie())
	{
		Log->error("The energy calibration file %s is not opened correctly.",EcalibFile.c_str());
		return false;
	}

	Ecal_au	= -999.9;
	Ecal_ad	= -999.9;
	Ecal_bu	= -999.9;
	Ecal_bd	= -999.9;

	TObjArray *Array	= (TObjArray*)File->Get(EcalibArray.c_str());
	if (Array == NULL)
	{
		Log->error("The array %s was not found in the energy calibration file %s.", EcalibArray.c_str(), EcalibFile.c_str());
		return false;
	}
	TObjArray *SubArr		= (TObjArray*) Array->At(0);
	TVectorT<double> *Vect	= (TVectorT<double>*) Array->At(1);
	// vector<string> ParName;
	// vector<double> ParVal;
	TObjString *Tmp;
	string Name;
	double Val;

	// Let's convert the ugly TObjArray/TObjString/TVector organisation into something more convenient
	for ( int i = 0; i < SubArr->GetEntries(); i++)
	{
		// ParName.push_back(SubArr->At(i)->GetString()->Data());
		// ParVal.push_back(Vect->GetMatrixArray()[i]);
		Tmp		= (TObjString*)SubArr->At(i);
		Name	= Tmp->GetString().Data();
		Val		= Vect->GetMatrixArray()[i];

		if (Name == "au" or Name == "upstslope")
			Ecal_au	= Val;
		if (Name == "ad" or Name == "dnstslope")
			Ecal_ad	= Val;
		if (Name == "beta")
		{
			Ecal_bu	= Val;
			Ecal_bd	= Val;
		}
		if (Name == "upstintercept")
			Ecal_bu	= Val;
		if (Name == "dnstintercept")
			Ecal_bd	= Val;
	}


	File->Close();
	if ( -999.9 == Ecal_au || -999.9 == Ecal_bu || -999.9 == Ecal_ad || -999.9 == Ecal_bd)
	{
		Log->error("The array %s does not contain the good parameters.", EcalibArray.c_str());
		return false;
	}

	// TODO : Output the details of the energy calibration
	Log->info("The energy calibration will be applied.");

	return true;
}

