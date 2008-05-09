//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 

class StatusHistograms : public ModuleClass{
	public :
		StatusHistograms()		{};
		StatusHistograms(string N, string T)		{Name = N; Title = T;};
		~StatusHistograms()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string Name;
		string Title;

		bool PerWindowType;
};

bool StatusHistograms::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Status Histograms");

	if (not Conf.read<bool>("StatusHistograms/Do", 1))
	{
		Log->info( "StatusHistograms %s turned OFF",Title.c_str());
		return false ;
	}

	string N = Name;
	string T = Title;

	//	 --------- Parameters initialization ---------		//
	PerWindowType	= Conf.read<bool>("StatusHistograms/PerWindowType");
	
	//	 --------- Histograms initialization ---------		//
	//______________________ Event Branch __________________________ //
	H.DefineTH1D( "Status",	"EvtType_"+N,	"Event Type "+T+";Event type",32, -0.5,31.5);
	//// These two can be optained by projection of the 2D plot
	// H.DefineTH1D( "Status",	"TCAP_"+N,		"TCAP "+T+";TCAP [ns]",200, 0,200);
	// H.DefineTH1D( "Status",	"m12width_"+N,	"m12 width "+T+";m12 width [ns]",300, 0, 30000);

	H.DefineTH2D( "Status",	"m12VsTCAP_"+N,	"m12 vs TCAP "+T+";TCAP [ns];m12 width [ns]",200, 0,200,300, 0, 30000);
	H.DefineTH1D( "Status",	"NumWin_"+N,	"Number of windows "+T,11, -0.5,10.5);
	H.DefineTH1D( "Status",	"NumTrk_"+N,	"Number of tracks "+T,11, -0.5,10.5);

	//_________________________ Windows ________________________ //
	H.DefineTH1D( "Status",	"Win_Type_"+N,	"Window Type "+T,22, -0.5,21.5);
	H.DefineTH2D( "Status",	"Win_FlagVsType_"+N, "Windowing flag versus window type "+T,22, -0.5,21.5, 31, -0.5,30.5);
	H.DefineTH2D( "Status", "Win_PCWidthAvg_"+N, "Average PC width vs window type "+T+";Window type;Average PC width [ns]",22,-0.5,21.5,400,0,400);
	H.DefineTH2D( "Status", "Win_PCtof_"+N, "PC time of flight vs window type "+T+";Window type;PC time of flight [ns]",22,-0.5,21.5,300,-150,150);
	if (PerWindowType)
	{
		string wtstr;
		for (int W = 0; W < 22; W++)
		{
			wtstr	= IntToStr(W);
			H.DefineTH1D( "Status",	"Win_Time_PerWinType_"+wtstr+"_"+N,	"Window time for the window type "+wtstr+" "+T+";Window time [ns]",320, -6000,10000);
			H.DefineTH2D( "Status",	"Win_uvfirst_PerWinType_"+wtstr+"_"+N,	"First plane v vs u position for the window type "+wtstr+" "+T+";u position [cm];v position [cm]",80, -16.0,16.0,80, -16.0,16.0);
			H.DefineTH2D( "Status",	"Win_uvlast_PerWinType_"+wtstr+"_"+N,	"Last plane v vs u position for the window type "+wtstr+" "+T+";u position [cm];v position [cm]",80, -16.0,16.0,80, -16.0,16.0);
			H.DefineTH1D( "Status",	"Win_pfirst_PerWinType_"+wtstr+"_"+N,	"First plane for the window type "+wtstr+" "+T, 56, 0.5,56.5);
			H.DefineTH2D( "Status",	"Win_pupvlast_PerWinType_"+wtstr+"_"+N,	"Last planes v vs last plane u for the window type "+wtstr+" "+T,56, 0.5,56.5, 56, 0.5,56.5);
		}
	}

	//___________________ Helix fit tracks ______________________ //
	H.DefineTH1D( "Status",	"Chi2DkWinTrk_"+N,		"Reduced Chisquare of the tracks in the decay window "+T,1000, 0,10);
	H.DefineTH1D( "Status",	"Chi2DkWinTrkFull_"+N,	"Reduced Chisquare of the tracks in the decay window "+T,1000, 0,1000);
	
	H.DefineTH1D( "Status",	"Chi2SelTrk_"+N,		"Reduced Chisquare of the selected tracks "+T,1000, 0,10);
	H.DefineTH1D( "Status",	"Chi2SelTrkFull_"+N,	"Reduced Chisquare of the selected tracks "+T,1000, 0,1000);
	
	H.DefineTH1D( "Status",	"NdofDkWinTrk_"+N,		"Number of degrees of freedom of the tracks in the decay window "+T,1000, 0,10);
	H.DefineTH1D( "Status",	"NdofDkWinTrkFull_"+N,	"Number of degrees of freedom of the tracks in the decay window "+T,1000, 0,1000);
	
	H.DefineTH1D( "Status",	"NdofSelTrk_"+N,		"Number of degrees of freedom of the selected tracks "+T,1000, 0,10);
	H.DefineTH1D( "Status",	"NdofSelTrkFull_"+N,	"Number of degrees of freedom of the selected tracks "+T,1000, 0,1000);




	return true;
}

bool StatusHistograms::Process(EventClass &E, HistogramFactory &H)
{
	// ____ ____ //
		//_________ Event Branch
		H.Fill("EvtType_"+Name,E.type);
		H.Fill("NumWin_"+Name,E.nwin);
		H.Fill("NumTrk_"+Name,E.ntr);
		// H.Fill("m12width_"+Name,E.m12width);

		for ( int i = 0; i < 3; i++)
			// H.Fill("TCAP_"+Name,E.cptime[i]);
			H.Fill("m12VsTCAP_"+Name,E.cptime[i], E.m12width);

		//_________ Windows
		string wtstr;
		for ( int w = 0; w < E.nwin; w++)
		{
			wtstr	= E.win_type[w];
			H.Fill("Win_Type_"+Name,E.win_type[w]);
			H.Fill("Win_FlagVsType_"+Name,E.win_type[w],E.win_flag[w]);
			H.Fill("Win_PCWidthAvg_"+Name,E.win_type[w],E.win_pcwidthavg[w]);
			H.Fill("Win_PCtof_"+Name,E.win_type[w],E.win_pctof[w]);
			if (PerWindowType)
			{
				H.Fill("Win_Time_PerWinType_"+wtstr+"_"+Name,E.win_time[w]);
				H.Fill("Win_uvfirst_PerWinType_"+wtstr+"_"+Name,E.win_ufirst[w],E.win_vfirst[w]);
				H.Fill("Win_uvlast_PerWinType_"+wtstr+"_"+Name,E.win_ulast[w],E.win_vlast[w]);
				H.Fill("Win_pfirst_PerWinType_"+wtstr+"_"+Name,E.win_pfirst[w]);
				H.Fill("Win_pupvlast_PerWinType_"+wtstr+"_"+Name,E.win_pulast[w],E.win_pvlast[w]);
			}
		}

		//_________ Helix fit tracks
		for ( int t = 0; t < E.dkwintrack.size(); t++)
		{
			if (E.hefit_ndof[t] > 0)
			{
				H.Fill("Chi2DkWinTrk_"+Name,E.hefit_chi2[t] / E.hefit_ndof[t]);
				H.Fill("Chi2DkWinTrkFull_"+Name,E.hefit_chi2[t] / E.hefit_ndof[t]);
				H.Fill("NdofDkWinTrk_"+Name,E.hefit_ndof[t]);
				H.Fill("NdofDkWinTrkFull_"+Name,E.hefit_ndof[t]);
			}
		}

		for ( int t = 0; t < E.seltrack.size(); t++)
		{
			if (E.hefit_ndof[t] > 0)
			{
				H.Fill("Chi2SelTrk_"+Name,E.hefit_chi2[t] / E.hefit_ndof[t]);
				H.Fill("Chi2SelTrkFull_"+Name,E.hefit_chi2[t] / E.hefit_ndof[t]);
				H.Fill("NdofSelTrk_"+Name,E.hefit_ndof[t]);
				H.Fill("NdofSelTrkFull_"+Name,E.hefit_ndof[t]);
			}
		}

	return true;
}

