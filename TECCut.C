//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Glen Marshall, December 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here

class TECCut : public ModuleClass{
	public :
		TECCut()		{};
		~TECCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		log4cpp::Category *Log;

		string TECifailxName;
		string TECifailyName;
		string TECxName;
		string TECyName;
		string TECtanthxName;
		string TECtanthyName;
		string TECsigmaxName;
		string TECsigmayName;

		double TECifailx_min;
		double TECifailx_max;
		double TECifaily_min;
		double TECifaily_max;
		double TECx_min;
		double TECx_max;
		double TECy_min;
		double TECy_max;
		double TECtanthx_min;
		double TECtanthx_max;
		double TECtanthy_min;
		double TECtanthy_max;
		double TECsigmax_min;
		double TECsigmax_max;
		double TECsigmay_min;
		double TECsigmay_max;

		double tec_x;
		double tec_y;
		double tec_tanthx;
		double tec_tanthy;
};

bool TECCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register TEC Cuts");

	// The module is stopped early to avoid the list of extra TEC cuts
	// to appear in the global histograms of the standard analysis.
	// Also the module is stopped if there is no TEC variables in the first tree.
	if ( (not Conf.read<bool>("TECCut/Do")) || (not (E.Exists("xy0") && E.Exists("tanthxy"))))
	{
		Log->info( "TEC Cuts turned OFF");
		return false ;
	}

	//    -------- Name of the cuts ---------     //
	TECifailxName	= "TECifailx";
	TECifailyName	= "TECifaily";
	TECxName	= "TECx";
	TECyName	= "TECy";
	TECtanthxName	= "TECtanthx";
	TECtanthyName	= "TECtanthy";
	TECsigmaxName	= "TECsigmax";
	TECsigmayName	= "TECsigmay";


	//	 --------- Special list of cuts in this class for the global histograms ---------		//
	H.AddCut(TECifailxName);
	H.AddCut(TECifailyName);
	H.AddCut(TECxName);
	H.AddCut(TECyName);
	H.AddCut(TECtanthxName);
	H.AddCut(TECtanthyName);
	H.AddCut(TECsigmaxName);
	H.AddCut(TECsigmayName);

	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "TECCut","TECifailx_before",	"TEC ifailx before the cut",40, -19.5,20.5);
	H.DefineTH1D( "TECCut","TECifailx_after",	"TEC ifailx after the cut",40, -19.5,20.5);
	H.DefineTH1D( "TECCut","TECifaily_before",	"TEC ifaily before the cut",40, -19.5,20.5);
	H.DefineTH1D( "TECCut","TECifaily_after",	"TEC ifaily after the cut",40, -19.5,20.5);

	H.DefineTH1D( "TECCut","TECx_before",	"TEC x before the cut",160, -4,4);
	H.DefineTH1D( "TECCut","TECx_after",	"TEC x after the cut",160, -4,4);
	H.DefineTH1D( "TECCut","TECy_before",	"TEC y before the cut",160, -4,4);
	H.DefineTH1D( "TECCut","TECy_after",	"TEC y after the cut",160, -4,4);

	H.DefineTH1D( "TECCut","TECtanthx_before",      "TEC tanthx before the cut",200, -0.2,0.2);
	H.DefineTH1D( "TECCut","TECtanthx_after",       "TEC tanthx after the cut",200, -0.2,0.2);
	H.DefineTH1D( "TECCut","TECtanthy_before",      "TEC tanthy before the cut",200, -0.2,0.2);
	H.DefineTH1D( "TECCut","TECtanthy_after",	"TEC tanthy after the cut",200, -0.2,0.2);

	H.DefineTH1D( "TECCut","TECsigmax_before",	"TEC sigmax before the cut",100, 0,0.5);
	H.DefineTH1D( "TECCut","TECsigmax_after",	"TEC sigmax after the cut",100, 0,0.5);
	H.DefineTH1D( "TECCut","TECsigmay_before",	"TEC sigmay before the cut",100, 0,0.5);
	H.DefineTH1D( "TECCut","TECsigmay_after",	"TEC sigmay after the cut",100, 0,0.5);


	H.DefineTH2D( "TECCut","TECyvsx_before", "TEC y vs x before the cut",160, -4,4,160, -4,4);
	H.DefineTH2D( "TECCut","TECyvsx_after",	"TEC y vs x after the cut",160, -4,4,160, -4,4);
	H.DefineTH2D( "TECCut","TECtanthyvstanthx_before", "TEC tanthy vs tanthx before the cut",200, -0.2,0.2, 200, -0.2,0.2);
	H.DefineTH2D( "TECCut","TECtanthyvstanthx_after", "TEC tanthy vs tanthx after the cut",200, -0.2,0.2, 200, -0.2,0.2);

	//	 --------- Parameters initialization ---------		//
	TECx_min		= Conf.read<double>("TECCut/TECxMin");
	TECx_max		= Conf.read<double>("TECCut/TECxMax");
	TECy_min		= Conf.read<double>("TECCut/TECyMin");
	TECy_max		= Conf.read<double>("TECCut/TECyMax");

	TECtanthx_min		= Conf.read<double>("TECCut/TECtanthxMin");
	TECtanthx_max		= Conf.read<double>("TECCut/TECtanthxMax");
	TECtanthy_min		= Conf.read<double>("TECCut/TECtanthyMin");
	TECtanthy_max		= Conf.read<double>("TECCut/TECtanthyMax");

	TECsigmax_min		= Conf.read<double>("TECCut/TECsigmaxMin");
	TECsigmax_max		= Conf.read<double>("TECCut/TECsigmaxMax");
	TECsigmay_min		= Conf.read<double>("TECCut/TECsigmayMin");
	TECsigmay_max		= Conf.read<double>("TECCut/TECsigmayMax");

	TECifailx_min		= Conf.read<double>("TECCut/TECifailxMin");
	TECifailx_max		= Conf.read<double>("TECCut/TECifailxMax");
	TECifaily_min		= Conf.read<double>("TECCut/TECifailyMin");
	TECifaily_max		= Conf.read<double>("TECCut/TECifailyMax");

	return true;
}

bool TECCut::Process(EventClass &E, HistogramFactory &H)
{
	// ____ TEC cut ____ //

	//1D Histograms are filled successively at the time of their respective cuts
	//2D histograms are filled cumulatively, before and after all cuts

	tec_tanthx = E.tanthxy[0];
	tec_tanthy = E.tanthxy[1];
	tec_x = tec_tanthx*(-191.944)+E.xy0[0];
	tec_y = tec_tanthy*(-191.944)+E.xy0[1];

	H.Fill("TECyvsx_before",tec_x,tec_y);
	H.Fill("TECtanthyvstanthx_before",tec_tanthx,tec_tanthy);

	H.NbCandidateTracks(TECifailxName,E);
	H.Fill("TECifailx_before",E.ifail[0]);
	if ( E.ifail[0] > TECifailx_min && E.ifail[0] < TECifailx_max )
	{
		H.Fill("TECifailx_after",E.ifail[0]);
	}
	else
	{
		H.CutApplied(TECifailxName);
		return false;
	}

	H.NbCandidateTracks(TECifailyName,E);
	H.Fill("TECifaily_before",E.ifail[1]);
	if ( E.ifail[1] > TECifaily_min && E.ifail[1] < TECifaily_max )
	{
		H.Fill("TECifaily_after",E.ifail[1]);
	}
	else
	{
		H.CutApplied(TECifailyName);
		return false;
	}

	H.NbCandidateTracks(TECxName,E);
	H.Fill("TECx_before",tec_x);
	if ( tec_x > TECx_min && tec_x < TECx_max )
	{
		H.Fill("TECx_after",tec_x);
	}
	else
	{
		H.CutApplied(TECxName);
		return false;
	}

	H.NbCandidateTracks(TECyName,E);
	H.Fill("TECy_before",tec_y);
	if ( tec_y > TECy_min && tec_y < TECy_max )
	{
		H.Fill("TECy_after",tec_y);
	}
	else
	{
		H.CutApplied(TECyName);
		return false;
	}

	H.NbCandidateTracks(TECtanthxName,E);
	H.Fill("TECtanthx_before",tec_tanthx);
	if ( tec_tanthx > TECtanthx_min && tec_tanthx < TECtanthx_max )
	{
		H.Fill("TECtanthx_after",tec_tanthx);
	}
	else
	{
		H.CutApplied(TECtanthxName);
		return false;
	}

	H.NbCandidateTracks(TECtanthyName,E);
	H.Fill("TECtanthy_before",tec_tanthy);
	if ( tec_tanthy > TECtanthy_min && tec_tanthy < TECtanthy_max )
	{
		H.Fill("TECtanthy_after",tec_tanthy);
	}
	else
	{
		H.CutApplied(TECtanthyName);
		return false;
	}

	H.NbCandidateTracks(TECsigmaxName,E);
	H.Fill("TECsigmax_before",E.sigma[0]);
	if ( E.sigma[0] > TECsigmax_min && E.sigma[0] < TECsigmax_max )
	{
		H.Fill("TECsigmax_after",E.sigma[0]);
	}
	else
	{
		H.CutApplied(TECsigmaxName);
		return false;
	}

	H.NbCandidateTracks(TECsigmayName,E);
	H.Fill("TECsigmay_before",E.sigma[1]);
	if ( E.sigma[1] > TECsigmay_min && E.sigma[1] < TECsigmay_max )
	{
		H.Fill("TECsigmay_after",E.sigma[1]);
	}
	else
	{
		H.CutApplied(TECsigmayName);
		return false;
	}

	H.Fill("TECyvsx_after",tec_x,tec_y);
	H.Fill("TECtanthyvstanthx_after",tec_tanthx,tec_tanthy);

	return true;
}

