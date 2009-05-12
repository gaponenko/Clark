//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Includes here 
// /usr/include/gsl
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include "EventLib.h"

class PairMatchingCut : public ModuleClass{
	public :
		PairMatchingCut()		{};
		~PairMatchingCut()		{};
		bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) ;
		bool	Process(EventClass &E, HistogramFactory &H);

	private :
		void Get_Anti_Tracks( EventClass &E, int t);
		bool Get_CDA_FromTree( EventClass &E, int t, int a, double &z, double &cda, double &defl);
		bool Get_CDA_FromCalculation( EventClass &E, int t, int a, double &z, double &cda, int &iterations);

		log4cpp::Category *Log;

		string Name;
		double MinCDA;
		double MinDT;
		vector<int>AntiTracks;

		double dT;

		double tolerance_abs;
		double tolerance_rel;
		int max_iter;
		bool CalculateCDA;

		float MinCDA_Beame;
		float MinCDA_BrokTrk;
		vector<float> BrokTrk_Z_Range;
		vector<float> Beame_Z_Range;
};

bool PairMatchingCut::Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog)
{
	Log	= TmpLog;
	Log->info( "Register Pair Matching Cut");
	//    -------- Name of the cut ---------     //
	Name		= "pair matching";


	//	 --------- Special list of cut in this class for the global histograms ---------		//
	H.AddCut(Name);
	
	if (not Conf.read<bool>("PairMatchingCut/Do"))
	{
		Log->info( "Pair Matching Cut turned OFF");
		return false ;
	}
	
	//	 --------- Histograms initialization ---------		//
	H.DefineTH1D( "PairMatchingCut", "ntrks_cda_before",	"Number of tracks before the CDA cut",10, -0.5,9.5);
	H.DefineTH1D( "PairMatchingCut", "ntrks_cda_after",		"Number of tracks after the CDA cut", 10, -0.5,9.5);
	H.DefineTH1D( "PairMatchingCut", "cda_up_before",		"Closest Distance of Approach before the cut, upstream tracks;CDA [cm]",	400, 0.0,10.0);
	H.DefineTH1D( "PairMatchingCut", "cda_up_after",		"Closest Distance of Approach after the cut, upstream tracks;CDA [cm]",		400, 0.0,10.0);
	H.DefineTH1D( "PairMatchingCut", "cda_down_before",		"Closest Distance of Approach before the cut, downstream tracks;CDA [cm]",	400, 0.0,10.0);
	H.DefineTH1D( "PairMatchingCut", "cda_down_after",		"Closest Distance of Approach after the cut, downstream tracks;CDA [cm]",	400, 0.0,10.0);
	H.DefineTH2D( "PairMatchingCut", "dt_vs_cda_before",	"T(more us) - T(more ds) vs CDA before the cut, downstream tracks;CDA [cm];dt [ns]",200, 0.0,10.0, 100, -500., 500.);
	H.DefineTH2D( "PairMatchingCut", "dt_vs_cda_after",		"T(more us) - T(more ds) vs CDA after the cut, downstream tracks;CDA [cm];dt [ns]",	200, 0.0,10.0, 100, -500., 500.);

	// if( self.CalcVsMofia)
	// {
	// 	H.DefineTH1D( "PairMatchingCut", "Calc-MOFIA_up_cda",			"Difference in cda between calculation and MOFIA, pair matching cut, upstream tracks", 2, -0.5, 1.5);
	// 	H.DefineTH1D( "PairMatchingCut", "Calc-MOFIA_down_cda",			"Difference in cda between calculation and MOFIA, pair matching cut, downstream tracks", 2, -0.5, 1.5);
	// 	H.DefineTH1D( "PairMatchingCut", "Calc-MOFIA_up_cda_agree",		"Agreement between calculation and MOFIA, pair matching cut, upstream tracks;Agreement (agree=1)", 2, -0.5, 1.5);
	// 	H.DefineTH1D( "PairMatchingCut", "Calc-MOFIA_down_cda_agree",	"Agreement between calculation and MOFIA, pair matching cut, downstream tracks;Agreement (agree=1)", 2, -0.5, 1.5);
	// }

	H.DefineTH1D( "PairMatchingCut", "cda_z_min",		"Z of the minimum closest approach found",600, -60., 60.);
	H.DefineTH1D( "PairMatchingCut", "cda_z_matching",	"Z of the closest approach for cut tracks",600, -60., 60.);

	H.DefineTH1D( "PairMatchingCut", "cda_cuttype_found",	"Cut type found from the zmin",3, -0.5, 2.5);
	H.DefineTH1D( "PairMatchingCut", "cda_cuttype_applied",	"Cut type applied",3, -0.5, 2.5);

	// if (not(E.Exists("hefit_cda")))
	// 	H.DefineTH1D( "PairMatchingCut", "anti-track_cut",	"Cuts applied on anti-tracks", 3, 0.5, 3.5);

	//	 --------- Parameters initialization ---------		//
	MinCDA			= Conf.read<double>("PairMatchingCut/MinCDA");
	MinCDA_Beame	= Conf.read<double>("PairMatchingCut/MinCDA_Beame");
	MinCDA_BrokTrk	= Conf.read<double>("PairMatchingCut/MinCDA_BrokTrk");
	MinDT			= Conf.read<double>("PairMatchingCut/MinDT");

	tolerance_abs	= Conf.read<double>("PairMatchingCut/tolerance_abs");	// cm
	tolerance_rel	= Conf.read<double>("PairMatchingCut/tolerance_rel");
	max_iter		= Conf.read<int>("PairMatchingCut/max_iter");
	CalculateCDA	= Conf.read<bool>("PairMatchingCut/CalculateCDA");

	BrokTrk_Z_Range	= StrToFloatVect(Conf.read<string>("PairMatchingCut/BrokTrk_Z_Range"));
	Beame_Z_Range	= StrToFloatVect(Conf.read<string>("PairMatchingCut/Beame_Z_Range"));

	return true;
}

bool PairMatchingCut::Process(EventClass &E, HistogramFactory &H)
{
	// cout<<" In PairMatchingCut 1   event = "<<E.nevt<<endl;
	// ____ ____ //
	H.NbCandidateTracks(Name,E);
	
	dT = 100000.;
	double cda, zmin, deflang; // deflection angle
	int iter;	// number of iterations for the minimzation
	bool IsCDAok;
	float CDAcut;
	int CDAcuttype;
	// cout<<" In PairMatchingCut 2"<<endl;
	
	for(vector<int>::iterator t = E.seltrack.begin(); t != E.seltrack.end(); t++)
	{
		// First, get the antitracks			
		Get_Anti_Tracks( E, *t);
		// cout<<" In PairMatchingCut 3  AntiTracks size = "<<AntiTracks.size()<<endl;
		for(vector<int>::iterator a = AntiTracks.begin(); a != AntiTracks.end(); a++)
		{
			// Calculate Delta t between the track and the best anti-track
			dT	= E.hefit_time[*t] - E.hefit_time[*a];
			// Which one is more downstream ?
			// (From tta ... not sure why it is done. Maybe just for the histograms)
			if( E.dcmin[*t] >= E.dcmin[*a])
				dT = -1 * dT;

			// cout<<" In PairMatchingCut 4"<<endl;
			//
			//	Read from the tree or calculate the CDA.
			if ( ( ! E.Exists("dkwin_ncda") ) || CalculateCDA )
				IsCDAok = Get_CDA_FromCalculation(E, (*t), (*a), zmin, cda, iter);
			else
				IsCDAok = Get_CDA_FromTree(E, (*t), (*a), zmin, cda, deflang);

			if (not IsCDAok )
			{
				// Something went wrong !
				Log->warn("PairMatchingCut: No cda cut for event %i, track %i, anti-track%i. The time cut is still applied though.", E.nevt, *t, *a);
				if( fabs(dT) < MinDT)
				{
					E.seltrack.erase(t);	// First erase
					t--;					// then decrement to avoid to skip the following track
					H.Fill("cda_z_matching",zmin);
					break;	// from anti-track loop
				}
				continue; // go to next anti-track
			}
			// cout<<" In PairMatchingCut 5   cda = "<<cda<<"   z = "<<zmin<<"    dt = "<<dT<<endl;
		
			H.Fill("cda_z_min",zmin);
			if( E.is_upstreamdk)
				H.Fill("cda_up_before",cda);
			else
				H.Fill("cda_down_before",cda);
		
			H.Fill("dt_vs_cda_before",cda, dT);
		
			// The anti track number is less than zero if the cut must not be applied.
			
			if( BrokTrk_Z_Range[0] < fabs(zmin) && fabs(zmin) < BrokTrk_Z_Range[1] )
			{
				CDAcut = MinCDA_BrokTrk;
				CDAcuttype = 1;
			}
			else if( Beame_Z_Range[0] < fabs(zmin) && fabs(zmin) < Beame_Z_Range[1] )
			{
				CDAcut = MinCDA_Beame;
				CDAcuttype = 2;
			}
			else
			{
				CDAcut = MinCDA;
				CDAcuttype = 0;
			}

			H.Fill("cda_cuttype_found",CDAcuttype);
			if( cda < CDAcut and fabs(dT) < MinDT)
			{
				E.seltrack.erase(t);	// First erase
				t--;					// then decrement to avoid to skip the following track
				H.Fill("cda_z_matching",zmin);
				H.Fill("cda_cuttype_applied",CDAcuttype);
				break;	// from anti-track loop
			}
			else
			{
				if( E.is_upstreamdk)
					H.Fill("cda_up_after",cda);
				else
					H.Fill("cda_down_after",cda);
				H.Fill("dt_vs_cda_after",cda, dT);
			}
		}
	}

	if (E.seltrack.size() < 1 )
	{
		H.CutApplied(Name);
		return false;
	}
	return true;
}


// =================================================================== //

void PairMatchingCut::Get_Anti_Tracks( EventClass &E, int t)
{
	// Selection of the anti-tracks with the following criteria:
	// - All tracks in decay window with an ierror of 0
	// - Anti-track doesn't overlap track candidate
	// - Anti-track candidate is closer to the target than track candidate.
	AntiTracks.clear();
	// cout<<" In Get_Anti_Tracks 1  AntiTracks size = "<<AntiTracks.size()<<endl;
	// cout<<" In Get_Anti_Tracks 1  dkwintrack size = "<<E.dkwintrack.size()<<endl;
	// First, get the antitracks, i.e. all the tracks in the decay window
	// minus the analysed track
	for (int j = 0; j < E.dkwintrack.size(); j++)
	{
	 	// cout<<" In Get_Anti_Tracks 1.1  DkT = "<<E.dkwintrack[j]<<"   t = "<<t<<endl;
	 	if( E.dkwintrack[j] != t)
	 	{
	 		// cout<<" In Get_Anti_Tracks 1.2  Selected "<<t<<endl;
	 		AntiTracks.push_back(E.dkwintrack[j]);
	 	}
	}
	// for(vector<int>::iterator DkT = E.dkwintrack.begin(); DkT != E.dkwintrack.end(); DkT++)
	// {
	// 	cout<<" In Get_Anti_Tracks 1.1  DkT = "<<*DkT<<"   t = "<<t<<endl;
	// 	if( *DkT != t)
	// 	{
	// 		cout<<" In Get_Anti_Tracks 1.2  Selected "<<t<<endl;
	// 		AntiTracks.push_back(*DkT);
	// 	}
	// }
	
	// cout<<" In Get_Anti_Tracks 2  AntiTracks size = "<<AntiTracks.size()<<endl;
	// Secondly, is there an overlap ?
	// If so, remove anti-track
	int mindcmax, maxdcmin;
	for(vector<int>::iterator a = AntiTracks.begin(); a != AntiTracks.end(); a++)
	{
		maxdcmin = int(max(E.dcmin[t],E.dcmin[*a]));
		mindcmax = int(min(E.dcmax[t],E.dcmax[*a]));
		if( ( maxdcmin - mindcmax ) < 0)
		{
			AntiTracks.erase(a);	// First erase
			a--;					// then decrement to avoid to skip the following track
		}
	}
	// cout<<" In Get_Anti_Tracks 3  AntiTracks size = "<<AntiTracks.size()<<endl;

	// Remove the anti-tracks that are farther from the target than the track
	for(vector<int>::iterator a = AntiTracks.begin(); a != AntiTracks.end(); a++)
		if( DistanceToTarget( E, t) < DistanceToTarget( E, *a))
		{
			AntiTracks.erase(a);	// First erase
			a--;					// then decrement to avoid to skip the following track
		}
	// cout<<" In Get_Anti_Tracks 4  AntiTracks size = "<<AntiTracks.size()<<endl;

	// Helicity checked here eventually

}

// =================================================================== //



// =================================================================== //
// Convert the track tree index into good track (ierror = 0) index.		//
// =================================================================== //

int ConvertToGoodTrack( EventClass &E, int &In)
{
	for ( int i = 0; i < E.dkwintrack.size(); i++ )
	{
		if ( E.dkwintrack[i] == In )
			return i;
	}
	// If we reach that point, that's not good ...
	return -1;
}


// =================================================================== //
// Read the CDA from the tree.											//
// =================================================================== //

bool PairMatchingCut::Get_CDA_FromTree( EventClass &E, int t, int a, double &z, double &cda, double &defl)
{
	// 1) convert the track index into the dkwintrack index.
	int T1 = ConvertToGoodTrack( E, t);
	int T2 = ConvertToGoodTrack( E, a);
	if ( T1 == -1 || T2 == -1 )
	{
		Log->warn("PairMatchingCut: A track or anti-track is not a good track (ierror = 0) so there is no CDA in the tree for it.");
		return false;
	}
	// 2) order the indices
	if ( T1 > T2 )
	{
		int Temp = T1;
		T1 = T2;
		T2 = Temp;
	}
	// 3) get the cda index
	int CDA_i = T1 * E.dkwintrack.size() + T2 - 1;
	cda		= E.dkwin_cda[CDA_i];
	z		= E.dkwin_cdaz[CDA_i];
	defl	= E.dkwin_cdadefl[CDA_i];

	return true;
}


// =================================================================== //



// =================================================================== //
// Class wrapper for the minimization.                                  //
// =================================================================== //
class TrackDistanceUV {
	EventClass *Evt;
	int Trka;
	int Trkb;
	double zz1;
	double zz2;
	public:
	TrackDistanceUV(EventClass &E, int Ta, int Tb, double z1, double z2) : Evt(&E), Trka(Ta), Trkb(Tb), zz1(z1), zz2(z2) {}
	double operator()(double z) const;
	static double fwrapper(double z, void *);
};

double TrackDistanceUV::operator()(double z) const {
	double au, av, bu, bv;
	Get_uv_at(Evt, Trka, z, au, av);
	Get_uv_at(Evt, Trkb, z, bu, bv);
	return (sqrt( (au-bu)*(au-bu) + (av-bv)*(av-bv) ));
}

double TrackDistanceUV::fwrapper(double z, void *ptr) {
	return static_cast<TrackDistanceUV*>(ptr)->operator()(z);
}


struct MinBracket { 
  // minimum is around m, on [a,b]
  double a; double m; double b; 
  MinBracket(double A,double M,double B) : a(A), m(M), b(B) {}
};


// =================================================================== //
// Calculates the cda for the track t and the anti-track a.            //
// =================================================================== //

bool PairMatchingCut::Get_CDA_FromCalculation( EventClass &E, int t, int a, double &z, double &cda, int &iterations)
{
	//=====================================================================//

	double z1,z2;
	//------------------------------------------------------------------//
	// Find the range where to look for the CDA                         //
	//------------------------------------------------------------------//
	if( E.dcmin[t] < E.dcmin[a])
	{
		z1 = E.geo->zdplane[int(min(E.dcmax[t], E.dcmin[a])-1)];
		z2 = E.geo->zdplane[int(max(E.dcmax[t], E.dcmin[a])-1)];
	}
	else
	{
		z1 = E.geo->zdplane[int(min(E.dcmax[a], E.dcmin[t])-1)];
		z2 = E.geo->zdplane[int(max(E.dcmax[a], E.dcmin[t])-1)];
	}
	z1 = z1 -0.2;
	z2 = z2 +0.2;

	// L corresponding to the maximum frequency:
	// L1 --> omega1, L2 --> omega2, omega_max = w1 + w2 --> Lmin

	double lmin = fabs(E.wavelen[t])*fabs(E.wavelen[a])/(fabs(E.wavelen[t])+fabs(E.wavelen[a]));
	// cout<<" In Get_CDA_FromCalculation 2   lmin ="<<lmin<<endl;
	// cout<<" In Get_CDA_FromCalculation 2   wavelen ="<<E.wavelen[t]<<"   "<<E.wavelen[a]<<endl;
	// cout<<" In Get_CDA_FromCalculation 2   ptot ="<<E.ptot[t]<<"   "<<E.costh[a]<<endl;
	//------------------------------------------------------------------//
	// Sample the distance to isolate local minima                      //
	// First decide about the number of points in the sample            //
	//------------------------------------------------------------------//
	int nsample = int(12*(1.+fabs(z1-z2)/lmin));

	TrackDistanceUV distance(E, t, a, z1, z2);
	vector<double> dd(nsample);
	for(int i = 0; i<nsample; i++) {
		double z = z1 + (z2-z1)*i/double(nsample-1);
		dd[i] = distance(z);
	}

	double A,M,B;
	vector<MinBracket> brackets;
	for(int i=1; i<nsample-1; i++) {
		if( (dd[i] < dd[i-1]) && (dd[i] <= dd[i+1]) ) {
			A = z1 + (z2-z1)*(i-1)/double(nsample-1);
			M = z1 + (z2-z1)*i/double(nsample-1);
			B = z1 + (z2-z1)*(i+1)/double(nsample-1);
			brackets.push_back(MinBracket(A,M,B));
		}
	}

	//------------------------------------------------------------------//
	// Check the ends of the search interval                            //
	//------------------------------------------------------------------//

	// cout<<" In Get_CDA_FromCalculation 3"<<endl;
	// The distance of approach
	vector<double> Alldist;
	// Position of minimal distance of approach
	vector<double> Allz;
	vector<int> Alliter;

	if(dd[0] <= dd[1]) {
		double ztest = z1 + std::min(z2-z1, tolerance_abs);
		double ftest = distance(ztest);
		if( (ftest <= dd[0]) && (ftest <= dd[1])) {
			//std::cerr<<"                z1 bracket"<<std::endl;
			brackets.push_back(MinBracket(z1,ztest, z1 + (z2-z1)/double(nsample-1)));
		}
		else {
			//std::cerr<<"                z1 min"<<std::endl;
			// Extract the distances of approach here
			Alldist.push_back(dd[0]);
			// Position of minimal distance of approach
			Allz.push_back(z1);
			Alliter.push_back(0);
		}
	}

	if(dd[nsample-1] <= dd[nsample-2]) {
		double ztest = z2 - std::min(z2-z1, tolerance_abs);
		double ftest = distance(ztest);
		if( (ftest <= dd[nsample-1]) && (ftest <= dd[nsample-2])) {
			//std::cerr<<"                z2 bracket"<<std::endl;
			brackets.push_back(MinBracket(z1 + (z2-z1)*(nsample-2)/double(nsample-1), ztest, z2));
		}
		else {
			//std::cerr<<"                z2 min"<<std::endl;
			// Extract the distances of approach here
			Alldist.push_back(dd[nsample-1]);
			// Position of minimal distance of approach
			Allz.push_back(z2);
			Alliter.push_back(0);
		}
	}
	// cout<<" In Get_CDA_FromCalculation 4   brackets size ="<<brackets.size()<<endl;

	//------------------------------------------------------------------//
	// Loop over the brackets and for each one of them, find a minimum. //
	//------------------------------------------------------------------//

	for(unsigned i=0; i<brackets.size(); i++) {
		//std::cerr<<"Bracket{"<<brackets[i].a<<", "<<brackets[i].m<<", "<<brackets[i].b<<std::endl;
		const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
		gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);

		gsl_function F;

		F.function = &TrackDistanceUV::fwrapper;
		F.params = &distance;

		gsl_min_fminimizer_set(s, &F, brackets[i].m, brackets[i].a, brackets[i].b);

		int iter=0;
		while(1)
		{
			iter++;
			int status;
			if(GSL_SUCCESS != (status = gsl_min_fminimizer_iterate (s)) )
			{
				gsl_min_fminimizer_free(s);
				Log->warn("PairMatchingCut: Get_CDA_FromCalculation(): failure from gsl_min_fminimizer_iterate(): status==%i", status);
				return false;
			}

			M = gsl_min_fminimizer_x_minimum (s);
			A = gsl_min_fminimizer_x_lower (s);
			B = gsl_min_fminimizer_x_upper (s);

			status = gsl_min_test_interval (A, B, tolerance_abs, tolerance_rel);

			if (status == GSL_SUCCESS) { 
				//printf ("Converged:\n");
				//printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, a, b, m, b - a);

				// Extract the distances of approach here
				Alldist.push_back(gsl_min_fminimizer_f_minimum(s));
				// Position of minimal distance of approach
				Allz.push_back(M);
				Alliter.push_back(iter);
				break;
			}
			else if(status == GSL_CONTINUE) {
				if(iter>max_iter) {
					gsl_min_fminimizer_free(s);
					Log->warn("PairMatchingCut: Get_CDA_FromCalculation(): exceeded max_iter=");
					return false;
				}
			}
			else {
				gsl_min_fminimizer_free(s);
				Log->warn("PairMatchingCut: Get_CDA_FromCalculation(); failure from gsl_min_test_interval(): status=");
				return false;
			}	

		} // while(1) minimize

		gsl_min_fminimizer_free(s);

	} // for(brackets)

	// cout<<" In Get_CDA_FromCalculation 5    Alldist size ="<<Alldist.size()<<endl;
	//------------------------------------------------------------------//
	// The CDA is the minimum of all the distances of approach          //
	//------------------------------------------------------------------//
	int MinIndex = -1;
	cda = 500.;
	z	= 1000.;
	iterations = 0;
	// cout<<" In Get_CDA_FromCalculation 6   cda = "<<cda<<endl;
	for (int i = 0; i < Alldist.size(); i++)
	{
		// cout<<" In Get_CDA_FromCalculation 6.1   AllDist[0] = "<<Alldist[0]<<endl;
		if ( Alldist[i] < cda )
		{
			MinIndex = i;
			cda	= Alldist[i];
			// cout<<" In Get_CDA_FromCalculation 6.2   cda = "<<cda<<endl;
		}
	}
	if (MinIndex == -1)
	{
		Log->warn("PairMatchingCut: Get_CDA_FromCalculation(); The cda could not be found.");
		cda	= 500;	// Safer with a large number since we cut on the small cda.
		return false;
	}
	else
	{
		z = Allz[MinIndex];
		iterations = Alliter[MinIndex];
	}
	// cout<<" In Get_CDA_FromCalculation 7   cda = "<<cda<<endl;
	
	return true;
}


