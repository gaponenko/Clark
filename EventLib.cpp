#include "EventLib.h"

// Simply returns the distance in plane number between the track and the target.
// For a track starting at the plane 22, this function would return 0.
// For a track starting at the plane 21, this function would return 1.
int DistanceToTarget( EventClass &E, int T)
{
		return int(min ( fabs(E.dcmin[T] - (E.NumDC / 2) - 0.5), fabs(E.dcmax[T] - (E.NumDC / 2) - 0.5) ));
}


// Returns the u and v position of the track at z
void Get_uv_at( const EventClass *E, int Trk, double z, double &uu, double &vv)
{
	// L	= -2 * M_PI * 8.4 * 2.0 / E.BField * ( 1/50.0 ) * E.hefit_pz[Trk]	// From tta, TrackParConverter.h and HFTrack.h. Sign very tough to extract.
	// phi = ( 2 * M_PI ( z - z0 ) ) / L										// From tta, Helix.h
	// phi = (z0 - z) / ( 8.4 * 2.0 / BField * ( 1/50.0 ) * -pz )

	// Radius
	// double R	= 8.4 * ( 2.0 / E.BField ) * ( 1/50.0 ) * E.hefit_pt[Trk];	// From tta, TrackParConverter.h
	// Wavelength
	// double L	= 2.0 * M_PI * R * (-1*E.hefit_q[Trk]*E.hefit_pz[Trk] / E.hefit_pt[Trk]);

	double u0	= E->hefit_u[Trk] - E->radius[Trk] * cos(E->hefit_phi[Trk]);
	double v0	= E->hefit_v[Trk] - E->radius[Trk] * sin(E->hefit_phi[Trk]);
	double z0	= E->hefit_z[Trk] - E->hefit_phi[Trk] * E->wavelen[Trk] / ( 2.0 * M_PI);
	double phi	= 2.0 * M_PI * ( z - z0 ) / E->wavelen[Trk];

	uu	= u0 + E->radius[Trk] * cos(2.0 * M_PI * ( z - z0 ) / E->wavelen[Trk]);
	vv	= v0 + E->radius[Trk] * sin(2.0 * M_PI * ( z - z0 ) / E->wavelen[Trk]);
}

double MichelWeight(unsigned accflag) {
  switch(accflag) {
  case 0x1: case 0x4: case 0x10: case 0x40: case 0x100: case 0x400: case 0x1000:
    return +1; 
  case 0x3: case 0xc: case 0x30: case 0xc0: case 0x300: case 0xc00: case 0x3000:
    return -1; 
  case 0: // ignore event-decay cache exhausted
    return 0;
  default:
    throw "MichelWeight: unknown accflag";
  }
}

// double DistanceOfApproach( double z, EventClass &E, int t, int a, double zmin, double zmax):
// {
// 	double ut, vt, ua, va;
// 	// Track
// 	Get_uv_at(E, t, z, ut, vt);
// 	// Anti-track
// 	Get_uv_at(E, a, z, ua, va);
// 
// 	// This is an exponential punishment to ensure that the minimum is found in the range we are looking in.
// 	if( z < zmin)
// 		return (sqrt( (ut-ua)*(ut-ua) + (vt-va)*(vt-va) ) + exp ( abs(zmin - z) ));
// 	else if( z > zmax)
// 		return (sqrt( (ut-ua)*(ut-ua) + exp ( abs(zmax - z) )));
// 	else
// 		return (sqrt( (ut-ua)*(ut-ua) + (vt-va)*(vt-va) ));
// }


