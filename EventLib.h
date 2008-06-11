//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Include de C++
using namespace std;

#include <string>
#include <vector>
#include <iostream>

#include "EventClass.h"

#ifndef EventLib_h
#define EventLib_h

int DistanceToTarget( EventClass &E, int T);
void Get_uv_at( const EventClass *E, int Trk, double z, double &uu, double &vv);
int MichelWeight( unsigned accflag);

#endif

