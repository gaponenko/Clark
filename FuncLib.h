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
#include <sstream>
#include <boost/regex.hpp>

#ifndef FuncLib_h
#define FuncLib_h

// Split a string wrt sep and store the result in a vector of various types
string IntToStr( int &t );
vector<int>		StrToIntVect	(string In, char sep = ',');
vector<float>	StrToFloatVect	(string In, char sep = ',');
vector<string>	StrToStrVect	(string In, char sep = ',');

string ReadParam( string Regexp, string InStr);

#endif
