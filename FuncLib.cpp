#include "FuncLib.h"

// template<class T>
// string IntToStr( const T& t )
string IntToStr( int &t )
{
	// Convert from a T to a string
	// Type T must support << operator
	std::ostringstream ost;
	ost << t;
	return ost.str();
}

vector<int> StrToIntVect(string In, char sep)
{
	vector<int>	Output;
	vector<string>TmpList = StrToStrVect(In, sep);

	for(vector<string>::iterator T=TmpList.begin(); T != TmpList.end(); T++)
		Output.push_back(atoi( (*T).c_str() ));

	return Output;
}

vector<float> StrToFloatVect(string In, char sep)
{
	vector<float>	Output;
	vector<string>TmpList = StrToStrVect(In, sep);

	for(vector<string>::iterator T=TmpList.begin(); T != TmpList.end(); T++)
		Output.push_back(atof( (*T).c_str() ));

	return Output;
}

vector<string> StrToStrVect(string In, char sep)
{
	vector<string>	Output;

	while (In.find(sep, 0) != string::npos)
	{
		Output.push_back(In.substr( 0, In.find(sep, 0)));
		In	= In.substr( In.find(sep, 0)+1);
	}
	Output.push_back(In);

	return Output;
}

string ReadParam( string Regexp, string InStr)
{
 	boost::regex e(Regexp);
	boost::smatch what;
	if(boost::regex_search(InStr, what, e, boost::match_default))
	{
		string Tmp = what[1];
		return Tmp;
	}
	else
	{
		return "";
	}
}
