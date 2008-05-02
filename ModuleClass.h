//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ModuleClass_h
#define ModuleClass_h
// Include de C++
using namespace std;

#include <iostream>
#include <string>
#include <vector>
#include <map>

// Include Program
#include "EventClass.h"
#include "HistogramFactory.h"
#include "ConfigFile.h"
// Include log4cpp
#include "log4cpp/Category.hh"

class ModuleClass {
	public :
		ModuleClass()	{};
		~ModuleClass()	{};
		virtual bool	Init(EventClass &E, HistogramFactory &H, ConfigFile &Conf, log4cpp::Category *TmpLog) {};
		virtual bool	Process(EventClass &E, HistogramFactory &H) {};

		// map<const char*, TH1D*> H1D;
	private :

};



#endif
