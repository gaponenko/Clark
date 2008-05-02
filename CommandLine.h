/******************************************************************************************************************/
// Fonction de lecture des arguments de la ligne de commande							  */
/******************************************************************************************************************/

#ifndef CommandLine_h
#define CommandLine_h
//
// Include de C++
#include <iostream>
#include <getopt.h>
using namespace std;
#include <math.h>
#include <string.h>
#include <sstream>
#include <fstream>

#include <TSystem.h>
#include <TROOT.h>

#include "ConfigFile.h"

class CommandLine
{
	public:
		CommandLine			();
		~CommandLine		(){};
		int		Init		(int argc, char **argv);
		void	Print_Help	();
		void	SaveToConfig(ConfigFile &Conf);

		string			GetConfFile	()			{return ConfFile;				};
		vector<string>	GetTrees	()			{return Trees;					};

	private:
		vector<string>	Trees;
		string			ConfFile;
		string			Executable;
		bool			Debug;
		string			EcalFile;
		string			EcalArray;
		string			LogFile;
		string			OutputFile;
};

#endif
