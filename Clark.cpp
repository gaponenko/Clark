//////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Anthony Hillairet, April 2008
///
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Include de C++
using namespace std;

#include <iostream>
#include <string.h>
#include <vector>

// Include Program
#include "TreeClass.h"
#include "CommandLine.h"
#include "ConfigFile.h"

// Include log4cpp
#include "log4cpp/Category.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/PatternLayout.hh"

// See Default.C
void SetDefault(ConfigFile &Conf);

// See Users.C
void LoadAnalysisClasses( TreeClass *AnaObj);

int main(int argc, char **argv)
{
	// ================ COMMAND LINE ===================== //

	CommandLine	Com;

	Com.Init (argc, argv);

	// ================ CONFIG FILE ===================== //

	string ConfFilename	= Com.GetConfFile();

	// Check that the config file can be opened correctly
	try
	{
	  ConfigFile TmpConf(ConfFilename.c_str());
	}
	catch( const ConfigFile::file_not_found &e)
	{
		cerr<<" Configuration file \""<<e.filename<<"\" NOT FOUND."<<endl;
		exit(1);
	}

	ConfigFile Conf(ConfFilename.c_str());
	// Name of the analysis from the name of the config filename
	Conf.add( "Name", ConfFilename );

	Com.SaveToConfig(Conf);

	// Set the default values in DefaultConfig.C
	SetDefault(Conf);

	// ================ LOG FILE ===================== //
	
	string LogFile = Conf.read<string>("CommandLine/LogFile");

	// The default for the logfile name is the configuration filename plus ".log"
	if (LogFile == "")
		LogFile	= ConfFilename + ".log";

	// 1 instantiate an appender object that 
	// will append to a log file or stdout
	log4cpp::Appender* Fapp = new log4cpp::FileAppender("FileAppender", LogFile.c_str(), false);
	log4cpp::Appender* Capp = new log4cpp::FileAppender("ConsoleAppender", "/dev/stdout", true);
	// 2. Instantiate a layout object
	log4cpp::PatternLayout* playout = new log4cpp::PatternLayout();
	playout->setConversionPattern("%d{%Y-%m-%d %H:%M:%S} %p: %m%n");

	// 3. attach the layout object to the appender object
	Fapp->setLayout(playout);
	Capp->setLayout(playout);

	// 4. Instantiate the category object
	log4cpp::Category &LogAll = log4cpp::Category::getRoot();

	// 5. Set the priorities and add the appenders
	if (LogFile != "stdout")
	{
		Fapp->setThreshold(log4cpp::Priority::INFO);
		Capp->setThreshold(log4cpp::Priority::WARN);
		LogAll.addAppender(Fapp);
	}
	else
		Capp->setThreshold(log4cpp::Priority::INFO);
	// PROBLEM: can't use that ... Clark crashes at the end.
	// Must be commented for now.
	//	LogAll.addAppender(Capp);


	LogAll.info("######## Starting treesumming with %s ########\n", argv[0]);

	// ================ TREE ANALYSIS ===================== //

	vector<string> TreeFiles = Com.GetTrees();

	TreeClass AnaObj;

	LoadAnalysisClasses( &AnaObj);


	LogAll.info("Initialization of the TreeAnalysis class");
	if (not AnaObj.InitAll(TreeFiles[0], Conf, &LogAll))
	{
		LogAll.error("Problem at the initialization with the first tree. Exit program.");
		exit(1);
	}

	// Save the configuration.
	LogAll << log4cpp::Priority::INFO <<"\n"<< Conf;

	for(vector<string>::iterator F=TreeFiles.begin(); F != TreeFiles.end(); F++)
	{
		LogAll.info("Processing tree %s",(*F).c_str());
		if (not AnaObj.InitTree(*F))
		{
			LogAll.warn("Skip this run");
			continue;
		}
		AnaObj.LoopTree();
		AnaObj.CloseTree();
	}

	// This is for the total nthrown and similar values stored in the output root file
	AnaObj.StoreExtraValues();

	AnaObj.SaveHistos();

	// clean up and flush all appenders
	LogAll.shutdown();
	cout<<"Done!"<<endl;
	return 0;
}
