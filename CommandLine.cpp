#include "CommandLine.h"


CommandLine::CommandLine()
{
	Debug		= false;
	ConfFile	= "";
	Executable	= "";
	EcalFile	= "";
	LogFile		= "";
	OutputFile	= "";
	EcalArray	= "";
}


int CommandLine::Init(int argc, char **argv)
{
	int		next_option;
	string	TmpStr;
	int		nb_option = -1;
	// Short option definition "-h -o ARG"
	const char* const short_options = "hl:e:";
	// Long options are like : "--help --output=ARG"
	// 0 means no ARG and 1 means need an ARG
	const struct option long_options[] = {
		{ "help",			0, NULL, 'h'	},
		{ "ecalib",			1, NULL, 'e'	},
		{ "logfile",		1, NULL, 'l'	},
		{ NULL,				0, NULL, 0		},
	};

	string tmp;
	do
	{
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
		nb_option++;
	// Check of the return value
		switch (next_option)
		{
		case 'h':
			Print_Help();
			return 1;
		case 'd':
			Debug	= true;
			break;
		case 'e':
			tmp	= optarg;
			EcalFile	= tmp.substr( 0, tmp.find(":", 0));
			EcalArray	= tmp.substr( tmp.find(":", 0)+1);
			break;
		case 'l':
			LogFile		= optarg;
			break;
		case 'o':
			OutputFile	= optarg;
			break;
		case '?':
			cout<<" WARNING : Option -"<<next_option<< " invalid and ignored"<<endl;
		case -1:
			break;
		default:
			break;
		}
	}
	while (next_option != -1);
	// We grab the argument that are not options
	char const *const *Tmp;
	int Nb_InFiles	= argc - optind;
	Tmp = (char const *const *) argv + optind;

	Executable	= argv[0];

	ConfFile	= Tmp[0];

	for ( int i = 1; i < Nb_InFiles; i++ )
		Trees.push_back(Tmp[i]);
	

	// Did the user give a correct command line ?
	// if ( ! IsOk() )
	// {
	// 	cerr<<" *** ERROR *** Please check that you put all the required arguments. See the --help option"<<endl;
	// 	return 1;
	// }

	return 0;
}

void	CommandLine::SaveToConfig(ConfigFile &Conf)
{
	Conf.add( "CommandLine/Debug",		Debug);
	Conf.add( "CommandLine/LogFile",	LogFile);
	Conf.add( "CommandLine/EcalFile",	EcalFile);
	Conf.add( "CommandLine/EcalArray",	EcalArray);
}

void	CommandLine::Print_Help()
{
	cout<<"Usage:"<<endl
		<<"options:"<<endl
		<<endl<<"examples:"<<endl
		<<endl;

}

