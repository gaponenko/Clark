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
	const char* const short_options = "hl:e:o:";
	// Long options are like : "--help --output=ARG"
	// 0 means no ARG and 1 means need an ARG
	const struct option long_options[] = {
		{ "help",			0, NULL, 'h'	},
		{ "ecalib",			1, NULL, 'e'	},
		{ "logfile",		1, NULL, 'l'	},
		{ "outputfile",		1, NULL, 'o'	},
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

	// Did the user give a correct command line ?
	if ( argc < 3 )
	{
		cerr<<" *** ERROR *** Please check that you put all the required arguments. See the --help option"<<endl;
		return 1;
	}

	// We grab the argument that are not options
	char const *const *Tmp;
	int Nb_InFiles	= argc - optind;
	Tmp = (char const *const *) argv + optind;

	Executable	= argv[0];

	ConfFile	= Tmp[0];

	for ( int i = 1; i < Nb_InFiles; i++ )
		Trees.push_back(Tmp[i]);
	


	return 0;
}

void	CommandLine::SaveToConfig(ConfigFile &Conf)
{
	Conf.add( "CommandLine/Debug",		Debug);
	Conf.add( "CommandLine/LogFile",	LogFile);
	Conf.add( "CommandLine/OutputFile",	OutputFile);
	Conf.add( "CommandLine/EcalFile",	EcalFile);
	Conf.add( "CommandLine/EcalArray",	EcalArray);
}

void	CommandLine::Print_Help()
{
	cout<<endl<<"Usage:"<<endl
		<<" Clark [-e ECALIBFILE:ECALIBARRAY] [-l LOG] [-o ROOTOUTPUT] CONFIGFILE TREEFILE1 TREEFILE2 ..."<<endl
		<<" Clark [-e ECALIBFILE:ECALIBARRAY] [-l LOG] [-o ROOTOUTPUT] CONFIGFILE TREEDIR1 TREEDIR2 ..."<<endl
		<<endl<<"options:"<<endl
		<<"  -e, --ecalib       define the energy calibration file and the name of the TArray containing the calibration"<<endl
		<<"  -h, --help         print this help"<<endl
		<<"  -l, --logfile      define the log filename"<<endl
		<<"  -o, --outputfile   define the root output filename"<<endl
		<<endl<<"examples:"<<endl
		<<"  To sum a couple of set75a12 trees:"<<endl
		<<"    Clark set75a12s1_clk /twist/tw00y/systematics/data/set75/anal12/root/run44444/tree044444.root /twist/tw00y/systematics/data/set75/anal12/root/run44445/tree044445.root"<<endl
		<<"  To sum all the good runs of set75a12:"<<endl
		<<"    Clark set75a12s1_clk /twist/tw00y/systematics/data/set75/anal12/goodlinks/*"<<endl
		<<endl;

}

