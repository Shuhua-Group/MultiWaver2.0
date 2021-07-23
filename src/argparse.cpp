/*
 * argparse.cpp
 *
 *  Created on: Jan 28, 2018
 *      Author: yuankai
 */

# include "argparse.h"
# include "models.h"
# include "exception.h"

# include <getopt.h>
# include <cstdlib>
# include <unistd.h>
# include <boost/math/distributions/chi_squared.hpp>

par SoftPar;

double cv_chisq(int df, double alpha)
{
	boost::math::chi_squared dist(df);
	return boost::math::quantile(dist, 1 - alpha);
}

void print_help()
{
	std::cout << "MultiWaver 2.0" << std::endl;
	std::cout << "MultiWaver series software designed to infer population admixture history"\
			" in case of various and complex scenarios. The earlier version of MultiWaver "\
			"considered only the discrete admixture models. In the newly developed version,"\
			" we implemented a more flexible framework to automatically select an optimal "\
			"admixture model among discrete models and continuous models." << std::endl;
	std::cout << " -i/--input     <string>    Input of the ancestral tracks [required]" \
			<< std::endl;
	std::cout << " -o/--output    <string>    Output file [required]" << std::endl;
	std::cout << " -g/--log       <string>    Log file of software [optional, no default]" \
			<< std::endl;
	std::cout << " -l/--lower     [double]    Lower bound to discard short tracks "\
			<< "[optional, default 0]" << std::endl;
	std::cout << " -a/--alpha     [double]    Significance level to reject null " \
			<< "hypothesis in LRT (Discrete model) [optional, default 0.001]" << std::endl;
	std::cout << " -e/--epsilon   [double]    Epsilon to check whether a parameter" \
			<< " converge or not (Discrete model) [optional, default 1.0e-6]" << std::endl;
	std::cout << " -p/--minProp   [double]    Minimum survival proportion for a wave" \
			<< " at the final generation (Discrete model) [optional, default 0.05]" \
			<< std::endl;
	std::cout << " -m/--maxIter   [integer]   Maximum number of iterations to scan for " \
			<< "waves of admixture events [optional, default 10000]" << std::endl;
	std::cout << " -w/--waves     [integer]   Maximum number of discrete model waves "\
			"[optional, default 2]" << std::endl;
	std::cout << " -c/--ci        [double]    Value for the bootstrapping "\
			<< "confident interval [optional, default 0.95]" << std::endl;
	std::cout << " -b/--bootstrap [double]    Number of bootstrapping [optional, default 0]" \
			<< std::endl;
	std::cout << " -t/--thread    [integer]   Number of thread [optional, default 1]" \
			<< std::endl;
	std::cout << " -M/--mode      [DMode/CMode/MixMode]" << std::endl;
	std::cout << "                            Model parameter [optional, default MixMode]" \
			<< std::endl;
	std::cout << "                            DMode    Discrete model only. " \
			<< "(Same as MultiWaver-1.0)" << std::endl;
	std::cout << "                            CMode    Continuous model only. " \
			<< "(Same as AdmixInfer-1.0)" << std::endl;
	std::cout << "                            MixMode  Consider all models together." \
			<< std::endl;
	std::cout << " -s/--simple    Run in simple mode (Discrete model)" << std::endl;
	std::cout << " -h/--help      Print this help" << std::endl;
	exit(1);
}

void argparse(int argc, char **argv)
{
	if(argc == 1)
		print_help();
	static struct option long_options[] =
	{
			{"input",	required_argument,	0,	'i'},
			{"output",	required_argument,	0,	'o'},
			{"lower",	required_argument,	0,	'l'},
			{"alpha",	required_argument,	0,	'a'},
			{"epsilon",	required_argument,	0,	'e'},
			{"minProp",	required_argument,	0,	'p'},
			{"maxIter",	required_argument,	0,	'm'},
			{"waves",	required_argument,	0,	'w'},
			{"simple",	no_argument,		0,	's'},
			{"help",	no_argument,		0,	'h'},
			{"ci",		required_argument,	0,	'c'},
			{"bootstrap",	required_argument,	0,	'b'},
			{"mode",	required_argument,	0,	'M'},
			{"thread",	required_argument,	0,	't'},
			{"log",		required_argument,	0,	'g'},
			{0, 0, 0, 0}
	};

	int ch;
	int option_index(0);

	while((ch = getopt_long(argc, argv,"i:o:l:a:e:p:m:w:shc:b:M:t:g:", long_options, &option_index)) != -1)
	{
		switch(ch)
		{
			case 'i' :
			{
				SoftPar.inpath = optarg;
				break;
			}
			case 'o' :
			{
				SoftPar.outpath = optarg;
				break;
			}
			case 'l' :
			{
				SoftPar.minLen = atof(optarg);
				if(SoftPar.minLen < 0)
				{
					except e;
					e << "Error: The minimum ancestral segment length should be positive!\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'a' :
			{
				models::critical = atof(optarg);
				if(models::critical > 1 || models::critical < 0)
				{
					except e;
					e << "Error: Bad value for the significance in LRT!\n";
					e << "    The value should be ∈ (0,1)\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'e' :
			{
				models::epsilon = atof(optarg);
				if(models::epsilon < 0)
				{
					except e;
					e << "Error: The EM converge parameter should be positive!\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'p' :
			{
				models::minProp = atof(optarg);
				if(models::minProp > 1 || models::minProp < 0)
				{
					except e;
					e << "Error: Bad value for the minimum survival proportion \
							for a wave!\n";
					e << "    The value should be ∈ (0,1)\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'm' :
			{
				models::MaxInter = atoi(optarg);
				if(models::MaxInter < 100)
				{
					except e;
					e << "Error: The number of EM iterator should be larger than 100!\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'w' :
			{
				models::MultiMaxWave = atoi(optarg);
				if(models::MultiMaxWave <= 0)
				{
					except e;
					e << "Error: The maximum number of waves should be larger than 0\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 's' :
			{
				models::MultiMaxWave = 1;
				break;
			}
			case 'h' :
			{
				print_help();
				break;
			}
			case 'c' :
			{
				SoftPar.ci = atof(optarg);
				if(SoftPar.ci < 0 || SoftPar.ci > 1)
				{
					except e;
					e << "Error: Bad value for the bootstraping confident interval!\n";
					e << "    The value should be ∈ (0,1)\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'b' :
			{
				SoftPar.nbootstrap = atoi(optarg);
				if(SoftPar.nbootstrap < 0)
				{
					except e;
					e << "Error: The number of bootstraping should be positive!\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'M' :
			{
				std::string mode = optarg;
				if(mode == "CMode")
				{
					models::isHI = true;
					models::isGA = true;
					models::isCGF = true;
					models::isMulti = false;
				}
				else if(mode == "DMode")
				{
					models::isHI = false;
					models::isGA = false;
					models::isCGF = false;
					models::isMulti = true;
				}
				else if(mode == "MixMode")
				{
					models::isHI = true;
					models::isGA = true;
					models::isCGF = true;
					models::isMulti = true;
				}
				else
				{
					std::cerr << "Cannot identify mode option " << mode << std::endl;
					std::cerr << "Set as defuatl: MixMode" << std::endl;
					models::isHI = true;
					models::isGA = true;
					models::isCGF = true;
					models::isMulti = true;
				}
				break;
			}
			case 't' :
			{
				SoftPar.nthread = atoi(optarg);
				if(SoftPar.nthread <= 0)
				{
					except e;
					e << "Error: The number of threads should be larger than 0\n";
					e << "    The input is \"" << optarg << '\"';
					throw e;
				}
				break;
			}
			case 'g' :
			{
				SoftPar.logpath = optarg;
				break;
			}
			default:
			{
				print_help();
				break;
			}
		}
	}
	models::critical = cv_chisq(2, models::critical);
}
