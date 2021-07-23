/*
 * MultiWaveInfer2.cpp
 *
 *  Created on: Jan 29, 2018
 *      Author: yuankai
 */

# include "models.h"
# include "argparse.h"
# include "exception.h"
# include "ancSegs.h"

# include <fstream>

# ifndef NOOMP
# include <omp.h>
# endif//NOOMP

using namespace std;

int main(int argc, char **argv)
{
	try
	{
		argparse(argc, argv);
		ofstream fpo(SoftPar.outpath.c_str());
		ofstream fplog;
		if(!fpo)
		{
			except e;
			e << "Error: Cannot open the output file!\n";
			e << "    The input path is \"" << SoftPar.outpath << '\"';
			throw e;
		}
		if(SoftPar.logpath != "")
		{
			fplog.open(SoftPar.logpath.c_str());
			if(!fplog)
			{
				except e;
				e << "Error: Cannot open the log file!\n";
				e << "    The input path is \"" << SoftPar.logpath << '\"';
				throw e;
			}
		}
		ancSegs dat(SoftPar.inpath, SoftPar.minLen);
		models m(dat);
		m.solve();
		std::vector<models*> bootstrap;
		std::vector<ancBootstrap> bot;
		bootstrap.resize(SoftPar.nbootstrap, NULL);
		bot.resize(SoftPar.nbootstrap);
		srand(time(0));
		for(int i = 0 ; i < SoftPar.nbootstrap; ++i)
		{
			bot[i].bootstrap(dat);
			bootstrap[i] = new models(bot[i]);
		}
#ifndef NOOMP
		omp_set_num_threads(SoftPar.nthread);

#pragma omp parallel for
#endif//NOOMP
		for(int i = 0 ; i < SoftPar.nbootstrap; ++i)
			bootstrap[i]->solve();
		if(fplog)
		{
			fplog << "Results:\t" << m.info() << std::endl;
			for(int i = 0 ; i < SoftPar.nbootstrap; ++i)
				fplog << "Bootstrapping " << i + 1 << "\t" << \
						bootstrap[i]->info() << std::endl;
		}
		Summary(m, bootstrap, dat.popLabels, SoftPar.ci, fpo, fplog);
	}
	catch(except & e)
	{
		e.err();
	}
}


