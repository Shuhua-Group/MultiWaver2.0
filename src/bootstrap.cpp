/*
 * bootstrap.cpp
 *
 *  Created on: Jan 27, 2018
 *      Author: yuankai
 */

# include "ancSegs.h"

# include <cstdlib>

void ancBootstrap::bootstrap(const ancSegs & dat)
{
	totNumSegs = dat.totNumSegs;
	minLen = dat.minLen;
	totLen = 0;
	segs.resize(dat.segs.size());
	for(std::vector<std::vector<double> >::iterator it = segs.begin(); \
			it != segs.end(); ++it)
		it->reserve(dat.totNumSegs);
	ancProp.resize(dat.ancProp.size(), 0);
	tancProp.resize(dat.tancProp.size(), 0);
	popLabels.insert(popLabels.begin(), dat.popLabels.begin(), dat.popLabels.end());

	double ttotLen(0);
	for(int i = 0 ; i < dat.totNumSegs; ++i)
	{
		int index = (double) rand() / RAND_MAX * dat.totNumSegs;
		double tseg = dat.segs4bootstrap[index];
		int lab = dat.popLabels4bootstrap[index];
		if(tseg > dat.minLen)
		{
			segs[lab].push_back(tseg);
			ancProp[lab] += tseg;
			totLen += tseg;
		}
		tancProp[lab] += tseg;
		ttotLen += tseg;
	}

	for(unsigned int i = 0 ; i < ancProp.size(); ++i)
	{
		ancProp[i] /= totLen;
		tancProp[i] /= ttotLen;
	}
}

void ancBootstrap::clean()
{
	segs.clear();
	std::vector<std::vector<double> > ssd;
	segs = ssd;
}

