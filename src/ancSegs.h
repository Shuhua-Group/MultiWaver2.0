/*
 * ancSegs.h
 *
 *  Created on: Jan 26, 2018
 *      Author: yuankai
 */

#ifndef ANCSEGS_H_
#define ANCSEGS_H_

# include <iostream>
# include <vector>

# include "bootstrap.h"

class ancSegs
{
public:
	ancSegs(const std::string & path, double cutoff);

	friend void ancBootstrap::bootstrap(const ancSegs & dat);

	std::vector<std::vector<double> > segs;
	std::vector<double> ancProp, tancProp;
	std::vector<std::string> popLabels;

	int totNumSegs;
	double minLen, totLen;

private:
	std::vector<double> segs4bootstrap;
	std::vector<int> popLabels4bootstrap;
};

#endif /* ANCSEGS_H_ */
