/*
 * bootstrap.h
 *
 *  Created on: Jan 27, 2018
 *      Author: yuankai
 */

#ifndef BOOTSTRAP_H_
#define BOOTSTRAP_H_

class ancSegs;

class ancBootstrap
{
public:
	void bootstrap(const ancSegs & dat);
	void clean();

	std::vector<std::vector<double> > segs;
	std::vector<double> ancProp, tancProp;
	std::vector<std::string> popLabels;

	int totNumSegs;
	double minLen, totLen;
};



#endif /* BOOTSTRAP_H_ */
