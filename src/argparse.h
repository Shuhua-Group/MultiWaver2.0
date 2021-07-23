/*
 * argparse.h
 *
 *  Created on: Jan 28, 2018
 *      Author: yuankai
 */

#ifndef ARGPARSE_H_
#define ARGPARSE_H_

# include <iostream>

void print_help();

void argparse(int argc, char **argv);

class par
{
public:
	par() : inpath(""), outpath(""), logpath(""), \
			minLen(0), ci(0.95), nbootstrap(0), nthread(1){};
	std::string inpath, outpath, logpath;
	double minLen, ci;
	int nbootstrap, nthread;
};

extern par SoftPar;

#endif /* ARGPARSE_H_ */
