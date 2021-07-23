/*
 * exception.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: yuankai
 */

# include "exception.h"

# include <sys/utsname.h>
# include <sys/types.h>
# include <unistd.h>
# include <time.h>
# include <cstdlib>
# include <sstream>
# include <iostream>

except::except()
{
    struct utsname name;
    uname(&name);
    std::string SvrName( name.nodename );
    for(unsigned int i = 0 ; i < SvrName.size() ; ++i)
        if(SvrName[i] == ' ')
                SvrName[i] = '_';
    long pid = getpid();

    std::ostringstream os;
    os << "SvrName:\t" << SvrName << std::endl;
    os << "PID:\t" << pid << std::endl;

    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );

    os << "Time:\t" << tmp << std::endl;

    errInfo = os.str();
}

void except::err()
{
	std::cerr << "\033[1m\033[31m";
	std::cout << errInfo << std::endl;
	std::cerr << "\033[0m";
	exit(1);
}


