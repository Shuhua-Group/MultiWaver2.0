/*
 * exception.h
 *
 *  Created on: Jan 26, 2018
 *      Author: yuankai
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

# include <exception>
# include <iostream>
# include <sstream>

class except
{
public:
	except();
	template <typename T>
	except& operator<<(const T & in)
	{
		std::ostringstream str;
		str << in;
		errInfo += str.str();
		return *this;
	}

	void err();
private:
	std::string errInfo;
};



#endif /* EXCEPTION_H_ */
