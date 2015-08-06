/*
 * Utils.h
 *
 *  Created on: Jul 13, 2015
 *      Author: joppich
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <inttypes.h>
#include <vector>
#include <stdio.h>

template <typename T1, typename T2>
class Utils {

public:

	static uint32_t find(std::vector< std::pair<T1,T2> >* pVec, T1* pElem1, T2* pElem2)
	{

		if (pVec == NULL)
			return -1;

		if (pElem1 != NULL)
		{

			for (uint32_t i = 0; i < pVec->size(); ++i)
			{
				if (pVec->at(i).first == *pElem1)
					return i;
			}

		} else {

			for (uint32_t i = 0; i < pVec->size(); ++i)
			{
				if (pVec->at(i).second == *pElem2)
					return i;
			}

		}

		return -1;

	}

};

#endif /* UTILS_H_ */
