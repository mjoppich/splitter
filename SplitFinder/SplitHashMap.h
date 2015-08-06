/*
 * SplitHashMap.h
 *
 *  Created on: Jul 20, 2015
 *      Author: joppich
 */

#ifndef SPLITHASHMAP_H_
#define SPLITHASHMAP_H_

#include "../Utils/ThreadedHashMap.h"

class Split;

class SplitHashMap : public ThreadedHashMap<uint32_t, Split*>{
public:
	SplitHashMap(uint32_t iBins)
		: 	ThreadedHashMap(iBins)
	{

	}

	std::vector<Split*> *findElements(std::vector<Split*>* pLookIn, uint32_t oKey)
	{
		return NULL;
	}

	uint32_t getIndexForKey(uint32_t okey)
	{
		return 0;
	}

	uint32_t getKeyForElement(Split* pSplit)
	{
		return 0;
	}

	Split* add(Split* pSplit)
	{

		uint32_t iKey = this->getKeyForElement(pSplit);

		this->insert(iKey, pSplit);

		return pSplit;

	}

};

#endif /* SPLITHASHMAP_H_ */
