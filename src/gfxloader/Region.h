/*
 * Region.h
 *
 *  Created on: Jul 30, 2015
 *      Author: joppich
 */

#ifndef REGION_H_
#define REGION_H_

#include <inttypes.h>
#include <vector>

class Region {
public:
	Region(uint32_t iStart, uint32_t iEnd)
	{
		m_iStart = iStart;
		m_iEnd = iEnd;
	}

	std::vector<Region*>* split(uint32_t iPosition)
	{

		if ((iPosition <= m_iStart) || (iPosition >= m_iEnd))
			return 0;

		std::vector<Region*>* pNewRegions = new std::vector<Region*>();

		Region* pRegion1 = new Region(m_iStart, iPosition);
		Region* pRegion2 = new Region(iPosition+1, m_iEnd);

		pNewRegions->push_back(pRegion1);
		pNewRegions->push_back(pRegion2);

		return pNewRegions;
	}

private:

	uint32_t m_iStart, m_iEnd;

};

#endif /* REGION_H_ */
