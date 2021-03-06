#ifndef GENOMICREGION_H_
#define GENOMICREGION_H_

#include <inttypes.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdlib.h>


class GenomicRegion
{
public:

	GenomicRegion(uint32_t iStart, uint32_t iEnd)
	{
		m_iStart = iStart;
		m_iEnd = iEnd;
	}

	uint32_t getStart()
	{
		return this->m_iStart;
	}

	uint32_t getEnd()
	{
		return this->m_iEnd;
	}

	uint32_t getLength()
	{

		return this->m_iEnd - this->m_iStart + 1;

	}

	void setStart(int iStart)
	{
		this->m_iStart = iStart;
	}
	void setEnd(int iEnd)
	{
		this->m_iEnd = iEnd;
	}

	bool contains(uint32_t iPosition)
	{
		return ((iPosition >= m_iStart) && (iPosition <= m_iEnd));
	}

	bool contains(GenomicRegion& oOther)
	{
            
		bool bStartContained = this->contains(oOther.m_iStart);
		bool bEndContained = this->contains(oOther.m_iEnd);

		return bStartContained && bEndContained;
	}

	bool contains(GenomicRegion* pOther)
	{

		bool bStartContained = this->contains(pOther->m_iStart);
		bool bEndContained = this->contains(pOther->m_iEnd);

		return bStartContained && bEndContained;
	}
        
        std::string toString()
        {
            
            std::stringstream oStringStream;
            oStringStream << "[ " << this->getStart() << " " << this->getEnd() << " ]" << this->getLength() << std::endl;
            
            return oStringStream.str();
            
        }

protected:

	uint32_t m_iStart, m_iEnd;


};


#endif /*GENOMICREGION_H_*/
