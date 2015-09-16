/*
 * Split.h
 *
 *  Created on: Jul 7, 2015
 *      Author: joppich
 */

#ifndef SPLITEVENT_H_
#define SPLITEVENT_H_

#include <inttypes.h>
#include <gtxloader/GffEntry.h>
#include <htslib/sam.h>
#include <vector>
#include <sstream>


class SplitEvent : public GenomicRegion {
public:

	struct Entry {
		std::string* pReadName;
		uint32_t iLeftCover;
		uint32_t iRightCover;

		~Entry()
		{
			delete pReadName;
		}
	};

	SplitEvent( uint32_t iChromStart, uint32_t iChromEnd, uint32_t iChromID)
	: GenomicRegion(iChromStart, iChromEnd)
	{

		m_iChromID = iChromID;

		m_pEntries = new std::vector<Entry*>();
	}

	~SplitEvent()
	{

		for (uint32_t i = 0; i < m_pEntries->size(); ++i)
			delete m_pEntries->at(i);

		delete m_pEntries;
	}


	void addEntry( bam1_t* pRead, uint32_t iLeft, uint32_t iRight )
	{
		Entry* pNewEntry = new Entry();
		pNewEntry->pReadName = new std::string( bam_get_qname(pRead) );
                
                m_pEntries->push_back(pNewEntry);
	}

	uint32_t getChromID()
	{
		return m_iChromID;
	}

	uint32_t getEntryCount()
	{
		return m_pEntries->size();
	}

	std::vector<Entry*>* getEntries()
	{
		return m_pEntries;
	}
        
        std::string toString()
        {
            std::ostringstream oStringStream;
            
            oStringStream << m_iChromID << " " << m_iStart << " " << m_iEnd << " " << m_pEntries->size();
            
            return oStringStream.str();
        }

private:

	std::vector<Entry*>* m_pEntries;

	uint32_t m_iChromID;

};

#endif /* SPLITEVENT_H_ */
